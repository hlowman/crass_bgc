## This script loads the munged `solutes` data frame and fits MARSS models

##-----------------------
## necessary user inputs
##-----------------------

yr_first <- 1984
yr_last <- 2010


## load libraries
library(dplyr)
## library(tidyr)
library(MARSS)

## load solutes data
df <- readr::read_csv(file.path(here::here("data"), "tbl_solutes_unmanaged_mon.csv"))
# df <- readr::read_csv(file.path(here::here("data"), "tbl_solutes_mon.csv"))


## CAVEATS

## 1) remove SLP
df <- df %>%
  filter(site != "SLP" & site != "SEF")

## 2) trim data set to common time frame
df <- df %>%
  filter(dec_water_yr >= yr_first & dec_water_yr <= yr_last)

## 3) convert 0 to NA following mtg and emails (ie, "no sample")
df <- df %>%
   purrr::modify_at(-c(1:4), ~na_if(., 0))

## dummy covariates for season
n_months <- df$dec_water_yr %>%
  unique() %>%
  length()
seas_1 <- sin(2 * pi * seq(n_months) / 12)
seas_2 <- cos(2 * pi * seq(n_months) / 12)


## names of solutes
solutes <- df %>%
  select(starts_with("FWA")) %>%
  colnames

## empty lists for model fits
mod_set_RW <- vector("list", length(solutes))
mod_set_RW_b <- vector("list", length(solutes))
mod_set_one_RW <- vector("list", length(solutes))
mod_set_one_RW_b <- vector("list", length(solutes))
mod_set_region_RW <- vector("list", length(solutes))
mod_set_region_RW_b <- vector("list", length(solutes))
mod_set_site_RW <- vector("list", length(solutes))
mod_set_site_RW_b <- vector("list", length(solutes))

## loop over solutes
for(i in 1:length(solutes)) {
  
  ## assign df to temp file
  df_tmp <- df
  
  ## screen for lack of TDP data
  if(solutes[i] == "FWATDPmgL") {
    df_tmp <- df_tmp %>%
      filter(site != "BBWM" & site != "SLP")
  }
  
  ## select solute and pivot wider
  dat_sol <- df_tmp %>%
    select(region:dec_water_yr, all_of(solutes[i])) %>%
    tidyr::pivot_wider(names_from = c(region, site, catchment),
                       values_from = solutes[i]) %>%
    select(-dec_water_yr) %>%
    log() %>%
    scale(scale = FALSE) %>%
    t()
  
  ## SET 1: all catchments are unique processes
  ## (a) RW without bias
  ## model setup
  mod_list <- list(
    B = "identity",
    U = "zero",
    C = "unconstrained",
    c = rbind(seas_1, seas_2),
    Q = "diagonal and unequal",
    Z = "identity",
    A = "zero",
    R = "diagonal and equal"
  )
  ## fit model
  mod_set_RW[[i]] <- MARSS(dat_sol, model = mod_list,
                            control = list(maxit = 5000), method = "BFGS")

  ## (b) RW with bias
  mod_list$U = "unconstrained"
  ## fit model
  mod_set_RW_b[[i]] <- MARSS(dat_sol, model = mod_list,
                              control = list(maxit = 5000), method = "BFGS")

  ## SET 2: all catchments are shared processes within N Am
  ## (a) RW without bias
  ## model setup
  mod_list <- list(
    B = "identity",
    U = "zero",
    C = "unconstrained",
    c = rbind(seas_1, seas_2),
    Q = "diagonal and unequal",
    Z = matrix(1, nrow = nrow(dat_sol), ncol = 1),
    A = "zero",
    R = "diagonal and equal"
  )
  ## fit model
  mod_set_one_RW[[i]] <- MARSS(dat_sol, model = mod_list,
                               control = list(maxit = 5000), method = "BFGS")

  ## (b) RW with bias
  mod_list$U = "unconstrained"
  ## fit model
  mod_set_one_RW_b[[i]] <- MARSS(dat_sol, model = mod_list,
                                 control = list(maxit = 5000), method = "BFGS")
  
  ## SET 3: all catchments share a process within a region
  ## (a) RW without bias
  ## regions for all catchments
  names_regions <- dat_sol %>%
    rownames %>%
    gsub(pattern = "(.{2})(_.*)", replacement = "\\1")
  ##  unique regions 
  regions <- unique(names_regions)
  ## number of regions
  n_regions <- length(regions)
  ## empty Z matrix
  ZZ <- matrix(0, nrow = nrow(dat_sol), ncol = n_regions)
  colnames(ZZ) <- regions
  ## fill in Z with correct regions
  for(j in 1:n_regions) {
    ZZ[names_regions == regions[j], j] <- 1
  }
  ## model setup
  mod_list <- list(
    B = "identity",
    U = "zero",
    C = "unconstrained",
    c = rbind(seas_1, seas_2),
    Q = "diagonal and unequal",
    Z = ZZ,
    A = "zero",
    R = "diagonal and equal"
  )
  ## fit model
  mod_set_region_RW[[i]] <- MARSS(dat_sol, model = mod_list,
                               control = list(maxit = 5000), method = "BFGS")
  
  ## (b) RW with bias
  ## U for biased RW's
  mod_list$U <- matrix(regions, ncol = 1)
  ## fit model
  mod_set_region_RW_b[[i]] <- MARSS(dat_sol, model = mod_list,
                                 control = list(maxit = 5000), method = "BFGS")
  
  ## SET 4: all catchments share a process within a site
  ## (a) RW without bias
  ## sites for all catchments
  names_sites <- dat_sol %>%
    rownames %>%
    gsub(pattern = "(.{2}_)(.{3,4})(_.*)", replacement = "\\2")
  ##  unique regions 
  sites <- unique(names_sites)
  ## number of regions
  n_sites <- length(sites)
  ## empty Z matrix
  ZZ <- matrix(0, nrow = nrow(dat_sol), ncol = n_sites)
  colnames(ZZ) <- sites
  ## fill in Z with correct regions
  for(j in 1:n_sites) {
    ZZ[names_sites == sites[j], j] <- 1
  }
  ## model setup
  mod_list <- list(
    B = "identity",
    U = "zero",
    C = "unconstrained",
    c = rbind(seas_1, seas_2),
    Q = "diagonal and unequal",
    Z = ZZ,
    A = "zero",
    R = "diagonal and equal"
  )
  ## fit model
  mod_set_site_RW[[i]] <- MARSS(dat_sol, model = mod_list,
                                  control = list(maxit = 5000), method = "BFGS")
  
  # ## (b) RW with bias
  ## U for biased RW's
  mod_list$U <- matrix(sites, ncol = 1)
  ## fit model
  mod_set_site_RW_b[[i]] <- MARSS(dat_sol, model = mod_list,
                                    control = list(maxit = 5000), method = "BFGS")
  
}

## save results
## unique states
saveRDS(mod_set_RW,
        file = file.path(here::here("analysis"), "fitted_seas_unique_states_RW.rds"))
saveRDS(mod_set_RW_b,
        file = file.path(here::here("analysis"), "fitted_seas_unique_states_RW_b.rds"))
## site states
saveRDS(mod_set_site_RW,
        file = file.path(here::here("analysis"), "fitted_seas_site_state_RW.rds"))
saveRDS(mod_set_site_RW_b,
        file = file.path(here::here("analysis"), "fitted_seas_site_state_RW_b.rds"))
## regional states
saveRDS(mod_set_region_RW,
        file = file.path(here::here("analysis"), "fitted_seas_region_state_RW.rds"))
saveRDS(mod_set_region_RW_b,
        file = file.path(here::here("analysis"), "fitted_seas_region_state_RW_b.rds"))
## global N Am state
saveRDS(mod_set_one_RW,
        file = file.path(here::here("analysis"), "fitted_seas_one_state_RW.rds"))
saveRDS(mod_set_one_RW_b,
        file = file.path(here::here("analysis"), "fitted_seas_one_state_RW_b.rds"))

## model selection table
tbl_mod_aic <- data.frame(model = c("unique states", "unique states + bias",
                                    "site states", "site states + bias",
                                    "regional states", "regional states + bias",
                                    "one N Am state", "one N Am state + bias"),
                          Ca = rep(NA, 8),
                          DOC = rep(NA, 8),
                          NH4 = rep(NA, 8),
                          NO3 = rep(NA, 8),
                          TDP = rep(NA, 8),
                          SO4 = rep(NA, 8))

tbl_mod_aic[1,-1] <- sapply(mod_set_RW, AIC)
tbl_mod_aic[2,-1] <- sapply(mod_set_RW_b, AIC)
tbl_mod_aic[3,-1] <- sapply(mod_set_site_RW, AIC)
tbl_mod_aic[4,-1] <- sapply(mod_set_site_RW_b, AIC)
tbl_mod_aic[5,-1] <- sapply(mod_set_region_RW, AIC)
tbl_mod_aic[6,-1] <- sapply(mod_set_region_RW_b, AIC)
tbl_mod_aic[7,-1] <- sapply(mod_set_one_RW, AIC)
tbl_mod_aic[8,-1] <- sapply(mod_set_one_RW_b, AIC)

for(i in 1:length(solutes)) {
  tbl_mod_aic[,i+1] <- round(tbl_mod_aic[,i+1] - min(tbl_mod_aic[,i+1]), 1)
}

tbl_mod_aic

## bootstrap biased RW models by site
bias_bootstrap <- lapply(mod_set_site_RW_b, MARSSparamCIs, nboot = 1000)

## create table of bias estimates (+/- CI) 
tmp <- list()
for(i in 1:length(solutes)) {
  bias_ID <- grep("U.*", names(bias_bootstrap[[i]]$parMean))
  tmp$solute <- rep(sub("(FWA)(.*)(mgL)", "\\2", solutes[i]), length(bias_ID))
  tmp$site <- sub("(U\\.)(.*)", "\\2", names(bias_bootstrap[[i]]$parMean[bias_ID]))
  tmp$bias <- bias_bootstrap[[i]]$parMean[bias_ID]
  tmp$loCI <- bias_bootstrap[[i]]$par.lowCI$U
  tmp$upCI <- bias_bootstrap[[i]]$par.upCI$U
  if(i == 1) {
    tbl_bias_bootstrap <- data.frame(tmp)
  } else {
    tbl_bias_bootstrap <- rbind(tbl_bias_bootstrap, data.frame(tmp))
  }
}
tbl_bias_bootstrap[, -c(1:2)] <- signif(tbl_bias_bootstrap[, -c(1:2)], 3)
rownames(tbl_bias_bootstrap) <- NULL
tbl_bias_bootstrap


