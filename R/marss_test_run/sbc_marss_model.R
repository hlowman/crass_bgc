# MARSS model @ Santa Barbara Coastal LTER Sites
# October 24, 2021
# Heili Lowman

# The following script will run a MARSS analysis at SBC LTER sites for the CRASS project.
# M - Multi-variate
# AR - Auto-regressive
# S - State
# S - Space

# Note: As of November 3, I've been encountering issues with missing data so this first go-round, trying to make it work with just a few sites at which I've verified the precip and fire covariate data is complete.

#### Setup ####

# Load packages
library(tidyverse)
library(lubridate)
library(here)
library(MARSS)

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))

# Load datasets
# Stream Chemistry - all sites
chem <- readRDS("data_working/SBchem_edited_110721.rds")
# Precipitation - all sites
precip <- readRDS("data_working/SBprecip_edited_110721.rds")
# Fire Events - HO00 and RG01 sites only so far
fire <- readRDS("data_working/SBfire_edited_111021.rds")
# Site Location information
location <- read_csv("data_raw/sbc_sites_stream_hydro.csv")

#### Filter & Joining Data ####

# First, to avoid strange gaps in data, I'm creating a dummy column with the dates of interest.
# And I need to do this for each site, since otherwise it'll allow gaps.
# Just being really careful since this was causing issues previously.
dates1 <- data.frame(seq(as.Date("2002/9/1"), by = "month", length.out = 166)) %>%
  dplyr::rename(Date = 'seq.as.Date..2002.9.1....by....month...length.out...166.') %>%
  mutate(site = "HO00")

dates2 <- data.frame(seq(as.Date("2002/9/1"), by = "month", length.out = 166)) %>%
  dplyr::rename(Date = 'seq.as.Date..2002.9.1....by....month...length.out...166.') %>%
  mutate(site = "RG01")
## by month from 9/1/2002 to 7/1/2016

dates <- rbind(dates1, dates2)

# And to join the datasets, I need to identify the matching stream sites for the precip data
precip_ed <- precip %>%
  mutate(sitecode_match = factor(case_when(
    sitecode == "GV202" ~ "GV01",
    sitecode == "HO201" ~ "HO00",
    sitecode == "RG202" ~ "RG01",
    sitecode == "TECA" ~ "TO02",
    sitecode == "GLAN" ~ "BC02",
    sitecode == "SMPA" ~ "SP02",
    sitecode == "GOFS" ~ "DV01",
    sitecode == "GORY" ~ "AT07",
    sitecode == "CAWTP" ~ "AB00",
    sitecode == "SBEB" ~ "MC00",
    sitecode == "ELDE" ~ "RS02"))) %>%
  dplyr::rename(site_precip = site,
         sitecode_precip = sitecode) %>%
  mutate(Day = 1) %>% # new column of "days"
  mutate(Date = make_date(Year, Month, Day))

# First, join precip with the dates.
dates_precip <- left_join(dates, precip_ed, by = c("Date", "site" = "sitecode_match"))
# Then, left join precip with chemistry so as not to lose any data.
precip_chem <- left_join(dates_precip, chem, by = c("Year", "Month", "site"))
# And again left join with fire so as not to accidentally drop data.
dat <- left_join(precip_chem, fire, by = c("Date", "site"))

# Adding in dummy covariates by season
n_months <- dat$Date %>%
  unique() %>%
  length()
seas_1 <- sin(2 * pi * seq(n_months) / 12)
seas_2 <- cos(2 * pi * seq(n_months) / 12)

dat <- dat %>%
  mutate(Season1 = rep(seas_1, 2),
         Season2 = rep(seas_2, 2))

# And, for this first attempt to try and get the MARSS model working, I'm going to filter down to only HO00 and RG01 sites.
dat_2 <- dat %>%
  filter(site %in% c("HO00", "RG01"))

# AJW: replace NaNs with NAs
dat_2[is.nan(dat_2)] = NA

# And, inspect dataset for missing covariate data.
sum(is.na(dat_2$fire)) # 0
sum(is.na(dat_2$cumulative_precip_mm)) # 0
# Great!

# Note for future me - be VERY careful with the joining above. Something weird was happening previously where precip data that IS present was simply dropping off.

#### Model fit ####

# Data: Stream Chemistry analytes (NH4, NO3, TDN, TPN, PO4, TDP, TP, TPP, TPC, TSS, SCond)
# Covariates: Month, Precip, Fire

# Starting with NH4 for test run.

# Note: Not scaling for now, but this should also be added in later.
dat_nh4 <- dat_2 %>%
  select(site, Year, Month, Season1, Season2, mean_nh4_uM, cumulative_precip_mm, HO00_Gaviota, HO00_Sherpa, HO00_Whittier, RG01_Gaviota, RG01_Sherpa, RG01_Whittier) %>% # keeping year and month in so I don't get the "duplicates arise" error
  pivot_wider(names_from = site, values_from = c(mean_nh4_uM, cumulative_precip_mm, HO00_Gaviota, HO00_Sherpa, HO00_Whittier, RG01_Gaviota, RG01_Sherpa, RG01_Whittier)) %>%
  select(-c(HO00_Gaviota_RG01, HO00_Sherpa_RG01, HO00_Whittier_RG01, RG01_Gaviota_HO00, RG01_Sherpa_HO00, RG01_Whittier_HO00))

par(mfrow=c(2,1))
plot(dat_nh4$mean_nh4_uM_HO00, type="b")
plot(dat_nh4$mean_nh4_uM_RG01, type="b")
par(mfrow=c(1,1))
  
# Pull out only NH4 data
dat_dep <- t(dat_nh4[,5:6])

# Make covariate inputs
# dat_cov <- dat_nh4[,c(3:4,7:14)]
# dat_cov <- t(scale(dat_cov))
# removed Whittier fire because it is all zeros and identical covars messes MARSS up!
dat_cov <- dat_nh4[,c(3:4,7:10,12:13)]
dat_cov <- t(scale(dat_cov))

# Replace NaNs with 0, which don't play nice with the scale()
# dat_cov[is.nan(dat_cov)] <- 0

# make C matrix
# CC <- matrix(list("Season1", "Season1", "Season2", "Season2", "HO00_precip", 0, 0, "RG01_precip", "HO00_Gaviota_fire", 0, "HO00_Sherpa_fire", 0,"HO00_Whittier_fire", 0, 0, "RG01_Gaviota_fire", 0, "RG01_Sherpa_fire", 0,"RG01_Whittier_fire"),2,10)
# this matrix controls what covars predict what response vars; in contrast to...
# unconstrained: seperately estimates ALL correlations of predictor vars with x, i.e. how predictor vars drive ts dynamics. I think this structure is ok given that covars are unique to each site, but it would not be ok if there is a mix of shared and unique covars (are year and month used? bc if so, these are shared)
# removed Whittier fire because it is all zeros and identical covars messes MARSS up!
CC <- matrix(list("Season1", "Season1", "Season2", "Season2", "HO00_precip", 0, 0, "RG01_precip", "HO00_Gaviota_fire", 0, "HO00_Sherpa_fire", 0, 0, "RG01_Gaviota_fire", 0, "RG01_Sherpa_fire"),2,8)


# Model setup
mod_list <- list(
  B = "diagonal and unequal",
  U = "zero", # zero: does NOT allow a drift term in process model to be estimated # removing due to lack of anticipated monotonic trend
  C = CC, # see Alex's matrix above
  c = dat_cov, # we should probably de-mean and scale covariates to units of sd 
  Q = "diagonal and unequal", # diagonal and unequal: allows for and estimates the covariance matrix of process errors
  Z = "identity", # identity: estimated state processes are unique to each site
  A = "zero",
  R = "diagonal and equal" # diagonal and equal: allows for and estimates the covariance matrix of observations errors (may want to provide a number for this from method precision etc if possible)
)

# Fit model
# fit <- MARSS(y = dat_dep, model = mod_list,
#                 control = list(maxit = 5000), method = "BFGS")
# see pg 5 in MARSS manual for notes on method BFGS vs method EM: EM algorithm gives more robust estimation for datasets replete with missing values and for high-dimensional models with various constraints. BFGS is faster and is good enough for some datasets. Typically, both should be tried.

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE)

#### Plotting Results ####

## check for hidden errors
fit[["errors"]]
# no hidden errors

### Plot coef and coef estimates ###

## estimates

# non-parametric:
#est = MARSSparamCIs(fit, method = "innovations", alpha = 0.05, nboot = 70, silent=F) 
# ^ Errors were caught in MARSSboot :Innovations bootstrapping uses the innovations resampling and can only be done if there are no missing values in the data.
# Error: Stopped in MARSSboot() due to problem(s) with function arguments.

# parametric:
# est_fit = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F) # nboot should be ~ 2000 for final results
# note - this code takes a while to run, so can import "CIs_fit.csv" from the
# data_working/marss_test_run folder to bypass this and the next few steps

# hessian method is much fast but not ideal for final results
est_fit = MARSSparamCIs(fit)

CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)

## save CI table ##
# write_csv(CIs_fit, "data_working/marss_test_run/CIs_fit.csv")

### PLOT HO00 ###
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])
CIs_HO00$parm_name = c("Season1","Season2","Precip. (monthly cumm.)", "Gaviota Fire", "Sherpa Fire")
# plot
RESULTS_HO00 = 
  ggplot(CIs_HO00, aes(parm_name, Est.)) + 
  geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
  geom_point(position=position_dodge(width=0.3), size=2) + 
  theme_bw()+
  theme(plot.title = element_text(size = 8)) +
  theme(axis.text = element_text(size = 8)) +
  geom_hline(aes(yintercept=0), linetype="dashed")+
  coord_flip()+ ggtitle("HO00 NH4+")+ 
  ylab("") +
  theme(plot.margin=unit(c(.2,-.2,.05,.01),"cm"))
RESULTS_HO00

### PLOT RG01 ###
CIs_RG01 = rbind(CIs_fit[1:2,], CIs_fit[grepl("RG01", CIs_fit$parm),])
CIs_RG01$parm_name = c("Season1","Season2","Precip. (monthly cumm.)", "Gaviota Fire", "Sherpa Fire")
# plot
RESULTS_RG01 = 
  ggplot(CIs_RG01, aes(parm_name, Est.)) + 
  geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
  geom_point(position=position_dodge(width=0.3), size=2) + 
  theme_bw()+
  theme(plot.title = element_text(size = 8)) +
  theme(axis.text = element_text(size = 8)) +
  geom_hline(aes(yintercept=0), linetype="dashed")+
  coord_flip()+ ggtitle("RG01 NH4+")+ 
  ylab("") +
  theme(plot.margin=unit(c(.2,-.2,.05,.01),"cm"))
RESULTS_RG01

gridExtra::grid.arrange(RESULTS_HO00, RESULTS_RG01, nrow=1)
g <- gridExtra::arrangeGrob(RESULTS_HO00, RESULTS_RG01, nrow=1)
ggsave("figures/MARSS_NH4_HO00_RG01.pdf", g)


#### Script for diagnoses ####

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- residuals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity", 
  A = "zero",
  R = "diagonal and equal" 
)
mod.null <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE)

bbmle::AICtab(fit, mod.null)

#           dAIC df
# fit       0   15
# mod.null 11   7 
# RESULT: covariate model is better than null model - suggests covars add information

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
forecast::Acf(resids$model.residuals[1,], main="HO00 model residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$state.residuals[1,], main="HO00 state residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$model.residuals[2,], main="RG01 model residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$state.residuals[2,], main="RG01 state residuals", na.action=na.pass, lag.max = 24)
mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
# RESULT: yes, at lag 1 for RG01 and in multiple locations for HO00

### Do resids have temporal trend? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
plot(resids$model.residuals[1,], ylab="model residual", xlab="", main="HO00 model residuals")
abline(h=0)
plot(resids$state.residuals[1,], ylab="state residual", xlab="", main="HO00 state residuals")
abline(h=0)
plot(resids$model.residuals[2,], ylab="model residual", xlab="", main="RC01 model residuals")
abline(h=0)
plot(resids$state.residuals[2,], ylab="state residual", xlab="", main="RC01 state residuals")
abline(h=0)
mtext("Do resids have temporal trend?", outer = TRUE, cex = 1.5)

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
qqnorm(resids$model.residuals[1,], main="HO00 model residuals", pch=16, 
       xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[1,])[1]))
qqline(resids$model.residuals[1,])
qqnorm(resids$state.residuals[1,], main="HO00 state residuals", pch=16, 
       xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[1,])[1]))
qqline(resids$state.residuals[1,])
qqnorm(resids$model.residuals[2,], main="RG01 model residuals", pch=16, 
       xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[2,])[1]))
qqline(resids$model.residuals[2,])
qqnorm(resids$state.residuals[2,], main="RG01 state residuals", pch=16, 
       xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[2,])[1]))
qqline(resids$state.residuals[2,])
mtext("Are resids normal?", outer = TRUE, cex = 1.5)

### residuals vs fitted ###
par(mfrow=c(2,2),oma = c(2, 0, 2, 0))
scatter.smooth(as.vector(kf$xtT[1,]) ~ resids$model.residuals[1,],
               main = "HO00 model resids vs fitted (y conditioned)")
abline(lm(as.vector(kf$xtT[1,]) ~ resids$model.residuals[1,]), col="blue")
scatter.smooth(as.vector(kf$xtT[1,]) ~ resids$state.residuals[1,],
               main = "HO00 state resids vs fitted (y conditioned)")
abline(lm(as.vector(kf$xtT[1,]) ~ resids$state.residuals[1,]), col="blue")
scatter.smooth(as.vector(kf$xtT[2,]) ~ resids$model.residuals[2,],
               main = "RC01 model resids vs fitted (y conditioned)")
abline(lm(as.vector(kf$xtT[2,]) ~ resids$model.residuals[2,]), col="blue")
scatter.smooth(as.vector(kf$xtT[2,]) ~ resids$state.residuals[2,],
               main = "RC01 state resids vs fitted (y conditioned)")
abline(lm(as.vector(kf$xtT[2,]) ~ resids$state.residuals[2,]), col="blue")
mtext("Trends in fitted values vs residuals?", outer = TRUE, cex = 1.5)

### residuals vs observed ###
par(mfrow=c(2,2),oma = c(2, 0, 2, 0))
scatter.smooth(as.vector(dat[1,]) ~ resids$model.residuals[1,],
               main = "HO00 model resids vs observed")
abline(lm(as.vector(kf$xtT[1,]) ~ resids$model.residuals[1,]))
scatter.smooth(as.vector(dat[1,]) ~ resids$state.residuals[1,], 
               main = "HO00 state resids vs observed") 
abline(lm(as.vector(kf$xtT[1,]) ~ resids$state.residuals[1,]))
scatter.smooth(as.vector(dat[2,]) ~ resids$model.residuals[2,],
               main = "RC01 model resids vs observed")
abline(lm(as.vector(kf$xtT[2,]) ~ resids$model.residuals[2,]))
scatter.smooth(as.vector(dat[2,]) ~ resids$state.residuals[2,], 
               main = "RC01 state resids vs observed") 
abline(lm(as.vector(kf$xtT[2,]) ~ resids$state.residuals[2,]))
mtext("Trends in observed values vs residuals?", outer = TRUE, cex = 1.5)

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

# End of script.

