# Nutrient MARSS models
# with fire x ppt interactions and legacy effects
# as well as 4 state structure for CA sites
# and observation error
# Script started August 9, 2023 by Heili Lowman

# This script will run 12 NH4 MARSS models.
# Note, each model fit will remove all stored data, until you reach the IC
# portion of this script.

#### Setup ####

# Load packages.
library(tidyverse)
library(lubridate)
library(MARSS)
library(naniar) 
library(here)
library(bbmle)
library(broom)
library(stats4)

#### Summary Stats ####

# load data with fire x ppt interactions and legacy effects
dat <- readRDS("data_working/marss_data_sb_092123.rds")

# Create summary table.
dat_summary <- dat %>%
  summarize(min_NH4 = min(vwm_nh4, na.rm = TRUE),
            max_NH4 = max(vwm_nh4, na.rm = TRUE),
            mean_NH4 = mean(vwm_nh4, na.rm = TRUE),
            sd_NH4 = sd(vwm_nh4, na.rm = TRUE),
            med_NH4 = median(vwm_nh4, na.rm = TRUE)) %>%
  ungroup()

#### 0y legacy, 4 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_092123.rds")

# pivot wider for MARSS format
# selecting interaction term with no lag
dat_nh4 <- dat %>%
  select(
    site, index, 
    vwm_nh4, 
    cumulative_precip_mm, 
    fire_perc_ws, fire_perc_ws_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_nh4, cumulative_precip_mm, 
                    fire_perc_ws, fire_perc_ws_ppt))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 218
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,
  0,0,"cumulative_precip_mm_HO00",0,
  0,0,0,"cumulative_precip_mm_RS02",
  # fire_perc: fire effect in 2 m window
  "fire_perc_ws_AB00",0,0,0,
  0,"fire_perc_ws_GV01",0,0,
  0,0,"fire_perc_ws_HO00",0,
  0,0,0,"fire_perc_ws_RS02",
  # fire_perc_ws_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_AB00",0,0,0,
  0,"fire_perc_ws_ppt_GV01",0,0,
  0,0,"fire_perc_ws_ppt_HO00",0,
  0,0,0,"fire_perc_ws_ppt_RS02"), 4, 12)

##### Model setup for MARSS 

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", # 4 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # observation error turned on
  ### initial conditions ###
  #x0 = matrix(x0_fixed),
  V0="zero" ,
  tinitx=0
)

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_092123_4state_nh4_mBFGS.rds")

##### Diagnoses 

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Compare to null model ###
# make sure this matches the fitted model in all ways besides the inclusion of C and c
mod_list_null <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  #C = CC, 
  #c = dat_cov,
  Q = "diagonal and unequal", # 4 states
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # observation error turned on
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit, output = c("AICc")) # 1200.804 
MARSSaic(null.fit, output = c("AICc")) # 1307.385

### Autoplot diagnoses ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Looks alright

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Eh, HO00
# doesn't look great, but not much we can do there.

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 0y legacy, 1 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_092123.rds")

# pivot wider for MARSS format
# selecting interaction term with no lag
dat_nh4 <- dat %>%
  select(
    site, index, 
    vwm_nh4, 
    cumulative_precip_mm, 
    fire_perc_ws, fire_perc_ws_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_nh4, cumulative_precip_mm, 
                    fire_perc_ws, fire_perc_ws_ppt))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 218
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,
  0,0,"cumulative_precip_mm_HO00",0,
  0,0,0,"cumulative_precip_mm_RS02",
  # fire_perc: fire effect in 2 m window
  "fire_perc_ws_AB00",0,0,0,
  0,"fire_perc_ws_GV01",0,0,
  0,0,"fire_perc_ws_HO00",0,
  0,0,0,"fire_perc_ws_RS02",
  # fire_perc_ws_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_AB00",0,0,0,
  0,"fire_perc_ws_ppt_GV01",0,0,
  0,0,"fire_perc_ws_ppt_HO00",0,
  0,0,0,"fire_perc_ws_ppt_RS02"), 4, 12)

# Make Q matrix
# see page 28 for Q matrix structure:
# https://cran.r-project.org/web/packages/MARSS/vignettes/UserGuide.pdf
QQ <- matrix(list("s1","b12","b13","b14",
                  "b12","s2","b23","b24",
                  "b13","b23","s3","b34",
                  "b14","b24","b34","s4"),4,4)

##### Model setup for MARSS 

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = QQ, # 1 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # observation error turned on 
  ### initial conditions ###
  #x0 = matrix(x0_fixed),
  V0="zero" ,
  tinitx=0
)

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_092123_1state_nh4_mBFGS.rds")

##### Diagnoses 

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Compare to null model ###
# make sure this matches the fitted model in all ways besides the inclusion of C and c
mod_list_null <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  #C = CC, 
  #c = dat_cov,
  Q = QQ, 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # observation error turned on 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit, output = c("AICc")) # 1187.583 
MARSSaic(null.fit, output = c("AICc")) # 1275.519

### Autoplot diagnoses ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Looks fine

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Eh, HO00
# still doesn't look great but better!

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? Some
# weirdness happening at HO00.

#### 1y legacy, 4 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_092123.rds")

# pivot wider for MARSS format
dat_nh4 <- dat %>%
  select(
    site, index, 
    vwm_nh4, 
    cumulative_precip_mm, 
    fire_perc_ws_1ylegacy, fire_perc_ws_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_nh4, cumulative_precip_mm, 
                    fire_perc_ws_1ylegacy, fire_perc_ws_ppt_1ylegacy))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 218
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,
  0,0,"cumulative_precip_mm_HO00",0,
  0,0,0,"cumulative_precip_mm_RS02",
  # fire_perc: fire effect in 1 y window
  "fire_perc_ws_1ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",
  # fire_perc_ws_ppt: interaction of cum. ppt with fire in 1 y window
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02"), 4, 12)

##### Model setup for MARSS 

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", # 4 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # observation error turned on
  ### initial conditions ###
  #x0 = matrix(x0_fixed),
  V0="zero" ,
  tinitx=0
)

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_092123_4state_nh4_1ylegacy_mBFGS.rds")

##### Diagnoses 

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Compare to null model ###
# make sure this matches the fitted model in all ways besides the inclusion of C and c
mod_list_null <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  #C = CC, 
  #c = dat_cov,
  Q = "diagonal and unequal", # 4 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # observation error turned on
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit, output = c("AICc")) # 1228.187 
MARSSaic(null.fit, output = c("AICc")) # 1307.385

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? HO00 weird again

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Eh.
# HO00 still wonky.

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 1y legacy, 1 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_092123.rds")

# pivot wider for MARSS format
dat_nh4 <- dat %>%
  select(
    site, index, 
    vwm_nh4, 
    cumulative_precip_mm, 
    fire_perc_ws_1ylegacy, fire_perc_ws_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_nh4, cumulative_precip_mm, 
                    fire_perc_ws_1ylegacy, fire_perc_ws_ppt_1ylegacy))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 218
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,
  0,0,"cumulative_precip_mm_HO00",0,
  0,0,0,"cumulative_precip_mm_RS02",
  # fire_perc: fire effect in 1 y window
  "fire_perc_ws_1ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",
  # fire_perc_ws_ppt: interaction of cum. ppt with fire in 1 y window
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02"), 4, 12)

# Make Q matrix
# see page 28 for Q matrix structure:
# https://cran.r-project.org/web/packages/MARSS/vignettes/UserGuide.pdf
QQ <- matrix(list("s1","b12","b13","b14",
                  "b12","s2","b23","b24",
                  "b13","b23","s3","b34",
                  "b14","b24","b34","s4"), 4, 4)

##### Model setup for MARSS 

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = QQ, # 1 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # observation error turned on
  ### initial conditions ###
  #x0 = matrix(x0_fixed),
  V0="zero" ,
  tinitx=0
)

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_092123_1state_nh4_1ylegacy_mBFGS.rds")

##### Diagnoses 

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Compare to null model ###
# make sure this matches the fitted model in all ways besides the inclusion of C and c
mod_list_null <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  #C = CC, 
  #c = dat_cov,
  Q = QQ, # 1 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # observation error turned on
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit, output = c("AICc")) # 1221.692
MARSSaic(null.fit, output = c("AICc")) # 1275.519

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Very few

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Yes.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? 
# For the most part.

#### 2y legacy, 4 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_092123.rds")

# pivot wider for MARSS format
dat_nh4 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_nh4, 
    cumulative_precip_mm, 
    fire_perc_ws_2ylegacy, fire_perc_ws_ppt_2ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_nh4, cumulative_precip_mm, 
                    fire_perc_ws_2ylegacy, fire_perc_ws_ppt_2ylegacy))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 218
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,
  0,0,"cumulative_precip_mm_HO00",0,
  0,0,0,"cumulative_precip_mm_RS02",
  # fire_perc: fire effect in 2 y window
  "fire_perc_ws_2ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_2ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_2ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_2ylegacy_RS02",
  # fire_perc_ws_ppt: interaction of cum. ppt with fire in 2 y window
  "fire_perc_ws_ppt_2ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_ppt_2ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_ppt_2ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_ppt_2ylegacy_RS02"), 4, 12)

##### Model setup for MARSS 

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", # 4 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # observation error turned on 
  ### initial conditions ###
  #x0 = matrix(x0_fixed),
  V0="zero" ,
  tinitx=0
)

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_092123_4state_nh4_2ylegacy_mBFGS.rds")

##### Diagnoses 

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Compare to null model ###
# make sure this matches the fitted model in all ways besides the inclusion of C and c
mod_list_null <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  #C = CC, 
  #c = dat_cov,
  Q = "diagonal and unequal", # 4 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # observation error turned on 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit, output = c("AICc")) # 1265.013 
MARSSaic(null.fit, output = c("AICc")) # 1307.385

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Ho00 looks bad again.

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? HO00 still bad.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 2y legacy, 1 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_092123.rds")

# pivot wider for MARSS format
dat_nh4 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_nh4, 
    cumulative_precip_mm, 
    fire_perc_ws_2ylegacy, fire_perc_ws_ppt_2ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_nh4, cumulative_precip_mm, 
                    fire_perc_ws_2ylegacy, fire_perc_ws_ppt_2ylegacy))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 218
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,
  0,0,"cumulative_precip_mm_HO00",0,
  0,0,0,"cumulative_precip_mm_RS02",
  # fire_perc: fire effect in 2 y window
  "fire_perc_ws_2ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_2ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_2ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_2ylegacy_RS02",
  # fire_perc_ws_ppt: interaction of cum. ppt with fire in 2 y window
  "fire_perc_ws_ppt_2ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_ppt_2ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_ppt_2ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_ppt_2ylegacy_RS02"), 4, 12)

# Make Q matrix
# see page 28 for Q matrix structure:
# https://cran.r-project.org/web/packages/MARSS/vignettes/UserGuide.pdf
QQ <- matrix(list("s1","b12","b13","b14",
                  "b12","s2","b23","b24",
                  "b13","b23","s3","b34",
                  "b14","b24","b34","s4"), 4, 4)

##### Model setup for MARSS 

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = QQ, # 1 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # observation error turned on
  ### initial conditions ###
  #x0 = matrix(x0_fixed),
  V0="zero" ,
  tinitx=0
)

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_092123_1state_nh4_2ylegacy_mBFGS.rds")

##### Diagnoses 

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Compare to null model ###
# make sure this matches the fitted model in all ways besides the inclusion of C and c
mod_list_null <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  #C = CC, 
  #c = dat_cov,
  Q = QQ, # 1 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # observation error turned on
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit, output = c("AICc")) # 1253.181 
MARSSaic(null.fit, output = c("AICc")) # 1275.519

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Looks ok!

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Looks ok!

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 3y legacy, 4 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_092123.rds")

# pivot wider for MARSS format
dat_nh4 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_nh4, 
    cumulative_precip_mm, 
    fire_perc_ws_3ylegacy, fire_perc_ws_ppt_3ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_nh4, cumulative_precip_mm, 
                    fire_perc_ws_3ylegacy, fire_perc_ws_ppt_3ylegacy))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 218
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,
  0,0,"cumulative_precip_mm_HO00",0,
  0,0,0,"cumulative_precip_mm_RS02",
  # fire_perc: fire effect in 3 y window
  "fire_perc_ws_3ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_3ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_3ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_3ylegacy_RS02",
  # fire_perc_ws_ppt: interaction of cum. ppt with fire in 3 y window
  "fire_perc_ws_ppt_3ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_ppt_3ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_ppt_3ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_ppt_3ylegacy_RS02"), 4, 12)

##### Model setup for MARSS 

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", # 4 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # observation error turned on 
  ### initial conditions ###
  #x0 = matrix(x0_fixed),
  V0="zero" ,
  tinitx=0
)

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_092123_4state_nh4_3ylegacy_mBFGS.rds")

##### Diagnoses 

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Compare to null model ###
# make sure this matches the fitted model in all ways besides the inclusion of C and c
mod_list_null <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  #C = CC, 
  #c = dat_cov,
  Q = "diagonal and unequal", # 4 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # observation error turned on 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit, output = c("AICc")) # 1264.527
MARSSaic(null.fit, output = c("AICc")) # 1307.385

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes most fall within CIs.

# Plot 4 (std.model.resids.ytT): DDo resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes.

# Plot 5 (std.state.resids.xtT): Any outliers? Looks alright!

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal? HO00 still worst but fine.

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No.

#### 3y legacy, 1 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_092123.rds")

# pivot wider for MARSS format
dat_nh4 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_nh4, 
    cumulative_precip_mm, 
    fire_perc_ws_3ylegacy, fire_perc_ws_ppt_3ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_nh4, cumulative_precip_mm, 
                    fire_perc_ws_3ylegacy, fire_perc_ws_ppt_3ylegacy))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 218
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,
  0,0,"cumulative_precip_mm_HO00",0,
  0,0,0,"cumulative_precip_mm_RS02",
  # fire_perc: fire effect in 3 y window
  "fire_perc_ws_3ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_3ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_3ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_3ylegacy_RS02",
  # fire_perc_ws_ppt: interaction of cum. ppt with fire in 3 y window
  "fire_perc_ws_ppt_3ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_ppt_3ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_ppt_3ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_ppt_3ylegacy_RS02"), 4, 12)

# Make Q matrix
# see page 28 for Q matrix structure:
# https://cran.r-project.org/web/packages/MARSS/vignettes/UserGuide.pdf
QQ <- matrix(list("s1","b12","b13","b14",
                  "b12","s2","b23","b24",
                  "b13","b23","s3","b34",
                  "b14","b24","b34","s4"), 4, 4)

##### Model setup for MARSS 

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = QQ, # 1 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # observation error turned on 
  ### initial conditions ###
  #x0 = matrix(x0_fixed),
  V0="zero" ,
  tinitx=0
)

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_092123_1state_nh4_3ylegacy_mBFGS.rds")

##### Diagnoses 

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Compare to null model ###
# make sure this matches the fitted model in all ways besides the inclusion of C and c
mod_list_null <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  #C = CC, 
  #c = dat_cov,
  Q = QQ, # 1 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit, output = c("AICc")) # 1252.593 
MARSSaic(null.fit, output = c("AICc")) # 1275.519

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes most fall within CIs.

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes most fall within CIs.

# Plot 5 (std.state.resids.xtT): Any outliers? Looks alright!

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal? Ok.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No.

#### 4y legacy, 4 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_092123.rds")

# pivot wider for MARSS format
dat_nh4 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_nh4, 
    cumulative_precip_mm, 
    fire_perc_ws_4ylegacy, fire_perc_ws_ppt_4ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_nh4, cumulative_precip_mm, 
                    fire_perc_ws_4ylegacy, fire_perc_ws_ppt_4ylegacy))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 218
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,
  0,0,"cumulative_precip_mm_HO00",0,
  0,0,0,"cumulative_precip_mm_RS02",
  # fire_perc: fire effect in 4 y window
  "fire_perc_ws_4ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_4ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_4ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_4ylegacy_RS02",
  # fire_perc_ws_ppt: interaction of cum. ppt with fire in 4 y window
  "fire_perc_ws_ppt_4ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_ppt_4ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_ppt_4ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_ppt_4ylegacy_RS02"), 4, 12)

##### Model setup for MARSS 

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", # 4 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # observation error turned on
  ### initial conditions ###
  #x0 = matrix(x0_fixed),
  V0="zero" ,
  tinitx=0
)

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_092123_4state_nh4_4ylegacy_mBFGS.rds")

##### Diagnoses 

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Compare to null model ###
# make sure this matches the fitted model in all ways besides the inclusion of C and c
mod_list_null <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  #C = CC, 
  #c = dat_cov,
  Q = "diagonal and unequal", # 4 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # obs. error turned on 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit, output = c("AICc")) # 1263.847 
MARSSaic(null.fit, output = c("AICc")) # 1307.385

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes 

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes 

# Plot 5 (std.state.resids.xtT): Any outliers? Looks ok!

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? HO00
# still meh, all others fine.

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 4y legacy, 1 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_092123.rds")

# pivot wider for MARSS format
dat_nh4 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_nh4, 
    cumulative_precip_mm, 
    fire_perc_ws_4ylegacy, fire_perc_ws_ppt_4ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_nh4, cumulative_precip_mm, 
                    fire_perc_ws_4ylegacy, fire_perc_ws_ppt_4ylegacy))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 218
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,
  0,0,"cumulative_precip_mm_HO00",0,
  0,0,0,"cumulative_precip_mm_RS02",
  # fire_perc: fire effect in 4 y window
  "fire_perc_ws_4ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_4ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_4ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_4ylegacy_RS02",
  # fire_perc_ws_ppt: interaction of cum. ppt with fire in 4 y window
  "fire_perc_ws_ppt_4ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_ppt_4ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_ppt_4ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_ppt_4ylegacy_RS02"), 4, 12)

# Make Q matrix
# see page 28 for Q matrix structure:
# https://cran.r-project.org/web/packages/MARSS/vignettes/UserGuide.pdf
QQ <- matrix(list("s1","b12","b13","b14",
                  "b12","s2","b23","b24",
                  "b13","b23","s3","b34",
                  "b14","b24","b34","s4"), 4, 4)

##### Model setup for MARSS 

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = QQ, # 1 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # obs. error
  ### initial conditions ###
  #x0 = matrix(x0_fixed),
  V0="zero" ,
  tinitx=0
)

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_092123_1state_nh4_4ylegacy_mBFGS.rds")

##### Diagnoses 

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Compare to null model ###
# make sure this matches the fitted model in all ways besides the inclusion of C and c
mod_list_null <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  #C = CC, 
  #c = dat_cov,
  Q = QQ, # 1 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # obs. error 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit, output = c("AICc")) # 1255.154 
MARSSaic(null.fit, output = c("AICc")) # 1275.519

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes 

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Looks ok!

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? HO00
# still meh, all others fine.

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 5y legacy, 4 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_092123.rds")

# pivot wider for MARSS format
dat_nh4 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_nh4, 
    cumulative_precip_mm, 
    fire_perc_ws_5ylegacy, fire_perc_ws_ppt_5ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_nh4, cumulative_precip_mm, 
                    fire_perc_ws_5ylegacy, fire_perc_ws_ppt_5ylegacy))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 218
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,
  0,0,"cumulative_precip_mm_HO00",0,
  0,0,0,"cumulative_precip_mm_RS02",
  # fire_perc: fire effect in 5 y window
  "fire_perc_ws_5ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_5ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_5ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_5ylegacy_RS02",
  # fire_perc_ws_ppt: interaction of cum. ppt with fire in 5 y window
  "fire_perc_ws_ppt_5ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_ppt_5ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_ppt_5ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_ppt_5ylegacy_RS02"), 4, 12)

##### Model setup for MARSS 

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", # 4 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # obs. error 
  ### initial conditions ###
  #x0 = matrix(x0_fixed),
  V0="zero" ,
  tinitx=0
)

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_092123_4state_nh4_5ylegacy_mBFGS.rds")

##### Diagnoses 

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Compare to null model ###
# make sure this matches the fitted model in all ways besides the inclusion of C and c
mod_list_null <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  #C = CC, 
  #c = dat_cov,
  Q = "diagonal and unequal", # 4 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # obs. error 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit, output = c("AICc")) # 1265.634
MARSSaic(null.fit, output = c("AICc")) # 1307.385

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes most fall within CIs.

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes most fall within CIs.

# Plot 5 (std.state.resids.xtT): Any outliers? Looks ok!

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Fine.

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 5y legacy, 1 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_092123.rds")

# pivot wider for MARSS format
dat_nh4 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_nh4, 
    cumulative_precip_mm, 
    fire_perc_ws_5ylegacy, fire_perc_ws_ppt_5ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_nh4, cumulative_precip_mm, 
                    fire_perc_ws_5ylegacy, fire_perc_ws_ppt_5ylegacy))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 218
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,
  0,0,"cumulative_precip_mm_HO00",0,
  0,0,0,"cumulative_precip_mm_RS02",
  # fire_perc: fire effect in 5 y window
  "fire_perc_ws_5ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_5ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_5ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_5ylegacy_RS02",
  # fire_perc_ws_ppt: interaction of cum. ppt with fire in 5 y window
  "fire_perc_ws_ppt_5ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_ppt_5ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_ppt_5ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_ppt_5ylegacy_RS02"), 4, 12)

# Make Q matrix
# see page 28 for Q matrix structure:
# https://cran.r-project.org/web/packages/MARSS/vignettes/UserGuide.pdf
QQ <- matrix(list("s1","b12","b13","b14",
                  "b12","s2","b23","b24",
                  "b13","b23","s3","b34",
                  "b14","b24","b34","s4"), 4, 4)

##### Model setup for MARSS 

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = QQ, # 1 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # obs. error 
  ### initial conditions ###
  #x0 = matrix(x0_fixed),
  V0="zero" ,
  tinitx=0
)

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_092123_1state_nh4_5ylegacy_mBFGS.rds")

##### Diagnoses 

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Compare to null model ###
# make sure this matches the fitted model in all ways besides the inclusion of C and c
mod_list_null <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  #C = CC, 
  #c = dat_cov,
  Q = QQ, # 1 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "diagonal and equal", # obs. error
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit, output = c("AICc")) # 1262.697 
MARSSaic(null.fit, output = c("AICc")) # 1275.519

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes most fall within CIs.

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes most fall within CIs.

# Plot 5 (std.state.resids.xtT): Any outliers? Looks ok!

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Ok.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No.

#### AICc Comparisons ####

# Compare all model fits for each legacy window to see which state 
# configuration was best.

# Presented here are the information criterion:
# AICc - Akaike Information Criterion adjusted for small sample sizes

# no legacy, 4 state
noleg_4state <- readRDS(file = "data_working/marss_fits/fit_092123_4state_nh4_mBFGS.rds")

# no legacy, 1 state
noleg_1state <- readRDS(file = "data_working/marss_fits/fit_092123_1state_nh4_mBFGS.rds")

MARSSaic(noleg_4state, output = c("AICc")) # 1200.804 
MARSSaic(noleg_1state, output = c("AICc")) # 1187.583

# So, based on lowest AICc values, 1-state wins.

###

# 1y legacy, 4 state
leg1_4state <- readRDS(file = "data_working/marss_fits/fit_092123_4state_nh4_1ylegacy_mBFGS.rds")

# 1y legacy, 1 state
leg1_1state <- readRDS(file = "data_working/marss_fits/fit_092123_1state_nh4_1ylegacy_mBFGS.rds")

MARSSaic(leg1_4state, output = c("AICc")) # 1228.187
MARSSaic(leg1_1state, output = c("AICc")) # 1221.692

# So, based on lowest AICc values, 1-state wins.

###

# 2y legacy, 4 state
leg2_4state <- readRDS(file = "data_working/marss_fits/fit_092123_4state_nh4_2ylegacy_mBFGS.rds")

# 2y legacy, 1 state
leg2_1state <- readRDS(file = "data_working/marss_fits/fit_092123_1state_nh4_2ylegacy_mBFGS.rds")

MARSSaic(leg2_4state, output = c("AICc")) # 1265.013
MARSSaic(leg2_1state, output = c("AICc")) # 1253.181

# So, based on lowest AICc values, 1-state wins.

###

# 3y legacy, 4 state
leg3_4state <- readRDS(file = "data_working/marss_fits/fit_092123_4state_nh4_3ylegacy_mBFGS.rds")

# 3y legacy, 1 state
leg3_1state <- readRDS(file = "data_working/marss_fits/fit_092123_1state_nh4_3ylegacy_mBFGS.rds")

MARSSaic(leg3_4state, output = c("AICc")) # 1264.527
MARSSaic(leg3_1state, output = c("AICc")) # 1252.593

# So, based on lowest AICc values, 1-state wins.

###

# 4y legacy, 4 state
leg4_4state <- readRDS(file = "data_working/marss_fits/fit_092123_4state_nh4_4ylegacy_mBFGS.rds")

# 4y legacy, 1 state
leg4_1state <- readRDS(file = "data_working/marss_fits/fit_092123_1state_nh4_4ylegacy_mBFGS.rds")

MARSSaic(leg4_4state, output = c("AICc")) # 1263.847
MARSSaic(leg4_1state, output = c("AICc")) # 1255.154

# So, based on lowest AICc values, 1-state wins.

###

# 5y legacy, 4 state
leg5_4state <- readRDS(file = "data_working/marss_fits/fit_092123_4state_nh4_5ylegacy_mBFGS.rds")

# 5y legacy, 1 state
leg5_1state <- readRDS(file = "data_working/marss_fits/fit_092123_1state_nh4_5ylegacy_mBFGS.rds")

MARSSaic(leg5_4state, output = c("AICc")) # 1265.634
MARSSaic(leg5_1state, output = c("AICc")) # 1262.697

# So, based on lowest AICc values, 1-state wins.

###

# So, based on AICc values, it would seem the 1 "state" model structure wins out.

#### Results Figure ####

# For presentation consistency, I will only be creating figures with a
# single state configuration, whichever yielded the most parsimonious
# models. So, in this case, all NH4 figures will represent the "1 state"
# scenario.

# Extract necessary confidence interval info
noleg_est <- MARSSparamCIs(noleg_1state, alpha = 0.01)
leg1y_est <- MARSSparamCIs(leg1_1state, alpha = 0.01)
leg2y_est <- MARSSparamCIs(leg2_1state, alpha = 0.01)
leg3y_est <- MARSSparamCIs(leg3_1state, alpha = 0.01)
leg4y_est <- MARSSparamCIs(leg4_1state, alpha = 0.01)
leg5y_est <- MARSSparamCIs(leg5_1state, alpha = 0.01)

# Format confidence intervals into dataframes
noleg_CI = data.frame(
  "Est." = noleg_est$par$U,
  "Lower" = noleg_est$par.lowCI$U,
  "Upper" = noleg_est$par.upCI$U)
noleg_CI$Parameter = rownames(noleg_CI)
noleg_CI[,1:3] = round(noleg_CI[,1:3], 3)
noleg_CI$Model = "immediate window"

leg1y_CI = data.frame(
  "Est." = leg1y_est$par$U,
  "Lower" = leg1y_est$par.lowCI$U,
  "Upper" = leg1y_est$par.upCI$U)
leg1y_CI$Parameter = rownames(leg1y_CI)
leg1y_CI[,1:3] = round(leg1y_CI[,1:3], 3)
leg1y_CI$Model = "1 year window"

leg2y_CI = data.frame(
  "Est." = leg2y_est$par$U,
  "Lower" = leg2y_est$par.lowCI$U,
  "Upper" = leg2y_est$par.upCI$U)
leg2y_CI$Parameter = rownames(leg2y_CI)
leg2y_CI[,1:3] = round(leg2y_CI[,1:3], 3)
leg2y_CI$Model = "2 year window"

leg3y_CI = data.frame(
  "Est." = leg3y_est$par$U,
  "Lower" = leg3y_est$par.lowCI$U,
  "Upper" = leg3y_est$par.upCI$U)
leg3y_CI$Parameter = rownames(leg3y_CI)
leg3y_CI[,1:3] = round(leg3y_CI[,1:3], 3)
leg3y_CI$Model = "3 year window"

leg4y_CI = data.frame(
  "Est." = leg4y_est$par$U,
  "Lower" = leg4y_est$par.lowCI$U,
  "Upper" = leg4y_est$par.upCI$U)
leg4y_CI$Parameter = rownames(leg4y_CI)
leg4y_CI[,1:3] = round(leg4y_CI[,1:3], 3)
leg4y_CI$Model = "4 year window"

leg5y_CI = data.frame(
  "Est." = leg5y_est$par$U,
  "Lower" = leg5y_est$par.lowCI$U,
  "Upper" = leg5y_est$par.upCI$U)
leg5y_CI$Parameter = rownames(leg5y_CI)
leg5y_CI[,1:3] = round(leg5y_CI[,1:3], 3)
leg5y_CI$Model = "5 year window"

# Bind all together
CIs = rbind(noleg_CI, leg1y_CI, leg2y_CI, leg3y_CI, leg4y_CI,leg5y_CI)

# Add column for site names
CIs$Stream = gsub("_","",str_sub(CIs$Parameter, start= -4))

# Simplify parameter names
CIs$Parm_simple = c(rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4))

# Add column to designate those sites at which effects are significant.
CIs <- CIs %>%
  mutate(sig = factor(case_when(`Est.` > 0 & 
                                  Lower > 0 & 
                                  Upper > 0 ~ "sig_pos",
                                `Est.` < 0 & 
                                  Lower < 0 & 
                                  Upper < 0 ~ "sig_neg",
                                TRUE ~ "not_sig"), 
                      levels = c("sig_pos", "not_sig", "sig_neg"))) %>%
  mutate(model = factor(Model, levels = c("immediate window",
                                          "1 year window",
                                          "2 year window",
                                          "3 year window",
                                          "4 year window",
                                          "5 year window")))

my_palette <- c("black", "white", "black")

# Plot results
(NH4_fig <- ggplot(CIs, aes(x = factor(Parm_simple, 
                                       levels = c("Ppt x Perc. burn",
                                                  "Perc. burn",
                                                  "Ppt")),
                            y = Est., fill = sig, shape = Stream)) + 
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
                  position=position_dodge(width = 0.5), width = 0) +
    geom_point(position=position_dodge(width = 0.5), 
               alpha = 0.8, size = 8) + 
    scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
    scale_fill_manual(values = my_palette) +
    theme_bw()+
    theme(axis.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          strip.text.x = element_text(size = 24),
          strip.text.y = element_text(size = 24),
          legend.title=element_text(size = 20), 
          legend.text=element_text(size = 20)) +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    coord_flip(ylim = c(-1, 1)) + 
    labs(y = "Effect Size", 
         x = "Covariates",
         fill = "Significance") +
    theme(plot.margin = unit(c(.2,.2,.05,.05),"cm")) + 
    guides(shape = guide_legend("Stream"), fill = "none") +
    facet_grid(.~model))

saveRDS(NH4_fig, "figures/MARSS_NH4_1state_101023.rds")

# Read in all three plots.
nh4 <- readRDS("figures/MARSS_NH4_1state_101023.rds")
no3 <- readRDS("figures/MARSS_NO3_1state_101023.rds")
po4 <- readRDS("figures/MARSS_PO4_1state_101023.rds")

# Combine into a single figure.
(nuts <- nh4 + no3 + po4 +
    plot_layout(nrow = 3) +
    plot_annotation(tag_levels = 'A'))

# Export plot.
# ggsave(("MARSS_nutrients_101023.png"),
#        path = "figures",
#        width = 65,
#        height = 50,
#        units = "cm"
# )

# End of script.
