# Nutrient MARSS models
# with fire x ppt interactions and legacy effects
# as well as 4 state structure for CA sites
# Script started August 9, 2023 by Heili Lowman

# This script will run 12 NO3 MARSS models.
# Note, each model fit will remove all stored data, until you reach the AIC
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

#### 0y legacy, 4 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_no3 <- dat %>%
  select(
    site, index, 
    vwm_no3, 
    cumulative_precip_mm, 
    fire_perc_ws, fire_perc_ws_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_no3, cumulative_precip_mm, 
                    fire_perc_ws, fire_perc_ws_ppt))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 213
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
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
  # fire_perc_ws_ppt: interaction of cum. ppt with fire in 6 m window
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
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
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
        file = "data_working/marss_fits/fit_080923_4state_no3_mBFGS.rds")

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
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#           dAIC df
# fit        0.0 24
# null.fit 103.4 12
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do all equal zero because we have nothing in the observation model? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Very few

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (a.k.a. straight line)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# No 

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 0y legacy, 1 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_no3 <- dat %>%
  select(
    site, index, 
    vwm_no3, 
    cumulative_precip_mm, 
    fire_perc_ws, fire_perc_ws_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_no3, cumulative_precip_mm, 
                    fire_perc_ws, fire_perc_ws_ppt))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 213
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
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
  # fire_perc_ws_ppt: interaction of cum. ppt with fire in 6 m window
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
  Q = QQ, 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
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
        file = "data_working/marss_fits/fit_080923_1state_no3_mBFGS.rds")

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
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#           dAIC df
# fit        0.0 30
# null.fit  46.7 18
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do all equal zero because we have nothing in the observation model? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Very few

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (a.k.a. straight line)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# No 

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 1y legacy, 4 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_no3 <- dat %>%
  select(
    site, index, 
    vwm_no3, 
    cumulative_precip_mm, 
    fire_perc_ws_1ylegacy, fire_perc_ws_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_no3, cumulative_precip_mm, 
                    fire_perc_ws_1ylegacy, fire_perc_ws_ppt_1ylegacy))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 213
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
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
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
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
        file = "data_working/marss_fits/fit_080923_4state_no3_1ylegacy_mBFGS.rds")

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
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#           dAIC df
# fit        0.0 24
# null.fit 124.6 12
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do all equal zero because we have nothing in the observation model? Yep!

# Plot 5 (std.state.resids.xtT): Any detectable outliers? Very few

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight line)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 1y legacy, 1 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_no3 <- dat %>%
  select(
    site, index, 
    vwm_no3, 
    cumulative_precip_mm, 
    fire_perc_ws_1ylegacy, fire_perc_ws_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_no3, cumulative_precip_mm, 
                    fire_perc_ws_1ylegacy, fire_perc_ws_ppt_1ylegacy))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 213
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
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
  Q = QQ, 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
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
        file = "data_working/marss_fits/fit_080923_1state_no3_1ylegacy_mBFGS.rds")

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
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#           dAIC df
# fit        0.0 30
# null.fit  66.7 18
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do all equal zero because we have nothing in the observation model? Yes

# Plot 5 (std.state.resids.xtT): Any detectable outliers? Very few

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 2y legacy, 4 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_no3 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_no3, 
    cumulative_precip_mm, 
    fire_perc_ws_2ylegacy, fire_perc_ws_ppt_2ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_no3, cumulative_precip_mm, 
                    fire_perc_ws_2ylegacy, fire_perc_ws_ppt_2ylegacy))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 213
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
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
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
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
        file = "data_working/marss_fits/fit_080923_4state_no3_2ylegacy_mBFGS.rds")

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
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#           dAIC df
# fit        0.0 24
# null.fit 132.6 12
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No 
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do all equal zero because we have nothing in the observation model? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Very few

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight line)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 2y legacy, 1 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_no3 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_no3, 
    cumulative_precip_mm, 
    fire_perc_ws_2ylegacy, fire_perc_ws_ppt_2ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_no3, cumulative_precip_mm, 
                    fire_perc_ws_2ylegacy, fire_perc_ws_ppt_2ylegacy))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 213
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
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
  Q = QQ, 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
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
        file = "data_working/marss_fits/fit_080923_1state_no3_2ylegacy_mBFGS.rds")

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
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#           dAIC df
# fit        0.0 30
# null.fit  75.7 18
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No 
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do all equal zero because we have nothing in the observation model? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Very few

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight line)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 3y legacy, 4 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_no3 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_no3, 
    cumulative_precip_mm, 
    fire_perc_ws_3ylegacy, fire_perc_ws_ppt_3ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_no3, cumulative_precip_mm, 
                    fire_perc_ws_3ylegacy, fire_perc_ws_ppt_3ylegacy))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 213
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
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
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
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
        file = "data_working/marss_fits/fit_080923_4state_no3_3ylegacy_mBFGS.rds")

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
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#           dAIC df
# fit        0.0 24
# null.fit 134.4 12
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero because we have nothing in the observation model? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Very few

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight line)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 3y legacy, 1 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_no3 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_no3, 
    cumulative_precip_mm, 
    fire_perc_ws_3ylegacy, fire_perc_ws_ppt_3ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_no3, cumulative_precip_mm, 
                    fire_perc_ws_3ylegacy, fire_perc_ws_ppt_3ylegacy))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 213
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
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
  Q = QQ, 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
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
        file = "data_working/marss_fits/fit_080923_1state_no3_3ylegacy_mBFGS.rds")

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
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#           dAIC df
# fit        0.0 30
# null.fit  74.5 18
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero because we have nothing in the observation model? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Very few

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight line)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 4y legacy, 4 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_no3 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_no3, 
    cumulative_precip_mm, 
    fire_perc_ws_4ylegacy, fire_perc_ws_ppt_4ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_no3, cumulative_precip_mm, 
                    fire_perc_ws_4ylegacy, fire_perc_ws_ppt_4ylegacy))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 213
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
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
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
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
        file = "data_working/marss_fits/fit_080923_4state_no3_4ylegacy_mBFGS.rds")

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
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#           dAIC df
# fit        0.0 24
# null.fit 131.5 12
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero because we have nothing in the observation model? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Very few

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight line)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 4y legacy, 1 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_no3 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_no3, 
    cumulative_precip_mm, 
    fire_perc_ws_4ylegacy, fire_perc_ws_ppt_4ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_no3, cumulative_precip_mm, 
                    fire_perc_ws_4ylegacy, fire_perc_ws_ppt_4ylegacy))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 213
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
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
  Q = QQ, 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
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
        file = "data_working/marss_fits/fit_080923_1state_no3_4ylegacy_mBFGS.rds")

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
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#           dAIC df
# fit        0.0 30
# null.fit  72.1 18
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero because we have nothing in the observation model? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Very few

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight line)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 5y legacy, 4 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_no3 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_no3, 
    cumulative_precip_mm, 
    fire_perc_ws_5ylegacy, fire_perc_ws_ppt_5ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_no3, cumulative_precip_mm, 
                    fire_perc_ws_5ylegacy, fire_perc_ws_ppt_5ylegacy))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 213
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
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
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
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
        file = "data_working/marss_fits/fit_080923_4state_no3_5ylegacy_mBFGS.rds")

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
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#           dAIC df
# fit        0.0 24
# null.fit 127.9 12
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero because we have nothing in the observation model? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Very few

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straigh lines)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 5y legacy, 1 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_no3 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_no3, 
    cumulative_precip_mm, 
    fire_perc_ws_5ylegacy, fire_perc_ws_ppt_5ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_no3, cumulative_precip_mm, 
                    fire_perc_ws_5ylegacy, fire_perc_ws_ppt_5ylegacy))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 213
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
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
  Q = QQ, 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
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
        file = "data_working/marss_fits/fit_080923_1state_no3_5ylegacy_mBFGS.rds")

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
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#           dAIC df
# fit        0.0 30
# null.fit  68.6 18
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero because we have nothing in the observation model? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Very few

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

### Overall ###
# None of these diagnoses look prohibitively bad.

#### IC Comparisons ####

# Compare all model fits for each legacy window to see which state 
# configuration was best.

# Presented here are three information criterion:
# AIC - Akaike Information Criterion
# AICc - Akaike Information Criterion adjusted for small sample sizes
# BIC - Bayesian Information Criterion

# Per the recommendations found in Brewer et al. 2016 (doi: 10.1111/2041-210X.12541) 
# we will be using BIC for most parsimonious model selection, although
# all three are displayed here for transparency/comparison.

# no legacy, 4 state
noleg_4state <- readRDS(file = "data_working/marss_fits/fit_080923_4state_no3_mBFGS.rds")

# no legacy, 1 state
noleg_1state <- readRDS(file = "data_working/marss_fits/fit_080923_1state_no3_mBFGS.rds")

bbmle::AICtab(noleg_4state, noleg_1state)

#              dAIC df
# noleg_1state  0.0 30
# noleg_4state  1.3 24

broom::glance(noleg_4state) # AICc 1251.948
stats4::BIC(noleg_4state) # BIC 1351.152

broom::glance(noleg_1state) # AICc 1252.052  
stats4::BIC(noleg_1state) # BIC 1375.268

# So, based on lowest AIC value, 1-state wins, but based on lowest AICc and BIC
# values, 4-state wins.

###

# 1y legacy, 4 state
leg1_4state <- readRDS(file = "data_working/marss_fits/fit_080923_4state_no3_1ylegacy_mBFGS.rds")

# 1y legacy, 1 state
leg1_1state <- readRDS(file = "data_working/marss_fits/fit_080923_1state_no3_1ylegacy_mBFGS.rds")

bbmle::AICtab(leg1_4state, leg1_1state)

#             dAIC df
# leg1_4state  0   24
# leg1_1state  0   30

#              dAIC df
# noleg_1state  0.0 30
# noleg_4state  1.3 24

broom::glance(leg1_4state) # AICc 1230.721
stats4::BIC(leg1_4state) # BIC 1329.924

broom::glance(leg1_1state) # AICc 1232.141 
stats4::BIC(leg1_1state) # BIC 1355.357

# So, based on lowest AIC value, 1-state wins, but based on lowest AICc and BIC
# values, 4-state wins.

###

# 2y legacy, 4 state
leg2_4state <- readRDS(file = "data_working/marss_fits/fit_080923_4state_no3_2ylegacy_mBFGS.rds")

# 2y legacy, 1 state
leg2_1state <- readRDS(file = "data_working/marss_fits/fit_080923_1state_no3_2ylegacy_mBFGS.rds")

bbmle::AICtab(leg2_4state, leg2_1state)

#             dAIC df
# leg2_1state  0   30
# leg2_4state  1   24

broom::glance(leg2_4state) # AICc 1222.743
stats4::BIC(leg2_4state) # BIC 1321.947

broom::glance(leg2_1state) # AICc 1223.121
stats4::BIC(leg2_1state) # BIC 1346.337

# So, based on lowest AIC value, 1-state wins, but based on lowest AICc and BIC
# values, 4-state wins.

###

# 3y legacy, 4 state
leg3_4state <- readRDS(file = "data_working/marss_fits/fit_080923_4state_no3_3ylegacy_mBFGS.rds")

# 3y legacy, 1 state
leg3_1state <- readRDS(file = "data_working/marss_fits/fit_080923_1state_no3_3ylegacy_mBFGS.rds")

bbmle::AICtab(leg3_4state, leg3_1state)

#             dAIC df
# leg3_4state  0.0 24
# leg3_1state  1.9 30

broom::glance(leg3_4state) # AICc 1220.994 
stats4::BIC(leg3_4state) # BIC 1320.197

broom::glance(leg3_1state) # AICc 1224.283
stats4::BIC(leg3_1state) # BIC 1347.499

# So, based on lowest AIC, AICc, and BIC values, 4-state wins.

###

# 4y legacy, 4 state
leg4_4state <- readRDS(file = "data_working/marss_fits/fit_080923_4state_no3_4ylegacy_mBFGS.rds")

# 4y legacy, 1 state
leg4_1state <- readRDS(file = "data_working/marss_fits/fit_080923_1state_no3_4ylegacy_mBFGS.rds")

bbmle::AICtab(leg4_4state, leg4_1state)

#             dAIC df
# leg4_4state  0.0 24
# leg4_1state  1.4 30

broom::glance(leg4_4state) # AICc 1223.863 
stats4::BIC(leg4_4state) # BIC 1323.067

broom::glance(leg4_1state) # AICc 1226.668
stats4::BIC(leg4_1state) # BIC 1349.884

# So, based on lowest AIC, AICc, and BIC values, 4-state wins.

###

# 5y legacy, 4 state
leg5_4state <- readRDS(file = "data_working/marss_fits/fit_080923_4state_no3_5ylegacy_mBFGS.rds")

# 5y legacy, 1 state
leg5_1state <- readRDS(file = "data_working/marss_fits/fit_080923_1state_no3_5ylegacy_mBFGS.rds")

bbmle::AICtab(leg5_4state, leg5_1state)

#             dAIC df
# leg5_4state  0.0 24
# leg5_1state  1.4 30

broom::glance(leg5_4state) # AICc 1227.41
stats4::BIC(leg5_4state) # BIC 1326.613

broom::glance(leg5_1state) # AICc 1230.192
stats4::BIC(leg5_1state) # BIC 1353.408

# So, based on lowest AIC, AICc, and BIC values, 4-state wins.

###

# So, it would seem the 4 "state" model structure wins out.

stats4::BIC(noleg_4state) # BIC 1351.152
stats4::BIC(leg1_4state) # BIC 1329.924
stats4::BIC(leg2_4state) # BIC 1321.947
stats4::BIC(leg3_4state) # BIC 1320.197
stats4::BIC(leg4_4state) # BIC 1323.067
stats4::BIC(leg5_4state) # BIC 1326.613

# And when comparing all models, the 3 year window/lag is most parsimonious.

#### Results Figure ####

# For presentation consistency, I will only be creating figures with a
# single state configuration, whichever yielded the most parsimonious
# models. So, in this case, all NO3 figures will represent the "4 state"
# scenario.

# Extract necessary confidence interval info
noleg_est <- MARSSparamCIs(noleg_4state)
leg1y_est <- MARSSparamCIs(leg1_4state)
leg2y_est <- MARSSparamCIs(leg2_4state)
leg3y_est <- MARSSparamCIs(leg3_4state)
leg4y_est <- MARSSparamCIs(leg4_4state)
leg5y_est <- MARSSparamCIs(leg5_4state)

# Format confidence intervals into dataframes
noleg_CI = data.frame(
  "Est." = noleg_est$par$U,
  "Lower" = noleg_est$par.lowCI$U,
  "Upper" = noleg_est$par.upCI$U)
noleg_CI$Parameter = rownames(noleg_CI)
noleg_CI[,1:3] = round(noleg_CI[,1:3], 3)
noleg_CI$Model = "0 year window"

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
                      levels = c("sig_pos", "not_sig", "sig_neg")))

my_palette <- c("black", "white", "black")

# Plot results
(NO3_fig <- ggplot(CIs, aes(x = factor(Parm_simple, 
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
    facet_grid(.~Model))

# Export plot.
# ggsave(("MARSS_NO3_080923.png"),
#        path = "figures",
#        width = 65,
#        height = 12,
#        units = "cm"
# )

# End of script.