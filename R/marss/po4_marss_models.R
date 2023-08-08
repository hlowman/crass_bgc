# Nutrient MARSS models
# with fire x ppt interactions and legacy effects
# as well as 5 state structure for CA sites
# Script started August 8, 2023 by Heili Lowman

# This script will run 12 PO4 MARSS models.
# Note, each model fit will remove all stored data, until you reach the AIC
# portion of this script.

#### Setup ####

# Load packages.
library(tidyverse)
library(lubridate)
library(MARSS)
library(naniar) 
library(here)

#### 0y legacy, 5 state ####

# Set up data for MARSS
# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_po4 <- dat %>%
  select(
    site, index, 
    vwm_po4, 
    cumulative_precip_mm, 
    fire_perc_ws, fire_perc_ws_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_po4, cumulative_precip_mm, 
                    fire_perc_ws, fire_perc_ws_ppt))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 213
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
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

# "XXXX_AB00",0,0,0,
# 0,"XXXX_GV01",0,0,
# 0,0,"XXXX_HO00",0,
# 0,0,0,"XXXX_RS02",

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
        file = "data_working/marss_fits/fit_080823_5state_po4_mBFGS.rds")

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
# null.fit  99.8 12
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yes!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# These look AHmazing!

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 0y legacy, 1 state ####

# Set up data for MARSS
# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_po4 <- dat %>%
  select(
    site, index, 
    vwm_po4, 
    cumulative_precip_mm, 
    fire_perc_ws, fire_perc_ws_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_po4, cumulative_precip_mm, 
                    fire_perc_ws, fire_perc_ws_ppt))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 213
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
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

# "XXXX_AB00",0,0,0,
# 0,"XXXX_GV01",0,0,
# 0,0,"XXXX_HO00",0,
# 0,0,0,"XXXX_RS02",

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
        file = "data_working/marss_fits/fit_080823_1state_po4_mBFGS.rds")

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
# null.fit    57 18
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes fall within CIs

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yes all zero

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks fine

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# These look great still!

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 1y legacy, 5 state ####

# Set up data for MARSS
# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_po4 <- dat %>%
  select(
    site, index, 
    vwm_po4, 
    cumulative_precip_mm, 
    fire_perc_ws_1ylegacy, fire_perc_ws_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_po4, cumulative_precip_mm, 
                    fire_perc_ws_1ylegacy, fire_perc_ws_ppt_1ylegacy))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 213
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
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
        file = "data_working/marss_fits/fit_080823_5state_po4_1ylegacy_mBFGS.rds")

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
# null.fit 108.6 12
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). All do.

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Look just fine.

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look fine.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 1y legacy, 1 state ####

# Set up data for MARSS
# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat <- readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_po4 <- dat %>%
  select(
    site, index, 
    vwm_po4, 
    cumulative_precip_mm, 
    fire_perc_ws_1ylegacy, fire_perc_ws_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_po4, cumulative_precip_mm, 
                    fire_perc_ws_1ylegacy, fire_perc_ws_ppt_1ylegacy))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 213
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
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
        file = "data_working/marss_fits/fit_080823_1state_po4_1ylegacy_mBFGS.rds")

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
# null.fit  70.4 18
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No; Yes

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). All zero.

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looking fine.

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Also look fine.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 2y legacy, 5 state ####

# Set up data for MARSS
# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_po4 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_po4, 
    cumulative_precip_mm, 
    fire_perc_ws_2ylegacy, fire_perc_ws_ppt_2ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_po4, cumulative_precip_mm, 
                    fire_perc_ws_2ylegacy, fire_perc_ws_ppt_2ylegacy))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 213
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
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
  # fire_perc_ws_ppt: interaction of cum. ppt with fire in 2y window
  "fire_perc_ws_ppt_2ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_ppt_2ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_ppt_2ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_ppt_2ylegacy_RS02"),4,12)

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
        file = "data_working/marss_fits/fit_080823_5state_po4_2ylegacy_mBFGS.rds")

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
# null.fit  98.9 12
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No and Yes most fall within CIs.

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yes.

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look good.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 2y legacy, 1 state ####

# Set up data for MARSS
# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_po4 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_po4, 
    cumulative_precip_mm, 
    fire_perc_ws_2ylegacy, fire_perc_ws_ppt_2ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_po4, cumulative_precip_mm, 
                    fire_perc_ws_2ylegacy, fire_perc_ws_ppt_2ylegacy))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 213
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
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
  # fire_perc_ws_ppt: interaction of cum. ppt with fire in 2y window
  "fire_perc_ws_ppt_2ylegacy_AB00",0,0,0,
  0,"fire_perc_ws_ppt_2ylegacy_GV01",0,0,
  0,0,"fire_perc_ws_ppt_2ylegacy_HO00",0,
  0,0,0,"fire_perc_ws_ppt_2ylegacy_RS02"),4,12)

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
        file = "data_working/marss_fits/fit_080823_1state_po4_2ylegacy_mBFGS.rds")

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
# null.fit  66.9 18
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No and Yes most fall within CIs.

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yes.

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Look ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look fine.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 3y legacy, 5 state ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_po4 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_po4, 
    cumulative_precip_mm, 
    fire_perc_ws_3ylegacy, fire_perc_ws_ppt_3ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_po4, cumulative_precip_mm, 
                    fire_perc_ws_3ylegacy, fire_perc_ws_ppt_3ylegacy))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 213
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
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
  0,0,0,"fire_perc_ws_ppt_3ylegacy_RS02"),4,12)

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
        file = "data_working/marss_fits/fit_080823_5state_po4_3ylegacy_mBFGS.rds")

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
# null.fit  96.4 12
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look good.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns. HO00 still a lag at 11?

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 3y legacy, 1 state ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_po4 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_po4, 
    cumulative_precip_mm, 
    fire_perc_ws_3ylegacy, fire_perc_ws_ppt_3ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_po4, cumulative_precip_mm, 
                    fire_perc_ws_3ylegacy, fire_perc_ws_ppt_3ylegacy))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 213
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
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
  0,0,0,"fire_perc_ws_ppt_3ylegacy_RS02"),4,12)

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
        file = "data_working/marss_fits/fit_080823_1state_po4_3ylegacy_mBFGS.rds")

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
# null.fit  66.3 18
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look good.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns. HO00 still a lag at 11?

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 4y legacy, 5 state ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_po4 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_po4, 
    cumulative_precip_mm, 
    fire_perc_ws_4ylegacy, fire_perc_ws_ppt_4ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_po4, cumulative_precip_mm, 
                    fire_perc_ws_4ylegacy, fire_perc_ws_ppt_4ylegacy))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 213
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
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
  0,0,0,"fire_perc_ws_ppt_4ylegacy_RS02"),4,12)

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
        file = "data_working/marss_fits/fit_080823_5state_po4_4ylegacy_mBFGS.rds")

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
# null.fit  94.9 12
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look fantastic, per usual.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 4y legacy, 1 state ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_po4 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_po4, 
    cumulative_precip_mm, 
    fire_perc_ws_4ylegacy, fire_perc_ws_ppt_4ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_po4, cumulative_precip_mm, 
                    fire_perc_ws_4ylegacy, fire_perc_ws_ppt_4ylegacy))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 213
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
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
  0,0,0,"fire_perc_ws_ppt_4ylegacy_RS02"),4,12)

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
        file = "data_working/marss_fits/fit_080823_1state_po4_4ylegacy_mBFGS.rds")

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
# null.fit  64.8 18
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look fantastic, per usual.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 5y legacy, 5 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_080823.rds")

# pivot wider for MARSS format
dat_po4 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_po4, 
    cumulative_precip_mm, 
    fire_perc_ws_5ylegacy, fire_perc_ws_ppt_5ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_po4, cumulative_precip_mm, 
                    fire_perc_ws_5ylegacy, fire_perc_ws_ppt_5ylegacy))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 213
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
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
  0,0,0,"fire_perc_ws_ppt_5ylegacy_RS02"),4,12)

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
        file = "data_working/marss_fits/fit_080823_5state_po4_5ylegacy_mBFGS.rds")

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
# null.fit 100.3 12
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look fantastic, per usual.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

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
dat_po4 <- dat %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    vwm_po4, 
    cumulative_precip_mm, 
    fire_perc_ws_5ylegacy, fire_perc_ws_ppt_5ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_po4, cumulative_precip_mm, 
                    fire_perc_ws_5ylegacy, fire_perc_ws_ppt_5ylegacy))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:5)
cov_cols = c(6:17)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 213
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
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
  0,0,0,"fire_perc_ws_ppt_5ylegacy_RS02"),4,12)

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
        file = "data_working/marss_fits/fit_080823_1state_po4_5ylegacy_mBFGS.rds")

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
# null.fit  64.7 18
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look fantastic, per usual.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

#### AIC Comparisons ####

# Compare all model fits for each legacy window to see which state 
# configuration was best.

# no legacy, 5 state
noleg_5state <- readRDS(file = "data_working/marss_fits/fit_080823_5state_po4_mBFGS.rds")

# no legacy, 1 state
noleg_1state <- readRDS(file = "data_working/marss_fits/fit_080823_1state_po4_mBFGS.rds")

bbmle::AICtab(noleg_5state, noleg_1state)

#              dAIC df
# noleg_1state  0.0 30
# noleg_5state 77.3 24

# 1y legacy, 5 state
leg1_5state <- readRDS(file = "data_working/marss_fits/fit_080823_5state_po4_1ylegacy_mBFGS.rds")

# 1y legacy, 1 state
leg1_1state <- readRDS(file = "data_working/marss_fits/fit_080823_1state_po4_1ylegacy_mBFGS.rds")

bbmle::AICtab(leg1_5state, leg1_1state)

#             dAIC df
# leg1_1state  0   30
# leg1_5state 82   24

# 2y legacy, 5 state
leg2_5state <- readRDS(file = "data_working/marss_fits/fit_080823_5state_po4_2ylegacy_mBFGS.rds")

# 2y legacy, 1 state
leg2_1state <- readRDS(file = "data_working/marss_fits/fit_080823_1state_po4_2ylegacy_mBFGS.rds")

bbmle::AICtab(leg2_5state, leg2_1state)

#             dAIC df
# leg2_1state  0.0 30
# leg2_5state 88.1 24

# 3y legacy, 5 state
leg3_5state <- readRDS(file = "data_working/marss_fits/fit_080823_5state_po4_3ylegacy_mBFGS.rds")

# 3y legacy, 1 state
leg3_1state <- readRDS(file = "data_working/marss_fits/fit_080823_1state_po4_3ylegacy_mBFGS.rds")

bbmle::AICtab(leg3_5state, leg3_1state)

#             dAIC df
# leg3_1state  0   30
# leg3_5state 90   24

# 4y legacy, 5 state
leg4_5state <- readRDS(file = "data_working/marss_fits/fit_080823_5state_po4_4ylegacy_mBFGS.rds")

# 4y legacy, 1 state
leg4_1state <- readRDS(file = "data_working/marss_fits/fit_080823_1state_po4_4ylegacy_mBFGS.rds")

bbmle::AICtab(leg4_5state, leg4_1state)

#             dAIC df
# leg4_1state  0   30
# leg4_5state 90   24

# 5y legacy, 5 state
leg5_5state <- readRDS(file = "data_working/marss_fits/fit_080823_5state_po4_5ylegacy_mBFGS.rds")

# 5y legacy, 1 state
leg5_1state <- readRDS(file = "data_working/marss_fits/fit_080823_1state_po4_5ylegacy_mBFGS.rds")

bbmle::AICtab(leg5_5state, leg5_1state)

#             dAIC df
# leg5_1state  0.0 30
# leg5_5state 84.5 24

# So, it would seem the 1 "state" model structure wins out every time.
# With little change between lag periods.

#### Results Figure ####

# End of script.