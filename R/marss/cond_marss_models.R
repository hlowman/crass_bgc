# Conductivity MARSS models
# with fire x ppt interactions and legacy effects
# as well as multiple state structure for CA/NM sites
# Script started August 14, 2023 by Heili Lowman

# This script will run 18 Conductivity MARSS models.
# Note, each model fit will remove all stored data, until you reach the AIC
# portion of this script.

#### Setup ####

# Load packages.
library(tidyverse)
library(lubridate)
library(MARSS)
library(naniar) 
library(here)

#### 0y legacy, 8 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat <- readRDS("data_working/marss_data_sb_vc_081423.rds")

# select sites
# include these sites only (8 total - these have the longest most complete ts 
# for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site) # Yay! All the timeseries are the same length.

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws, fire_perc_ws_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws, fire_perc_ws_ppt)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:9)
cov_cols = c(10:33)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols]))
sum(is.na(dat_cond_log[,resp_cols])) # 595
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) # 0
sum(is.na(dat_cond_log[,cov_cols])) # 0
sum(is.infinite(dat_cond_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]

# check for cols with all zeros
any(colSums(dat_cov)==0)

# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)

# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
#no

# make C matrix 

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,0,0,0,
  0,0,0,"cumulative_precip_mm_RS02",0,0,0,0,
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,0,
  0,0,0,0,0,0,"cumulative_precip_mm_RSA" ,0,
  0,0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_perc_ws: fire effect in 2 m window
  "fire_perc_ws_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_RED" ,
  # fire_perc_ws_ppt: interaction of cum. ppt with fire % burn in 6 m window
  "fire_perc_ws_ppt_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_RED" ), 8, 24)

# Model setup for MARSS 

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
  tinitx=0,
  V0="zero"
)

# Fit MARSS model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, file = "data_working/marss_fits/fit_081423_8state_cond_mBFGS.rds")

# Diagnoses
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
  #x0 = matrix(x0_fixed),
  tinitx=0,
  V0="zero"
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"
null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#          dAIC  df
# fit        0.0 48
# null.fit 229.9 24
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight line)? Yes for the most part

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 0y legacy, 2 states ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat <- readRDS("data_working/marss_data_sb_vc_081423.rds")

# select sites
# include these sites only (8 total - these have the longest most complete ts 
# for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site) # Yay! All the timeseries are the same length.

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws, fire_perc_ws_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws, fire_perc_ws_ppt)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:9)
cov_cols = c(10:33)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 595
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) # 0
sum(is.na(dat_cond_log[,cov_cols])) # 0
sum(is.infinite(dat_cond_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]

# check for cols with all zeros
any(colSums(dat_cov)==0)

# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)

# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
#no

# make C matrix 

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,0,0,0,
  0,0,0,"cumulative_precip_mm_RS02",0,0,0,0,
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,0,
  0,0,0,0,0,0,"cumulative_precip_mm_RSA" ,0,
  0,0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_perc_ws: fire effect in 2 m window
  "fire_perc_ws_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_RED" ,
  # fire_perc_ws_ppt: interaction of cum. ppt with fire % burn in 6 m window
  "fire_perc_ws_ppt_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_RED" ), 8, 24)

# make Q matrix
QQ <- matrix(list("s1","b12","b13","b14",0,    0,    0,    0,
                  "b12","s2","b23","b24",0,    0,    0,    0,
                  "b13","b23","s3","b34",0,    0,    0,    0,
                  "b14","b24","b34","s4",0,    0,    0,    0,
                  0,    0,    0,    0,"s5","b56","b57","b58",
                  0,    0,    0,    0,"b56","s6","b67","b68",
                  0,    0,    0,    0,"b57","b67","s7","b78",
                  0,    0,    0,    0,"b58","b68","b78","s8"), 8, 8)

# Model setup for MARSS 

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
  tinitx=0,
  V0="zero"
)

# Fit MARSS model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, file = "data_working/marss_fits/fit_081423_2state_cond_mBFGS.rds")

# Diagnoses
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
  #x0 = matrix(x0_fixed),
  tinitx=0,
  V0="zero"
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"
null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#          dAIC  df
# fit        0.0 60
# null.fit 125.9 36
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight line)? Yes for the most part - RED looking a bit strange

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 0y legacy, 1 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat <- readRDS("data_working/marss_data_sb_vc_081423.rds")

# select sites
# include these sites only (8 total - these have the longest most complete ts 
# for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site) # Yay! All the timeseries are the same length.

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws, fire_perc_ws_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws, fire_perc_ws_ppt)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:9)
cov_cols = c(10:33)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 595
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) # 0
sum(is.na(dat_cond_log[,cov_cols])) # 0
sum(is.infinite(dat_cond_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]

# check for cols with all zeros
any(colSums(dat_cov)==0)

# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)

# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
#no

# make C matrix 

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,0,0,0,
  0,0,0,"cumulative_precip_mm_RS02",0,0,0,0,
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,0,
  0,0,0,0,0,0,"cumulative_precip_mm_RSA" ,0,
  0,0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_perc_ws: fire effect in 2 m window
  "fire_perc_ws_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_RED" ,
  # fire_perc_ws_ppt: interaction of cum. ppt with fire % burn in 6 m window
  "fire_perc_ws_ppt_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_RED" ), 8, 24)

# make Q matrix
QQ <- matrix(list("s1","b12","b13","b14","b15","b16","b17","b18",
                  "b12","s2","b23","b24","b25","b26","b27","b28",
                  "b13","b23","s3","b34","b35","b36","b37","b38",
                  "b14","b24","b34","s4","b45","b46","b47","b48",
                  "b15","b25","b35","b45","s5","b56","b57","b58",
                  "b16","b26","b36","b46","b56","s6","b67","b68",
                  "b17","b27","b37","b47","b57","b67","s7","b78",
                  "b18","b28","b38","b48","b58","b68","b78","s8"), 8, 8)

# Model setup for MARSS 

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
  tinitx=0,
  V0="zero"
)

# Fit MARSS model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, file = "data_working/marss_fits/fit_081423_1state_cond_mBFGS.rds")

# Diagnoses
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
  #x0 = matrix(x0_fixed),
  tinitx=0,
  V0="zero"
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"
null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#          dAIC  df
# fit        0.0 76
# null.fit 127.4 52
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight line)? Yes for the most part although looking a bit strange

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

### Overall ###
# None of these diagnoses look prohibitively bad.

#### 1y legacy, 8 states ####

#### 1y legacy, 2 states ####

#### 1y legacy, 1 state ####

#### 2y legacy, 8 states ####

#### 2y legacy, 2 states ####

#### 2y legacy, 1 state ####

#### 3y legacy, 8 states ####

#### 3y legacy, 2 states ####

#### 3y legacy, 1 state ####

#### 4y legacy, 8 states ####

#### 4y legacy, 2 states ####

#### 4y legacy, 1 state ####

#### 5y legacy, 8 states ####

#### 5y legacy, 2 states ####

#### 5y legacy, 1 state ####

#### AIC Comparisons ####

# Compare all model fits for each legacy window to see which state 
# configuration was best.

# no legacy, 8 state
noleg_8state <- readRDS(file = "data_working/marss_fits/fit_081423_8state_cond_mBFGS.rds")

# no legacy, 2 state
noleg_2state <- readRDS(file = "data_working/marss_fits/fit_081423_2state_cond_mBFGS.rds")

# no legacy, 1 state
noleg_1state <- readRDS(file = "data_working/marss_fits/fit_081423_1state_cond_mBFGS.rds")

bbmle::AICtab(noleg_8state, noleg_2state, noleg_1state)

#              dAIC df
# noleg_2state  0.0 60
# noleg_1state  5.3 76
# noleg_8state 69.3 48

# 1y legacy, 4 state
leg1_4state <- readRDS(file = "data_working/marss_fits/fit_081023_4state_nh4_1ylegacy_mBFGS.rds")

# 1y legacy, 1 state
leg1_1state <- readRDS(file = "data_working/marss_fits/fit_081023_1state_nh4_1ylegacy_mBFGS.rds")

bbmle::AICtab(leg1_4state, leg1_1state)

#             dAIC df
# leg1_1state  0.0 30
# leg1_4state  2.7 24

# 2y legacy, 4 state
leg2_4state <- readRDS(file = "data_working/marss_fits/fit_081023_4state_nh4_2ylegacy_mBFGS.rds")

# 2y legacy, 1 state
leg2_1state <- readRDS(file = "data_working/marss_fits/fit_081023_1state_nh4_2ylegacy_mBFGS.rds")

bbmle::AICtab(leg2_4state, leg2_1state)

#             dAIC df
# leg2_1state  0.0 30
# leg2_4state 10.2 24

# 3y legacy, 4 state
leg3_4state <- readRDS(file = "data_working/marss_fits/fit_081023_4state_nh4_3ylegacy_mBFGS.rds")

# 3y legacy, 1 state
leg3_1state <- readRDS(file = "data_working/marss_fits/fit_081023_1state_nh4_3ylegacy_mBFGS.rds")

bbmle::AICtab(leg3_4state, leg3_1state)

#             dAIC df
# leg3_1state  0.0 30
# leg3_4state 10.9 24

# 4y legacy, 4 state
leg4_4state <- readRDS(file = "data_working/marss_fits/fit_081023_4state_nh4_4ylegacy_mBFGS.rds")

# 4y legacy, 1 state
leg4_1state <- readRDS(file = "data_working/marss_fits/fit_081023_1state_nh4_4ylegacy_mBFGS.rds")

bbmle::AICtab(leg4_4state, leg4_1state)

#             dAIC df
# leg4_1state  0.0 30
# leg4_4state  9.9 24

# 5y legacy, 4 state
leg5_4state <- readRDS(file = "data_working/marss_fits/fit_081023_4state_nh4_5ylegacy_mBFGS.rds")

# 5y legacy, 1 state
leg5_1state <- readRDS(file = "data_working/marss_fits/fit_081023_1state_nh4_5ylegacy_mBFGS.rds")

bbmle::AICtab(leg5_4state, leg5_1state)

#             dAIC df
# leg5_1state  0.0 30
# leg5_4state  5.5 24

# So, it would seem the 1 "state" model structure  wins out, and 
# increasingly so with increasingly lag windows.

bbmle::AICtab(noleg_1state, leg1_1state, leg2_1state,
              leg3_1state, leg4_1state, leg5_1state)

#              dAIC df
# noleg_1state  0.0 30
# leg3_1state   2.8 30
# leg4_1state   3.1 30
# leg1_1state   6.4 30
# leg5_1state  12.8 30
# leg2_1state  14.7 30

# And when comparing all models, the 0 year window/lag is most parsimonious.

#### Results Figure ####

# For presentation consistency, I will only be creating figures with a
# single state configuration, whichever yielded the most parsimonious
# models. So, in this case, all NH4 figures will represent the "1 state"
# scenario.

# Extract necessary confidence interval info
noleg_est <- MARSSparamCIs(noleg_1state)
leg1y_est <- MARSSparamCIs(leg1_1state)
leg2y_est <- MARSSparamCIs(leg2_1state)
leg3y_est <- MARSSparamCIs(leg3_1state)
leg4y_est <- MARSSparamCIs(leg4_1state)
leg5y_est <- MARSSparamCIs(leg5_1state)

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
    facet_grid(.~Model))

# Export plot.
# ggsave(("MARSS_NH4_081023.png"),
#        path = "figures",
#        width = 65,
#        height = 12,
#        units = "cm"
# )

# End of script.
