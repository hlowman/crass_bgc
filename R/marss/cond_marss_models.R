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
library(bbmle)
library(broom)

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
# saveRDS(fit, file = "data_working/marss_fits/fit_081423_8state_cond_mBFGS.rds")

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

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_081423.rds")

# select sites
# include these sites only (8 total - these have the longest 
# most complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws_1ylegacy, 
    fire_perc_ws_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_1ylegacy,
                    fire_perc_ws_ppt_1ylegacy)) 

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

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) #0
sum(is.na(dat_cond_log[,cov_cols])) #0
sum(is.infinite(dat_cond_log[,cov_cols])) #0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) #0
sum(is.na(dat_cov)) #0
sum(is.infinite(dat_cov)) #0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# check for cols with all zeros. this can cause model convergence issues
any(colSums(dat_cov)==0) # FALSE

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
  # fire_perc_ws_1ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_1ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 8, 24)

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
  V0="zero" ,
  tinitx=0
)

# Fit MARSS model 
# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, file = "data_working/marss_fits/fit_081423_8state_cond_1ylegacy_mBFGS.rds")

# DIAGNOSES

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

#          dAIC  df
# fit        0.0 48
# null.fit 208.2 24
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Fairly

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 1y legacy, 2 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_081423.rds")

# select sites
# include these sites only (8 total - these have the longest 
# most complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws_1ylegacy, 
    fire_perc_ws_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_1ylegacy,
                    fire_perc_ws_ppt_1ylegacy)) 

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

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) #0
sum(is.na(dat_cond_log[,cov_cols])) #0
sum(is.infinite(dat_cond_log[,cov_cols])) #0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) #0
sum(is.na(dat_cov)) #0
sum(is.infinite(dat_cov)) #0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# check for cols with all zeros. this can cause model convergence issues
any(colSums(dat_cov)==0) # FALSE

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
  # fire_perc_ws_1ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_1ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 8, 24)

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
  V0="zero" ,
  tinitx=0
)

# Fit MARSS model 
# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, file = "data_working/marss_fits/fit_081423_2state_cond_1ylegacy_mBFGS.rds")

# DIAGNOSES

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

#          dAIC  df
# fit        0.0 60
# null.fit 109.1 36
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No
# RED still looking a bit strange - but this is likely due to
# missingness

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Fairly

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 1y legacy, 1 state ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_081423.rds")

# select sites
# include these sites only (8 total - these have the longest 
# most complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws_1ylegacy, 
    fire_perc_ws_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_1ylegacy,
                    fire_perc_ws_ppt_1ylegacy)) 

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

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) #0
sum(is.na(dat_cond_log[,cov_cols])) #0
sum(is.infinite(dat_cond_log[,cov_cols])) #0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) #0
sum(is.na(dat_cov)) #0
sum(is.infinite(dat_cov)) #0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# check for cols with all zeros. this can cause model convergence issues
any(colSums(dat_cov)==0) # FALSE

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
  # fire_perc_ws_1ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_1ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 8, 24)

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
  V0="zero" ,
  tinitx=0
)

# Fit MARSS model 
# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, file = "data_working/marss_fits/fit_081423_1state_cond_1ylegacy_mBFGS.rds")

# DIAGNOSES

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

#          dAIC df
# fit        0  76
# null.fit 107  52
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No
# RED still looking a bit strange - but this is likely due to
# missingness

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Fairly

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 2y legacy, 8 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_081423.rds")

# select sites
# include these sites only (8 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws_2ylegacy, 
    fire_perc_ws_ppt_2ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_2ylegacy,
                    fire_perc_ws_ppt_2ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:9)
cov_cols = c(10:33)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) #0
sum(is.na(dat_cond_log[,resp_cols])) #595
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) #0
sum(is.na(dat_cond_log[,cov_cols])) #0
sum(is.infinite(dat_cond_log[,cov_cols])) #0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) #0
sum(is.na(dat_cov)) #0
sum(is.infinite(dat_cov)) #0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# check for cols with all zeros. this can cause model convergence issues
any(colSums(dat_cov)==0) # FALSE

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
  # fire_perc_ws_2ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_2ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 8, 24)

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
  V0="zero" ,
  tinitx=0
)

# Fit MARSS model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_081423_8state_cond_2ylegacy_mBFGS.rds")

# DIAGNOSES 
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

#          dAIC  df
# fit        0.0 48
# null.fit 208.7 24
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)?
# More so than in some diagnoses above.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 2y legacy, 2 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_081423.rds")

# select sites
# include these sites only (8 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws_2ylegacy, 
    fire_perc_ws_ppt_2ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_2ylegacy,
                    fire_perc_ws_ppt_2ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:9)
cov_cols = c(10:33)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) #0
sum(is.na(dat_cond_log[,resp_cols])) #595
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) #0
sum(is.na(dat_cond_log[,cov_cols])) #0
sum(is.infinite(dat_cond_log[,cov_cols])) #0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) #0
sum(is.na(dat_cov)) #0
sum(is.infinite(dat_cov)) #0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# check for cols with all zeros. this can cause model convergence issues
any(colSums(dat_cov)==0) # FALSE

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
  # fire_perc_ws_2ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_2ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 8, 24)

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
  V0="zero" ,
  tinitx=0
)

# Fit MARSS model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_081423_2state_cond_2ylegacy_mBFGS.rds")

# DIAGNOSES 
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

#          dAIC  df
# fit        0.0 60
# null.fit 123.7 36
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No (except RED)

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 2y legacy, 1 state ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_081423.rds")

# select sites
# include these sites only (8 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws_2ylegacy, 
    fire_perc_ws_ppt_2ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_2ylegacy,
                    fire_perc_ws_ppt_2ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:9)
cov_cols = c(10:33)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) #0
sum(is.na(dat_cond_log[,resp_cols])) #595
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) #0
sum(is.na(dat_cond_log[,cov_cols])) #0
sum(is.infinite(dat_cond_log[,cov_cols])) #0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) #0
sum(is.na(dat_cov)) #0
sum(is.infinite(dat_cov)) #0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# check for cols with all zeros. this can cause model convergence issues
any(colSums(dat_cov)==0) # FALSE

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
  # fire_perc_ws_2ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_2ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 8, 24)

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
  V0="zero" ,
  tinitx=0
)

# Fit MARSS model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_081423_1state_cond_2ylegacy_mBFGS.rds")

# DIAGNOSES 
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

#          dAIC  df
# fit        0.0 76
# null.fit 119.3 52
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No (except RED)

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 3y legacy, 8 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_081423.rds")

# select sites
# include these sites only (8 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws_3ylegacy, 
    fire_perc_ws_ppt_3ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_3ylegacy,
                    fire_perc_ws_ppt_3ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:9)
cov_cols = c(10:33)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) #0
sum(is.na(dat_cond_log[,resp_cols])) #595
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) #0
sum(is.na(dat_cond_log[,cov_cols])) #0
sum(is.infinite(dat_cond_log[,cov_cols])) #0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) #0
sum(is.na(dat_cov)) #0
sum(is.infinite(dat_cov)) #0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# check for cols with all zeros. this can cause model convergence issues
any(colSums(dat_cov)==0) # FALSE

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
  # fire_perc_ws_3ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_3ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 8, 24)

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
  V0="zero" ,
  tinitx=0
)

# Fit MARSS model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_081423_8state_cond_3ylegacy_mBFGS.rds")

# DIAGNOSES 
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

#          dAIC  df
# fit        0.0 48
# null.fit 205.8 24
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 3y legacy, 2 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_081423.rds")

# select sites
# include these sites only (8 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws_3ylegacy, 
    fire_perc_ws_ppt_3ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_3ylegacy,
                    fire_perc_ws_ppt_3ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:9)
cov_cols = c(10:33)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) #0
sum(is.na(dat_cond_log[,resp_cols])) #595
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) #0
sum(is.na(dat_cond_log[,cov_cols])) #0
sum(is.infinite(dat_cond_log[,cov_cols])) #0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) #0
sum(is.na(dat_cov)) #0
sum(is.infinite(dat_cov)) #0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# check for cols with all zeros. this can cause model convergence issues
any(colSums(dat_cov)==0) # FALSE

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
  # fire_perc_ws_3ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_3ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 8, 24)

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
  V0="zero" ,
  tinitx=0
)

# Fit MARSS model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_081423_2state_cond_3ylegacy_mBFGS.rds")

# DIAGNOSES 
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

#          dAIC  df
# fit        0.0 60
# null.fit 130.9 36
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Only at RED

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 3y legacy, 1 state ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_081423.rds")

# select sites
# include these sites only (8 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws_3ylegacy, 
    fire_perc_ws_ppt_3ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_3ylegacy,
                    fire_perc_ws_ppt_3ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:9)
cov_cols = c(10:33)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) #0
sum(is.na(dat_cond_log[,resp_cols])) #595
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) #0
sum(is.na(dat_cond_log[,cov_cols])) #0
sum(is.infinite(dat_cond_log[,cov_cols])) #0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) #0
sum(is.na(dat_cov)) #0
sum(is.infinite(dat_cov)) #0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# check for cols with all zeros. this can cause model convergence issues
any(colSums(dat_cov)==0) # FALSE

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
  # fire_perc_ws_3ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_3ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 8, 24)

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
  V0="zero" ,
  tinitx=0
)

# Fit MARSS model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_081423_1state_cond_3ylegacy_mBFGS.rds")

# DIAGNOSES 
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

#          dAIC  df
# fit        0.0 76
# null.fit 128.1 52
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Only at RED

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 4y legacy, 8 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_081423.rds")

# select sites
# include these sites only (8 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws_4ylegacy, 
    fire_perc_ws_ppt_4ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_4ylegacy,
                    fire_perc_ws_ppt_4ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:9)
cov_cols = c(10:33)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) #0
sum(is.na(dat_cond_log[,resp_cols])) #595
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) #0
sum(is.na(dat_cond_log[,cov_cols])) #0
sum(is.infinite(dat_cond_log[,cov_cols])) #0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) #0
sum(is.na(dat_cov)) #0
sum(is.infinite(dat_cov)) #0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# check for cols with all zeros. this can cause model convergence issues
any(colSums(dat_cov)==0) # FALSE

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
  # fire_perc_ws_4ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_4ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 8, 24)

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
  V0="zero" ,
  tinitx=0
)

# Fit MARSS model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_081423_8state_cond_4ylegacy_mBFGS.rds")

# DIAGNOSES 
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

#          dAIC df
# fit        0  48
# null.fit 209  24
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 4y legacy, 2 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_081423.rds")

# select sites
# include these sites only (8 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws_4ylegacy, 
    fire_perc_ws_ppt_4ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_4ylegacy,
                    fire_perc_ws_ppt_4ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:9)
cov_cols = c(10:33)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) #0
sum(is.na(dat_cond_log[,resp_cols])) #595
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) #0
sum(is.na(dat_cond_log[,cov_cols])) #0
sum(is.infinite(dat_cond_log[,cov_cols])) #0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) #0
sum(is.na(dat_cov)) #0
sum(is.infinite(dat_cov)) #0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# check for cols with all zeros. this can cause model convergence issues
any(colSums(dat_cov)==0) # FALSE

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
  # fire_perc_ws_4ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_4ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 8, 24)

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
  V0="zero" ,
  tinitx=0
)

# Fit MARSS model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_081423_2state_cond_4ylegacy_mBFGS.rds")

# DIAGNOSES 
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

#          dAIC df
# fit        0  60
# null.fit 131  36
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? RED is behaving strangely again

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Yes

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No but RED does have some questionable things here too.

#### 4y legacy, 1 state ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_081423.rds")

# select sites
# include these sites only (8 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws_4ylegacy, 
    fire_perc_ws_ppt_4ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_4ylegacy,
                    fire_perc_ws_ppt_4ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:9)
cov_cols = c(10:33)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) #0
sum(is.na(dat_cond_log[,resp_cols])) #595
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) #0
sum(is.na(dat_cond_log[,cov_cols])) #0
sum(is.infinite(dat_cond_log[,cov_cols])) #0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) #0
sum(is.na(dat_cov)) #0
sum(is.infinite(dat_cov)) #0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# check for cols with all zeros. this can cause model convergence issues
any(colSums(dat_cov)==0) # FALSE

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
  # fire_perc_ws_4ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_4ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 8, 24)

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
  V0="zero" ,
  tinitx=0
)

# Fit MARSS model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_081423_1state_cond_4ylegacy_mBFGS.rds")

# DIAGNOSES 
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

#          dAIC  df
# fit        0.0 76
# null.fit 126.7 52
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? RED is still behaving strangely

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Eh...

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? Lots of these look strange - perhaps this is due to the 1 state forcing.

#### 5y legacy, 8 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_081423.rds")

# select sites
# include these sites only (8 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws_5ylegacy, 
    fire_perc_ws_ppt_5ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_5ylegacy,
                    fire_perc_ws_ppt_5ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) #0
sum(is.na(dat_cond_log[,resp_cols])) #508
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) #0
sum(is.na(dat_cond_log[,cov_cols])) #0
sum(is.infinite(dat_cond_log[,cov_cols])) #0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) #0
sum(is.na(dat_cov)) #0
sum(is.infinite(dat_cov)) #0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# check for cols with all zeros. this can cause model convergence issues
any(colSums(dat_cov)==0) # FALSE

# make C matrix

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,0,0,
  0,0,0,"cumulative_precip_mm_RS02",0,0,0,
  0,0,0,0,"cumulative_precip_mm_RSAW",0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSA" ,0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_perc_ws_5ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_5ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 7, 21)

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
  V0="zero" ,
  tinitx=0
)

# Fit MARSS model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_081423_8state_cond_5ylegacy_mBFGS.rds")

# DIAGNOSES 
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

#          dAIC  df
# fit        0.0 42
# null.fit 212.1 21
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? For the most part

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No although some spikes at RED and RS02

#### 5y legacy, 2 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_081423.rds")

# select sites
# include these sites only (8 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws_5ylegacy, 
    fire_perc_ws_ppt_5ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_5ylegacy,
                    fire_perc_ws_ppt_5ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) #0
sum(is.na(dat_cond_log[,resp_cols])) #508
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) #0
sum(is.na(dat_cond_log[,cov_cols])) #0
sum(is.infinite(dat_cond_log[,cov_cols])) #0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) #0
sum(is.na(dat_cov)) #0
sum(is.infinite(dat_cov)) #0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# check for cols with all zeros. this can cause model convergence issues
any(colSums(dat_cov)==0) # FALSE

# make C matrix

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,0,0,
  0,0,0,"cumulative_precip_mm_RS02",0,0,0,
  0,0,0,0,"cumulative_precip_mm_RSAW",0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSA" ,0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_perc_ws_5ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_5ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 7, 21)

# make Q matrix
QQ <- matrix(list("s1","b12","b13","b14",0,    0,    0,    
                  "b12","s2","b23","b24",0,    0,    0,    
                  "b13","b23","s3","b34",0,    0,    0,    
                  "b14","b24","b34","s4",0,    0,    0,    
                  0,    0,    0,    0,"s5","b56","b57",
                  0,    0,    0,    0,"b56","s6","b67",
                  0,    0,    0,    0,"b57","b67","s7"), 7, 7)

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
  V0="zero" ,
  tinitx=0
)

# Fit MARSS model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_081423_2state_cond_5ylegacy_mBFGS.rds")

# DIAGNOSES 
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

#          dAIC  df
# fit        0.0 51
# null.fit 129.1 30
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No except RED

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? For the most part

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 5y legacy, 1 state ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_081423.rds")

# select sites
# include these sites only (8 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws_5ylegacy, 
    fire_perc_ws_ppt_5ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_5ylegacy,
                    fire_perc_ws_ppt_5ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) #0
sum(is.na(dat_cond_log[,resp_cols])) #508
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) #0
sum(is.na(dat_cond_log[,cov_cols])) #0
sum(is.infinite(dat_cond_log[,cov_cols])) #0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) #0
sum(is.na(dat_cov)) #0
sum(is.infinite(dat_cov)) #0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# check for cols with all zeros. this can cause model convergence issues
any(colSums(dat_cov)==0) # FALSE

# make C matrix

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,0,0,
  0,0,0,"cumulative_precip_mm_RS02",0,0,0,
  0,0,0,0,"cumulative_precip_mm_RSAW",0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSA" ,0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_perc_ws_5ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_5ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 7, 21)

# make Q matrix
QQ <- matrix(list("s1","b12","b13","b14","b15","b16","b17",
                  "b12","s2","b23","b24","b25","b26","b27",
                  "b13","b23","s3","b34","b35","b36","b37",
                  "b14","b24","b34","s4","b45","b46","b47",
                  "b15","b25","b35","b45","s5","b56","b57",
                  "b16","b26","b36","b46","b56","s6","b67",
                  "b17","b27","b37","b47","b57","b67","s7"), 7, 7)

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
  V0="zero" ,
  tinitx=0
)

# Fit MARSS model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_fits/fit_081423_1state_cond_5ylegacy_mBFGS.rds")

# DIAGNOSES 
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

#          dAIC  df
# fit        0.0 63
# null.fit 141.8 42
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): All equal zero? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? For the most part

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? Hmmm these aren't looking *great*

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

broom::glance(noleg_8state) # AICc 1808.588
stats4::BIC(noleg_8state) # BIC 2030.677

broom::glance(noleg_2state) # AICc 1742.701
stats4::BIC(noleg_2state) # BIC 2018.384

broom::glance(noleg_1state) # AICc 1753.883
stats4::BIC(noleg_1state) # BIC 2099.706

# So, based on lowest AIC, AICc, and BIC values, 2-state wins.

###

# 1y legacy, 8 state
leg1_8state <- readRDS(file = "data_working/marss_fits/fit_081423_8state_cond_1ylegacy_mBFGS.rds")

# 1y legacy, 2 state
leg1_2state <- readRDS(file = "data_working/marss_fits/fit_081423_2state_cond_1ylegacy_mBFGS.rds")

# 1y legacy, 1 state
leg1_1state <- readRDS(file = "data_working/marss_fits/fit_081423_1state_cond_1ylegacy_mBFGS.rds")

bbmle::AICtab(leg1_8state, leg1_2state, leg1_1state)

#             dAIC df
# leg1_2state  0.0 60
# leg1_1state  9.0 76
# leg1_8state 74.2 48

broom::glance(leg1_8state) # AICc 1830.287
stats4::BIC(leg1_8state) # BIC 2052.376

broom::glance(leg1_2state) # AICc 1759.479 
stats4::BIC(leg1_2state) # BIC 2035.162

broom::glance(leg1_1state) # AICc 1774.31
stats4::BIC(leg1_1state) # BIC 2120.133

# So, based on lowest AIC, AICc, and BIC values, 2-state wins.

###

# 2y legacy, 8 state
leg2_8state <- readRDS(file = "data_working/marss_fits/fit_081423_8state_cond_2ylegacy_mBFGS.rds")

# 2y legacy, 2 state
leg2_2state <- readRDS(file = "data_working/marss_fits/fit_081423_2state_cond_2ylegacy_mBFGS.rds")

# 2y legacy, 1 state
leg2_1state <- readRDS(file = "data_working/marss_fits/fit_081423_1state_cond_2ylegacy_mBFGS.rds")

bbmle::AICtab(leg2_8state, leg2_2state, leg2_1state)

#             dAIC df
# leg2_2state  0.0 60
# leg2_1state 11.2 76
# leg2_8state 88.3 48

broom::glance(leg2_8state) # AICc 1829.803
stats4::BIC(leg2_8state) # BIC 2051.893

broom::glance(leg2_2state) # AICc 1744.909
stats4::BIC(leg2_2state) # BIC 2020.592

broom::glance(leg2_1state) # AICc 1761.959
stats4::BIC(leg2_1state) # BIC 2107.782

# So, based on lowest AIC, AICc, and BIC values, 2-state wins.

###

# 3y legacy, 8 state
leg3_8state <- readRDS(file = "data_working/marss_fits/fit_081423_8state_cond_3ylegacy_mBFGS.rds")

# 3y legacy, 2 state
leg3_2state <- readRDS(file = "data_working/marss_fits/fit_081423_2state_cond_3ylegacy_mBFGS.rds")

# 3y legacy, 1 state
leg3_1state <- readRDS(file = "data_working/marss_fits/fit_081423_1state_cond_3ylegacy_mBFGS.rds")

bbmle::AICtab(leg3_8state, leg3_2state, leg3_1state)

#             dAIC df
# leg3_2state  0.0 60
# leg3_1state  9.5 76
# leg3_8state 98.4 48

broom::glance(leg3_8state) # AICc 1832.737
stats4::BIC(leg3_8state) # BIC 2054.827

broom::glance(leg3_2state) # AICc 1737.77 
stats4::BIC(leg3_2state) # BIC 2013.453

broom::glance(leg3_1state) # AICc 1753.133
stats4::BIC(leg3_1state) # BIC 2098.956

# So, based on lowest AIC, AICc, and BIC values, 2-state wins.

###

# 4y legacy, 8 state
leg4_8state <- readRDS(file = "data_working/marss_fits/fit_081423_8state_cond_4ylegacy_mBFGS.rds")

# 4y legacy, 2 state
leg4_2state <- readRDS(file = "data_working/marss_fits/fit_081423_2state_cond_4ylegacy_mBFGS.rds")

# 4y legacy, 1 state
leg4_1state <- readRDS(file = "data_working/marss_fits/fit_081423_1state_cond_4ylegacy_mBFGS.rds")

bbmle::AICtab(leg4_8state, leg4_2state, leg4_1state)

#             dAIC df
# leg4_2state  0.0 60
# leg4_1state 11.0 76
# leg4_8state 95.2 48

broom::glance(leg4_8state) # AICc 1829.456
stats4::BIC(leg4_8state) # BIC 2051.546

broom::glance(leg4_2state) # AICc 1737.66
stats4::BIC(leg4_2state) # BIC 2013.343

broom::glance(leg4_1state) # AICc 1754.531
stats4::BIC(leg4_1state) # BIC 2100.354

# So, based on lowest AIC, AICc, and BIC values, 2-state wins.

###

# 5y legacy, 8 state
leg5_8state <- readRDS(file = "data_working/marss_fits/fit_081423_8state_cond_5ylegacy_mBFGS.rds")

# 5y legacy, 2 state
leg5_2state <- readRDS(file = "data_working/marss_fits/fit_081423_2state_cond_5ylegacy_mBFGS.rds")

# 5y legacy, 1 state
leg5_1state <- readRDS(file = "data_working/marss_fits/fit_081423_1state_cond_5ylegacy_mBFGS.rds")

bbmle::AICtab(leg5_8state, leg5_2state, leg5_1state)

#             dAIC df
# leg5_1state  0.0 63
# leg5_2state  2.0 51
# leg5_8state 63.8 42

broom::glance(leg5_8state) # AICc 1726.657
stats4::BIC(leg5_8state) # BIC 1916.156

broom::glance(leg5_2state) # AICc 1667.376
stats4::BIC(leg5_2state) # BIC 1896.106

broom::glance(leg5_1state) # AICc 1669.438
stats4::BIC(leg5_1state) # BIC 1949.651

# So, based on lowest AIC value, 1-state wins, but based on lowest AICc and BIC values, 2-state wins.

###

# So, it would seem the 2 "state" model structure wins out.

stats4::BIC(noleg_2state) # BIC 2018.384
stats4::BIC(leg1_2state) # BIC 2035.162
stats4::BIC(leg2_2state) # BIC 2020.592
stats4::BIC(leg3_2state) # BIC 2013.453
stats4::BIC(leg4_2state) # BIC 2013.343
stats4::BIC(leg5_2state) # BIC 1896.106

# And when comparing all models, the 5 year window/lag is most parsimonious.

#### Results Figure ####

# For presentation consistency, I will only be creating figures with a
# single state configuration, whichever yielded the most parsimonious
# models. So, in this case, all conductivity figures will 
# represent the "2 state" scenario.

# Extract necessary confidence interval info
noleg_est <- MARSSparamCIs(noleg_2state)
leg1y_est <- MARSSparamCIs(leg1_2state)
leg2y_est <- MARSSparamCIs(leg2_2state)
leg3y_est <- MARSSparamCIs(leg3_2state)
leg4y_est <- MARSSparamCIs(leg4_2state)
leg5y_est <- MARSSparamCIs(leg5_2state)

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
CIs$Parm_simple = c(rep("Ppt",8),
                    rep("Perc. burn",8),
                    rep("Ppt x Perc. burn",8),
                    
                    rep("Ppt",8),
                    rep("Perc. burn",8),
                    rep("Ppt x Perc. burn",8),
                    
                    rep("Ppt",8),
                    rep("Perc. burn",8),
                    rep("Ppt x Perc. burn",8),
                    
                    rep("Ppt",8),
                    rep("Perc. burn",8),
                    rep("Ppt x Perc. burn",8),
                    
                    rep("Ppt",8),
                    rep("Perc. burn",8),
                    rep("Ppt x Perc. burn",8),
                    
                    rep("Ppt",7),
                    rep("Perc. burn",7),
                    rep("Ppt x Perc. burn",7))

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
  mutate(region = case_when(Stream %in% c("AB00", "GV01", 
                                          "HO00", "RS02") ~ "Santa Barbara",
                            TRUE ~ "Valles Caldera")) %>%
  mutate(site = factor(Stream,
                       levels = c("AB00", "GV01", "HO00", "RS02",
                                  "EFJ", "RED", "RSA", "RSAW")))

my_palette <- c("black", "white", "black")

# Plot results
(SpCond_fig <- ggplot(CIs, aes(x = factor(Parm_simple, 
                                       levels = c("Ppt x Perc. burn",
                                                  "Perc. burn",
                                                  "Ppt")),
                            y = Est., fill = sig, shape = site)) + 
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
                  position=position_dodge(width = 0.5), width = 0) +
    geom_point(position=position_dodge(width = 0.5), 
               alpha = 0.8, size = 8) + 
    scale_shape_manual(values = c(21, 22, 23, 24,
                                  21, 22, 23, 24)) +
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
    facet_grid(region~Model))

# Export plot.
# ggsave(("MARSS_SpCon_081423.png"),
#        path = "figures",
#        width = 65,
#        height = 24,
#        units = "cm"
# )

# End of script.
