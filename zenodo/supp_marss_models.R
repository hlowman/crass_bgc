# Conductivity MARSS models
# with fire x ppt interactions and legacy effects
# as well as multiple state structure for CA/NM sites
# Script started March 11,2024 by Heili Lowman

# This script will run 18 Conductivity MARSS models to address
# reviewer comments regarding using a binary variable in place
# of the percent of the watershed burned.

# So, the first part of this script will convert the existing
# dataset to include only 0/1s for the fire covariate and 
# recalculate the interaction term.

# Note, each model fit will remove all stored data, until you reach
# the AIC portion of this script.

#### Setup ####

# Load packages.
library(tidyverse)
library(lubridate)
library(MARSS)
library(naniar) 
library(here)
library(bbmle)
library(broom)

#### Binary Covariate ####

# load data with fire x ppt interactions and legacy effects
dat <- readRDS("data_working/marss_data_sb_vc_091123.rds")

# using existing binary covariate data ("pa" or presence/absence)
dat_ed <- dat %>%
  mutate(fire_pa_ppt = fire_pa_0.5ylegacy*cumulative_precip_mm,
         fire_pa_ppt_1ylegacy = fire_pa_1ylegacy*cumulative_precip_mm,
         fire_pa_ppt_2ylegacy = fire_pa_2ylegacy*cumulative_precip_mm,
         fire_pa_ppt_3ylegacy = fire_pa_3ylegacy*cumulative_precip_mm,
         fire_pa_ppt_4ylegacy = fire_pa_4ylegacy*cumulative_precip_mm,
         fire_pa_ppt_5ylegacy = fire_pa_5ylegacy*cumulative_precip_mm)

saveRDS(dat_ed, "data_working/marss_data_sb_vc_031124.rds")

#### 0y legacy, 7 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat <- readRDS("data_working/marss_data_sb_vc_031124.rds")

# select sites
# include these sites only (7 total - these have the longest most complete ts 
# for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site) # Yay! All the timeseries are the same length.

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa, fire_pa_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa, fire_pa_ppt)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols]))
sum(is.na(dat_cond_log[,resp_cols])) # 479
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
  "cumulative_precip_mm_AB00",0,0,0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,0,0,
  0,0,0,"cumulative_precip_mm_RS02",0,0,0,
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_pa: fire effect in 2 m window
  "fire_pa_AB00",0,0,0,0,0,0,
  0,"fire_pa_GV01",0,0,0,0,0,
  0,0,"fire_pa_HO00",0,0,0,0,
  0,0,0,"fire_pa_RS02",0,0,0,
  0,0,0,0,"fire_pa_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_RSAW",0,
  0,0,0,0,0,0,"fire_pa_RED" ,
  # fire_pa_ppt: interaction of cum. ppt with fire in 6 m window
  "fire_pa_ppt_AB00",0,0,0,0,0,0,
  0,"fire_pa_ppt_GV01",0,0,0,0,0,
  0,0,"fire_pa_ppt_HO00",0,0,0,0,
  0,0,0,"fire_pa_ppt_RS02",0,0,0,
  0,0,0,0,"fire_pa_ppt_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_ppt_RSAW",0,
  0,0,0,0,0,0,"fire_pa_ppt_RED" ), 7, 21)

# Make R matrix

RR <- matrix(list("v1",0,0,0,0,0,0,
                  0,"v1",0,0,0,0,0,
                  0,0,"v1",0,0,0,0,
                  0,0,0,"v1",0,0,0,
                  0,0,0,0,"v2",0,0,
                  0,0,0,0,0,"v2",0,
                  0,0,0,0,0,0,"v2"), 7, 7)

# Model setup for MARSS 

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", # 7 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = RR, # obs. error
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
saveRDS(fit, file = "data_working/marss_fits/fit_031124_7state_cond_binary_mBFGS.rds")

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
  R = RR, 
  ### initial conditions ###
  #x0 = matrix(x0_fixed),
  tinitx=0,
  V0="zero"
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"
null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit) # AICc 1616.778 
MARSSaic(null.fit) # AICc 1806.168   

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight line)? Yes for the most part

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 0y legacy, 2 states ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat <- readRDS("data_working/marss_data_sb_vc_031124.rds")

# select sites
# include these sites only (7 total - these have the longest most complete ts 
# for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site) # Yay! All the timeseries are the same length.

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa, fire_pa_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa, fire_pa_ppt)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 479
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
  "cumulative_precip_mm_AB00",0,0,0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,0,0,
  0,0,0,"cumulative_precip_mm_RS02",0,0,0,
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_pa: fire effect in 2 m window
  "fire_pa_AB00",0,0,0,0,0,0,
  0,"fire_pa_GV01",0,0,0,0,0,
  0,0,"fire_pa_HO00",0,0,0,0,
  0,0,0,"fire_pa_RS02",0,0,0,
  0,0,0,0,"fire_pa_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_RSAW",0,
  0,0,0,0,0,0,"fire_pa_RED" ,
  # fire_pa_ppt: interaction of cum. ppt with fire in 6 m window
  "fire_pa_ppt_AB00",0,0,0,0,0,0,
  0,"fire_pa_ppt_GV01",0,0,0,0,0,
  0,0,"fire_pa_ppt_HO00",0,0,0,0,
  0,0,0,"fire_pa_ppt_RS02",0,0,0,
  0,0,0,0,"fire_pa_ppt_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_ppt_RSAW",0,
  0,0,0,0,0,0,"fire_pa_ppt_RED" ), 7, 21)

# make Q matrix
QQ <- matrix(list("s1","b12","b13","b14",0,    0,    0,
                  "b12","s2","b23","b24",0,    0,    0,
                  "b13","b23","s3","b34",0,    0,    0,
                  "b14","b24","b34","s4",0,    0,    0,
                  0,    0,    0,    0,"s5","b56","b57",
                  0,    0,    0,    0,"b56","s6","b67",
                  0,    0,    0,    0,"b57","b67","s7"), 7, 7)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,0,
                  0,"v1",0,0,0,0,0,
                  0,0,"v1",0,0,0,0,
                  0,0,0,"v1",0,0,0,
                  0,0,0,0,"v2",0,0,
                  0,0,0,0,0,"v2",0,
                  0,0,0,0,0,0,"v2"), 7, 7)

# Model setup for MARSS 

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = QQ, # 2 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = RR, 
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
saveRDS(fit, file = "data_working/marss_fits/fit_031124_2state_cond_binary_mBFGS.rds")

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
  R = RR, 
  ### initial conditions ###
  #x0 = matrix(x0_fixed),
  tinitx=0,
  V0="zero"
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"
null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit) # AICc 1592.289   
MARSSaic(null.fit) # AICc 1666.072   

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight line)? Yes for the most part

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 0y legacy, 1 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat <- readRDS("data_working/marss_data_sb_vc_031124.rds")

# select sites
# include these sites only (7 total - these have the longest most complete ts 
# for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site) # Yay! All the timeseries are the same length.

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa, fire_pa_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa, fire_pa_ppt)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 479
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

# make C matrix 
CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,0,0,
  0,0,0,"cumulative_precip_mm_RS02",0,0,0,
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED",
  # fire_pa: fire effect in 2 m window
  "fire_pa_AB00",0,0,0,0,0,0,
  0,"fire_pa_GV01",0,0,0,0,0,
  0,0,"fire_pa_HO00",0,0,0,0,
  0,0,0,"fire_pa_RS02",0,0,0,
  0,0,0,0,"fire_pa_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_RSAW",0,
  0,0,0,0,0,0,"fire_pa_RED",
  # fire_pa_ppt: interaction of cum. ppt with fire in 6 m window
  "fire_pa_ppt_AB00",0,0,0,0,0,0,
  0,"fire_pa_ppt_GV01",0,0,0,0,0,
  0,0,"fire_pa_ppt_HO00",0,0,0,0,
  0,0,0,"fire_pa_ppt_RS02",0,0,0,
  0,0,0,0,"fire_pa_ppt_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_ppt_RSAW",0,
  0,0,0,0,0,0,"fire_pa_ppt_RED"), 7, 21)

# make Q matrix
QQ <- matrix(list("s1","b12","b13","b14","b15","b16","b17",
                  "b12","s2","b23","b24","b25","b26","b27",
                  "b13","b23","s3","b34","b35","b36","b37",
                  "b14","b24","b34","s4","b45","b46","b47",
                  "b15","b25","b35","b45","s5","b56","b57",
                  "b16","b26","b36","b46","b56","s6","b67",
                  "b17","b27","b37","b47","b57","b67","s7"), 7, 7)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,0,
                  0,"v1",0,0,0,0,0,
                  0,0,"v1",0,0,0,0,
                  0,0,0,"v1",0,0,0,
                  0,0,0,0,"v2",0,0,
                  0,0,0,0,0,"v2",0,
                  0,0,0,0,0,0,"v2"), 7, 7)

# Model setup for MARSS 

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
  R = RR, 
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
saveRDS(fit, file = "data_working/marss_fits/fit_031124_1state_cond_binary_mBFGS.rds")

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
  R = RR, 
  ### initial conditions ###
  #x0 = matrix(x0_fixed),
  tinitx=0,
  V0="zero"
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"
null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit) # AICc 1583.221   
MARSSaic(null.fit) # AICc 1653.978 

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight line)? Yes for the most part, VC sites don't look fantastic

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? Does not appear so

#### 1y legacy, 7 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_031124.rds")

# select sites
# include these sites only (7 total - these have the longest 
# most complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa_1ylegacy, 
    fire_pa_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa_1ylegacy,
                    fire_pa_ppt_1ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 479
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
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
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
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_pa_1ylegacy
  "fire_pa_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_1ylegacy_RED" ,
  # fire_pa_ppt_1ylegacy
  "fire_pa_ppt_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_ppt_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_ppt_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_ppt_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_ppt_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_ppt_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_ppt_1ylegacy_RED" ), 7, 21)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,0,
                  0,"v1",0,0,0,0,0,
                  0,0,"v1",0,0,0,0,
                  0,0,0,"v1",0,0,0,
                  0,0,0,0,"v2",0,0,
                  0,0,0,0,0,"v2",0,
                  0,0,0,0,0,0,"v2"), 7, 7)

# Model setup for MARSS

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC,
  c = dat_cov,
  Q = "diagonal and unequal", # 7 states
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = RR, # obs. error
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
saveRDS(fit, file = "data_working/marss_fits/fit_031124_7state_cond_1ylegacy_binary_mBFGS.rds")

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
  R = RR, 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit) # AICc 1618.708   
MARSSaic(null.fit) # AICc 1806.168

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Fairly

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 1y legacy, 2 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_031124.rds")

# select sites
# include these sites only (7 total - these have the longest 
# most complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa_1ylegacy, 
    fire_pa_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa_1ylegacy,
                    fire_pa_ppt_1ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 479
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) # 0
sum(is.na(dat_cond_log[,cov_cols])) # 0
sum(is.infinite(dat_cond_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
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
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_pa_1ylegacy
  "fire_pa_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_1ylegacy_RED" ,
  # fire_pa_ppt_1ylegacy
  "fire_pa_ppt_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_ppt_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_ppt_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_ppt_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_ppt_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_ppt_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_ppt_1ylegacy_RED" ), 7, 21)

# make Q matrix
QQ <- matrix(list("s1","b12","b13","b14",0,    0,    0,
                  "b12","s2","b23","b24",0,    0,    0,
                  "b13","b23","s3","b34",0,    0,    0,
                  "b14","b24","b34","s4",0,    0,    0,
                  0,    0,    0,    0,"s5","b56","b57",
                  0,    0,    0,    0,"b56","s6","b67",
                  0,    0,    0,    0,"b57","b67","s7"), 7, 7)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,0,
                  0,"v1",0,0,0,0,0,
                  0,0,"v1",0,0,0,0,
                  0,0,0,"v1",0,0,0,
                  0,0,0,0,"v2",0,0,
                  0,0,0,0,0,"v2",0,
                  0,0,0,0,0,0,"v2"), 7, 7)

# Model setup for MARSS

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = QQ, # 2 states
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = RR, # obs. error
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
saveRDS(fit, file = "data_working/marss_fits/fit_031124_2state_cond_1ylegacy_binary_mBFGS.rds")

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
  R = RR, 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit) # AICc 1587.028  
MARSSaic(null.fit) # AICc 1666.072 

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Fairly

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 1y legacy, 1 state ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_031124.rds")

# select sites
# include these sites only (7 total - these have the longest 
# most complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa_1ylegacy, 
    fire_pa_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa_1ylegacy,
                    fire_pa_ppt_1ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 479
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) # 0
sum(is.na(dat_cond_log[,cov_cols])) # 0
sum(is.infinite(dat_cond_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
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
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_pa_1ylegacy
  "fire_pa_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_1ylegacy_RED" ,
  # fire_pa_ppt_1ylegacy
  "fire_pa_ppt_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_ppt_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_ppt_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_ppt_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_ppt_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_ppt_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_ppt_1ylegacy_RED" ), 7, 21)

# make Q matrix
QQ <- matrix(list("s1","b12","b13","b14","b15","b16","b17",
                  "b12","s2","b23","b24","b25","b26","b27",
                  "b13","b23","s3","b34","b35","b36","b37",
                  "b14","b24","b34","s4","b45","b46","b47",
                  "b15","b25","b35","b45","s5","b56","b57",
                  "b16","b26","b36","b46","b56","s6","b67",
                  "b17","b27","b37","b47","b57","b67","s7"), 7, 7)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,0,
                  0,"v1",0,0,0,0,0,
                  0,0,"v1",0,0,0,0,
                  0,0,0,"v1",0,0,0,
                  0,0,0,0,"v2",0,0,
                  0,0,0,0,0,"v2",0,
                  0,0,0,0,0,0,"v2"), 7, 7)

# Model setup for MARSS

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
  R = RR, 
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
saveRDS(fit, file = "data_working/marss_fits/fit_031124_1state_cond_1ylegacy_binary_mBFGS.rds")

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
  R = RR, 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit) # AICc 1584.574   
MARSSaic(null.fit) # AICc 1653.978

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Fairly, VC sites still most wonky

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 2y legacy, 7 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_031124.rds")

# select sites
# include these sites only (7 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa_2ylegacy, 
    fire_pa_ppt_2ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa_2ylegacy,
                    fire_pa_ppt_2ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 479
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) # 0
sum(is.na(dat_cond_log[,cov_cols])) # 0
sum(is.infinite(dat_cond_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
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
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_pa_2ylegacy
  "fire_pa_2ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_2ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_2ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_2ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_2ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_2ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_2ylegacy_RED" ,
  # fire_pa_ppt_2ylegacy
  "fire_pa_ppt_2ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_ppt_2ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_ppt_2ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_ppt_2ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_ppt_2ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_ppt_2ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_ppt_2ylegacy_RED" ), 7, 21)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,0,
                  0,"v1",0,0,0,0,0,
                  0,0,"v1",0,0,0,0,
                  0,0,0,"v1",0,0,0,
                  0,0,0,0,"v2",0,0,
                  0,0,0,0,0,"v2",0,
                  0,0,0,0,0,0,"v2"), 7, 7)

# Model setup for MARSS

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", # 7 states
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = RR, # obs. error
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
        file = "data_working/marss_fits/fit_031124_7state_cond_2ylegacy_binary_mBFGS.rds")

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
  R = RR, 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit) # AICc 1620.29   
MARSSaic(null.fit) #AICc 1806.168

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)?
# Fairly.

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 2y legacy, 2 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_031124.rds")

# select sites
# include these sites only (7 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa_2ylegacy, 
    fire_pa_ppt_2ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa_2ylegacy,
                    fire_pa_ppt_2ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 479
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) # 0
sum(is.na(dat_cond_log[,cov_cols])) # 0
sum(is.infinite(dat_cond_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
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
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_pa_2ylegacy
  "fire_pa_2ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_2ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_2ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_2ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_2ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_2ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_2ylegacy_RED" ,
  # fire_pa_ppt_2ylegacy
  "fire_pa_ppt_2ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_ppt_2ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_ppt_2ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_ppt_2ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_ppt_2ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_ppt_2ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_ppt_2ylegacy_RED" ), 7, 21)

# make Q matrix
QQ <- matrix(list("s1","b12","b13","b14",0,    0,    0,
                  "b12","s2","b23","b24",0,    0,    0,
                  "b13","b23","s3","b34",0,    0,    0,
                  "b14","b24","b34","s4",0,    0,    0,
                  0,    0,    0,    0,"s5","b56","b57",
                  0,    0,    0,    0,"b56","s6","b67",
                  0,    0,    0,    0,"b57","b67","s7"), 7, 7)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,0,
                  0,"v1",0,0,0,0,0,
                  0,0,"v1",0,0,0,0,
                  0,0,0,"v1",0,0,0,
                  0,0,0,0,"v2",0,0,
                  0,0,0,0,0,"v2",0,
                  0,0,0,0,0,0,"v2"), 7, 7)

# Model setup for MARSS

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = QQ, # 2 states
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = RR, # obs. error
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
        file = "data_working/marss_fits/fit_031124_2state_cond_2ylegacy_binary_mBFGS.rds")

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
  R = RR, 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit) # AICc 1585.24  
MARSSaic(null.fit) #AICc 1666.072

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Not really

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Yes

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 2y legacy, 1 state ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_031124.rds")

# select sites
# include these sites only (7 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa_2ylegacy, 
    fire_pa_ppt_2ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa_2ylegacy,
                    fire_pa_ppt_2ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 479
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) # 0
sum(is.na(dat_cond_log[,cov_cols])) # 0
sum(is.infinite(dat_cond_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
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
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_pa_2ylegacy
  "fire_pa_2ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_2ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_2ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_2ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_2ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_2ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_2ylegacy_RED" ,
  # fire_pa_ppt_2ylegacy
  "fire_pa_ppt_2ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_ppt_2ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_ppt_2ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_ppt_2ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_ppt_2ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_ppt_2ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_ppt_2ylegacy_RED" ), 7, 21)

# make Q matrix
QQ <- matrix(list("s1","b12","b13","b14","b15","b16","b17",
                  "b12","s2","b23","b24","b25","b26","b27",
                  "b13","b23","s3","b34","b35","b36","b37",
                  "b14","b24","b34","s4","b45","b46","b47",
                  "b15","b25","b35","b45","s5","b56","b57",
                  "b16","b26","b36","b46","b56","s6","b67",
                  "b17","b27","b37","b47","b57","b67","s7"), 7, 7)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,0,
                  0,"v1",0,0,0,0,0,
                  0,0,"v1",0,0,0,0,
                  0,0,0,"v1",0,0,0,
                  0,0,0,0,"v2",0,0,
                  0,0,0,0,0,"v2",0,
                  0,0,0,0,0,0,"v2"), 7, 7)

# Model setup for MARSS

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
  R = RR, 
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
        file = "data_working/marss_fits/fit_031124_1state_cond_2ylegacy_binary_mBFGS.rds")

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
  R = RR, 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit) # AICc 1596.56   
MARSSaic(null.fit) #AICc 1653.978

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Fairly

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 3y legacy, 7 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_031124.rds")

# select sites
# include these sites only (7 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa_3ylegacy, 
    fire_pa_ppt_3ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa_3ylegacy,
                    fire_pa_ppt_3ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 479
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) # 0
sum(is.na(dat_cond_log[,cov_cols])) # 0
sum(is.infinite(dat_cond_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
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
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_pa_3ylegacy
  "fire_pa_3ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_3ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_3ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_3ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_3ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_3ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_3ylegacy_RED" ,
  # fire_pa_ppt_3ylegacy
  "fire_pa_ppt_3ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_ppt_3ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_ppt_3ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_ppt_3ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_ppt_3ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_ppt_3ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_ppt_3ylegacy_RED" ), 7, 21)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,0,
                  0,"v1",0,0,0,0,0,
                  0,0,"v1",0,0,0,0,
                  0,0,0,"v1",0,0,0,
                  0,0,0,0,"v2",0,0,
                  0,0,0,0,0,"v2",0,
                  0,0,0,0,0,0,"v2"), 7, 7)

# Model setup for MARSS

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", # 7 states
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = RR, # obs. error
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
        file = "data_working/marss_fits/fit_031124_7state_cond_3ylegacy_binary_mBFGS.rds")

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
  R = RR, 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit) # AICc 1612.722   
MARSSaic(null.fit) #AICc 1806.168

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Fairly

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 3y legacy, 2 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_031124.rds")

# select sites
# include these sites only (7 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa_3ylegacy, 
    fire_pa_ppt_3ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa_3ylegacy,
                    fire_pa_ppt_3ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 479
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) # 0
sum(is.na(dat_cond_log[,cov_cols])) # 0
sum(is.infinite(dat_cond_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
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
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_pa_3ylegacy
  "fire_pa_3ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_3ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_3ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_3ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_3ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_3ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_3ylegacy_RED" ,
  # fire_pa_ppt_3ylegacy
  "fire_pa_ppt_3ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_ppt_3ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_ppt_3ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_ppt_3ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_ppt_3ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_ppt_3ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_ppt_3ylegacy_RED" ), 7, 21)

# make Q matrix
QQ <- matrix(list("s1","b12","b13","b14",0,    0,    0,
                  "b12","s2","b23","b24",0,    0,    0,
                  "b13","b23","s3","b34",0,    0,    0,
                  "b14","b24","b34","s4",0,    0,    0,
                  0,    0,    0,    0,"s5","b56","b57",
                  0,    0,    0,    0,"b56","s6","b67",
                  0,    0,    0,    0,"b57","b67","s7"), 7, 7)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,0,
                  0,"v1",0,0,0,0,0,
                  0,0,"v1",0,0,0,0,
                  0,0,0,"v1",0,0,0,
                  0,0,0,0,"v2",0,0,
                  0,0,0,0,0,"v2",0,
                  0,0,0,0,0,0,"v2"), 7, 7)

# Model setup for MARSS

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = QQ, # 2 states
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = RR, 
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
        file = "data_working/marss_fits/fit_031124_2state_cond_3ylegacy_binary_mBFGS.rds")

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
  R = RR, 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit) # AICc 1578.385   
MARSSaic(null.fit) #AICc 1666.072

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plots 6 & 7(qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Yes

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 3y legacy, 1 state ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_031124.rds")

# select sites
# include these sites only (7 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa_3ylegacy, 
    fire_pa_ppt_3ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa_3ylegacy,
                    fire_pa_ppt_3ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 479
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) # 0
sum(is.na(dat_cond_log[,cov_cols])) # 0
sum(is.infinite(dat_cond_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
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
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_pa_3ylegacy
  "fire_pa_3ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_3ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_3ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_3ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_3ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_3ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_3ylegacy_RED" ,
  # fire_pa_ppt_3ylegacy
  "fire_pa_ppt_3ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_ppt_3ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_ppt_3ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_ppt_3ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_ppt_3ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_ppt_3ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_ppt_3ylegacy_RED" ), 7, 21)

# make Q matrix
QQ <- matrix(list("s1","b12","b13","b14","b15","b16","b17",
                  "b12","s2","b23","b24","b25","b26","b27",
                  "b13","b23","s3","b34","b35","b36","b37",
                  "b14","b24","b34","s4","b45","b46","b47",
                  "b15","b25","b35","b45","s5","b56","b57",
                  "b16","b26","b36","b46","b56","s6","b67",
                  "b17","b27","b37","b47","b57","b67","s7"), 7, 7)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,0,
                  0,"v1",0,0,0,0,0,
                  0,0,"v1",0,0,0,0,
                  0,0,0,"v1",0,0,0,
                  0,0,0,0,"v2",0,0,
                  0,0,0,0,0,"v2",0,
                  0,0,0,0,0,0,"v2"), 7, 7)

# Model setup for MARSS

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
  R = RR, 
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
        file = "data_working/marss_fits/fit_031124_1state_cond_3ylegacy_binary_mBFGS.rds")

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
  R = RR, 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit) # AICc 1587.914   
MARSSaic(null.fit) #AICc 1653.978

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Yes except RED

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 4y legacy, 7 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_031124.rds")

# select sites
# include these sites only (7 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa_4ylegacy, 
    fire_pa_ppt_4ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa_4ylegacy,
                    fire_pa_ppt_4ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 479
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) # 0
sum(is.na(dat_cond_log[,cov_cols])) # 0
sum(is.infinite(dat_cond_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
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
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_pa_4ylegacy
  "fire_pa_4ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_4ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_4ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_4ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_4ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_4ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_4ylegacy_RED" ,
  # fire_pa_ppt_4ylegacy
  "fire_pa_ppt_4ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_ppt_4ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_ppt_4ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_ppt_4ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_ppt_4ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_ppt_4ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_ppt_4ylegacy_RED" ), 7, 21)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,0,
                  0,"v1",0,0,0,0,0,
                  0,0,"v1",0,0,0,0,
                  0,0,0,"v1",0,0,0,
                  0,0,0,0,"v2",0,0,
                  0,0,0,0,0,"v2",0,
                  0,0,0,0,0,0,"v2"), 7, 7)

# Model setup for MARSS

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", # 7 states
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = RR, # obs. error
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
        file = "data_working/marss_fits/fit_031124_7state_cond_4ylegacy_binary_mBFGS.rds")

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
  R = RR, 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit) # AICc 1603.473   
MARSSaic(null.fit) #AICc 1806.168

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Ehhh they're ok

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 4y legacy, 2 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_031124.rds")

# select sites
# include these sites only (7 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa_4ylegacy, 
    fire_pa_ppt_4ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa_4ylegacy,
                    fire_pa_ppt_4ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 479
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) # 0
sum(is.na(dat_cond_log[,cov_cols])) # 0
sum(is.infinite(dat_cond_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
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
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_pa_4ylegacy
  "fire_pa_4ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_4ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_4ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_4ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_4ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_4ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_4ylegacy_RED" ,
  # fire_pa_ppt_4ylegacy
  "fire_pa_ppt_4ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_ppt_4ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_ppt_4ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_ppt_4ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_ppt_4ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_ppt_4ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_ppt_4ylegacy_RED" ), 7, 21)

# make Q matrix
QQ <- matrix(list("s1","b12","b13","b14",0,    0,    0,
                  "b12","s2","b23","b24",0,    0,    0,
                  "b13","b23","s3","b34",0,    0,    0,
                  "b14","b24","b34","s4",0,    0,    0,
                  0,    0,    0,    0,"s5","b56","b57",
                  0,    0,    0,    0,"b56","s6","b67",
                  0,    0,    0,    0,"b57","b67","s7"), 7, 7)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,0,
                  0,"v1",0,0,0,0,0,
                  0,0,"v1",0,0,0,0,
                  0,0,0,"v1",0,0,0,
                  0,0,0,0,"v2",0,0,
                  0,0,0,0,0,"v2",0,
                  0,0,0,0,0,0,"v2"), 7, 7)

# Model setup for MARSS

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = QQ, # 2 states
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = RR, # obs. error 
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
        file = "data_working/marss_fits/fit_031124_2state_cond_4ylegacy_binary_mBFGS.rds")

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
  R = RR, 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit) # AICc 1563.904 
MARSSaic(null.fit) # AICc 1666.072

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Yes, although RED a bit wonky

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 4y legacy, 1 state ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_031124.rds")

# select sites
# include these sites only (7 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa_4ylegacy, 
    fire_pa_ppt_4ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa_4ylegacy,
                    fire_pa_ppt_4ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 479
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols])) # 0
sum(is.na(dat_cond_log[,cov_cols])) # 0
sum(is.infinite(dat_cond_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
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
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_pa_4ylegacy
  "fire_pa_4ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_4ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_4ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_4ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_4ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_4ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_4ylegacy_RED" ,
  # fire_pa_ppt_4ylegacy
  "fire_pa_ppt_4ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_ppt_4ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_ppt_4ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_ppt_4ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_ppt_4ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_ppt_4ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_ppt_4ylegacy_RED" ), 7, 21)

# make Q matrix
QQ <- matrix(list("s1","b12","b13","b14","b15","b16","b17",
                  "b12","s2","b23","b24","b25","b26","b27",
                  "b13","b23","s3","b34","b35","b36","b37",
                  "b14","b24","b34","s4","b45","b46","b47",
                  "b15","b25","b35","b45","s5","b56","b57",
                  "b16","b26","b36","b46","b56","s6","b67",
                  "b17","b27","b37","b47","b57","b67","s7"), 7, 7)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,0,
                  0,"v1",0,0,0,0,0,
                  0,0,"v1",0,0,0,0,
                  0,0,0,"v1",0,0,0,
                  0,0,0,0,"v2",0,0,
                  0,0,0,0,0,"v2",0,
                  0,0,0,0,0,0,"v2"), 7, 7)

# Model setup for MARSS

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
  R = RR, # obs. error
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
        file = "data_working/marss_fits/fit_031124_1state_cond_4ylegacy_binary_mBFGS.rds")

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
  R = RR, 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit) # AICc 1569.659   
MARSSaic(null.fit) # AICc 1653.978

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No although RED is still a bit wonky

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Fairly

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No for the most part

#### 5y legacy, 7 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_031124.rds")

# select sites
# include these sites only (7 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa_5ylegacy, 
    fire_pa_ppt_5ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa_5ylegacy,
                    fire_pa_ppt_5ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 479
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
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_pa_5ylegacy
  "fire_pa_5ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_5ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_5ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_5ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_5ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_5ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_5ylegacy_RED" ,
  # fire_pa_ppt_5ylegacy
  "fire_pa_ppt_5ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_ppt_5ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_ppt_5ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_ppt_5ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_ppt_5ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_ppt_5ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_ppt_5ylegacy_RED" ), 7, 21)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,0,
                  0,"v1",0,0,0,0,0,
                  0,0,"v1",0,0,0,0,
                  0,0,0,"v1",0,0,0,
                  0,0,0,0,"v2",0,0,
                  0,0,0,0,0,"v2",0,
                  0,0,0,0,0,0,"v2"), 7, 7)

# Model setup for MARSS

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", # 7 states
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = RR, 
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
        file = "data_working/marss_fits/fit_031124_7state_cond_5ylegacy_binary_mBFGS.rds")

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
  R = RR, 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit) # AICc 1605.852   
MARSSaic(null.fit) #AICc 1806.168

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? For the most part

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 5y legacy, 2 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_031124.rds")

# select sites
# include these sites only (7 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa_5ylegacy, 
    fire_pa_ppt_5ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa_5ylegacy,
                    fire_pa_ppt_5ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 479
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
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_pa_5ylegacy
  "fire_pa_5ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_5ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_5ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_5ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_5ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_5ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_5ylegacy_RED" ,
  # fire_pa_ppt_5ylegacy
  "fire_pa_ppt_5ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_ppt_5ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_ppt_5ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_ppt_5ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_ppt_5ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_ppt_5ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_ppt_5ylegacy_RED" ), 7, 21)

# make Q matrix
QQ <- matrix(list("s1","b12","b13","b14",0,    0,   0, 
                  "b12","s2","b23","b24",0,    0,   0,
                  "b13","b23","s3","b34",0,    0,   0,
                  "b14","b24","b34","s4",0,    0,   0,
                  0,    0,    0,    0,"s5","b56","b57",
                  0,    0,    0,    0,"b56","s6","b67",
                  0,    0,    0,    0,"b57","b67","s7"), 7, 7)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,0,
                  0,"v1",0,0,0,0,0,
                  0,0,"v1",0,0,0,0,
                  0,0,0,"v1",0,0,0,
                  0,0,0,0,"v2",0,0,
                  0,0,0,0,0,"v2",0,
                  0,0,0,0,0,0,"v2"), 7, 7)

# Model setup for MARSS

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = QQ, # 2 state
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = RR, # obs. error
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
        file = "data_working/marss_fits/fit_031124_2state_cond_5ylegacy_binary_mBFGS.rds")

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
  R = RR, 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit) # AICc 1573.337   
MARSSaic(null.fit) # AICc 1666.072

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No except RED

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? For the most part

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 5y legacy, 1 state ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_031124.rds")

# select sites
# include these sites only (6 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa_5ylegacy, 
    fire_pa_ppt_5ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa_5ylegacy,
                    fire_pa_ppt_5ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:8)
cov_cols = c(9:29)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 479
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
  0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_pa_5ylegacy
  "fire_pa_5ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_5ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_5ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_5ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_5ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_5ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_5ylegacy_RED" ,
  # fire_pa_ppt_5ylegacy
  "fire_pa_ppt_5ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_pa_ppt_5ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_pa_ppt_5ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_pa_ppt_5ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_pa_ppt_5ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_pa_ppt_5ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_pa_ppt_5ylegacy_RED" ), 7, 21)

# make Q matrix
QQ <- matrix(list("s1","b12","b13","b14","b15","b16","b17",
                  "b12","s2","b23","b24","b25","b26","b27",
                  "b13","b23","s3","b34","b35","b36","b37",
                  "b14","b24","b34","s4","b45","b46","b47",
                  "b15","b25","b35","b45","s5","b56","b57",
                  "b16","b26","b36","b46","b56","s6","b67",
                  "b17","b27","b37","b47","b57","b67","s7"), 7, 7)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,0,
                  0,"v1",0,0,0,0,0,
                  0,0,"v1",0,0,0,0,
                  0,0,0,"v1",0,0,0,
                  0,0,0,0,"v2",0,0,
                  0,0,0,0,0,"v2",0,
                  0,0,0,0,0,0,"v2"), 7, 7)

# Model setup for MARSS

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
  R = RR, # obs. error
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
        file = "data_working/marss_fits/fit_031124_1state_cond_5ylegacy_binary_mBFGS.rds")

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
  R = RR, 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

MARSSaic(fit) # AICc 1570.87   
MARSSaic(null.fit) # AICc 1653.978 

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? For the most part

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### AICc Comparisons ####

# Compare all model fits for each legacy window to see which state 
# configuration was best.

# Presented here are the information criterion:
# AICc - Akaike Information Criterion adjusted for small sample sizes

# no legacy, 8 state
noleg_7state <- readRDS(file = "data_working/marss_fits/fit_031124_7state_cond_binary_mBFGS.rds")

# no legacy, 2 state
noleg_2state <- readRDS(file = "data_working/marss_fits/fit_031124_2state_cond_binary_mBFGS.rds")

# no legacy, 1 state
noleg_1state <- readRDS(file = "data_working/marss_fits/fit_031124_1state_cond_binary_mBFGS.rds")

MARSSaic(noleg_7state) # AICc 1616.778
MARSSaic(noleg_2state) # AICc 1592.289 
MARSSaic(noleg_1state) # AICc 1583.221 

# So, based on lowest AICc values, 1-state wins.

###

# 1y legacy, 7 state
leg1_7state <- readRDS(file = "data_working/marss_fits/fit_031124_7state_cond_1ylegacy_binary_mBFGS.rds")

# 1y legacy, 2 state
leg1_2state <- readRDS(file = "data_working/marss_fits/fit_031124_2state_cond_1ylegacy_binary_mBFGS.rds")

# 1y legacy, 1 state
leg1_1state <- readRDS(file = "data_working/marss_fits/fit_031124_1state_cond_1ylegacy_binary_mBFGS.rds")

MARSSaic(leg1_7state) # AICc 1618.708 
MARSSaic(leg1_2state) # AICc 1587.028  
MARSSaic(leg1_1state) # AICc 1584.574 

# So, based on lowest AICc values, 1-state wins.

###

# 2y legacy, 7 state
leg2_7state <- readRDS(file = "data_working/marss_fits/fit_031124_7state_cond_2ylegacy_binary_mBFGS.rds")

# 2y legacy, 2 state
leg2_2state <- readRDS(file = "data_working/marss_fits/fit_031124_2state_cond_2ylegacy_binary_mBFGS.rds")

# 2y legacy, 1 state
leg2_1state <- readRDS(file = "data_working/marss_fits/fit_031124_1state_cond_2ylegacy_binary_mBFGS.rds")

MARSSaic(leg2_7state) # AICc 1620.29  
MARSSaic(leg2_2state) # AICc 1585.24   
MARSSaic(leg2_1state) # AICc 1596.56   

# So, based on lowest AICc values, 2-state wins.

###

# 3y legacy, 8 state
leg3_7state <- readRDS(file = "data_working/marss_fits/fit_031124_7state_cond_3ylegacy_binary_mBFGS.rds")

# 3y legacy, 2 state
leg3_2state <- readRDS(file = "data_working/marss_fits/fit_031124_2state_cond_3ylegacy_binary_mBFGS.rds")

# 3y legacy, 1 state
leg3_1state <- readRDS(file = "data_working/marss_fits/fit_031124_1state_cond_3ylegacy_binary_mBFGS.rds")

MARSSaic(leg3_7state) # AICc 1612.722   
MARSSaic(leg3_2state) # AICc 1578.385   
MARSSaic(leg3_1state) # AICc 1587.914   

# So, based on lowest AICc values, 2-state wins.

###

# 4y legacy, 8 state
leg4_7state <- readRDS(file = "data_working/marss_fits/fit_031124_7state_cond_4ylegacy_binary_mBFGS.rds")

# 4y legacy, 2 state
leg4_2state <- readRDS(file = "data_working/marss_fits/fit_031124_2state_cond_4ylegacy_binary_mBFGS.rds")

# 4y legacy, 1 state
leg4_1state <- readRDS(file = "data_working/marss_fits/fit_031124_1state_cond_4ylegacy_binary_mBFGS.rds")

MARSSaic(leg4_7state) # AICc 1603.473  
MARSSaic(leg4_2state) # AICc 1563.904  
MARSSaic(leg4_1state) # AICc 1569.659   

# So, based on lowest AICc values, 2-state wins.

###

# 5y legacy, 7 state
leg5_7state <- readRDS(file = "data_working/marss_fits/fit_031124_7state_cond_5ylegacy_binary_mBFGS.rds")

# 5y legacy, 2 state
leg5_2state <- readRDS(file = "data_working/marss_fits/fit_031124_2state_cond_5ylegacy_binary_mBFGS.rds")

# 5y legacy, 1 state
leg5_1state <- readRDS(file = "data_working/marss_fits/fit_031124_1state_cond_5ylegacy_binary_mBFGS.rds")

MARSSaic(leg5_7state) # AICc 1605.852    
MARSSaic(leg5_2state) # AICc 1573.337   
MARSSaic(leg5_1state) # AICc 1570.87 

# So, based on lowest AICc values, 1-state wins.

###

# So, it would seem the 1 and 2 "state" model structures
# are equivalent, so I'll plot the 2 "state" model for
# consistency.

#### Results Figure ####

# Extract necessary confidence interval info
# First for 99th percentile credible intervals.
noleg_est99 <- MARSSparamCIs(noleg_2state, alpha = 0.01)
leg1y_est99 <- MARSSparamCIs(leg1_2state, alpha = 0.01)
leg2y_est99 <- MARSSparamCIs(leg2_2state, alpha = 0.01)
leg3y_est99 <- MARSSparamCIs(leg3_2state, alpha = 0.01)
leg4y_est99 <- MARSSparamCIs(leg4_2state, alpha = 0.01)
leg5y_est99 <- MARSSparamCIs(leg5_2state, alpha = 0.01)

# Format confidence intervals into dataframes
noleg_CI99 = data.frame(
  "Est." = noleg_est99$par$U,
  "Lower" = noleg_est99$par.lowCI$U,
  "Upper" = noleg_est99$par.upCI$U)
noleg_CI99$Parameter = rownames(noleg_CI99)
noleg_CI99[,1:3] = round(noleg_CI99[,1:3], 3)
noleg_CI99$Model = "immediate duration"

leg1y_CI99 = data.frame(
  "Est." = leg1y_est99$par$U,
  "Lower" = leg1y_est99$par.lowCI$U,
  "Upper" = leg1y_est99$par.upCI$U)
leg1y_CI99$Parameter = rownames(leg1y_CI99)
leg1y_CI99[,1:3] = round(leg1y_CI99[,1:3], 3)
leg1y_CI99$Model = "1 year duration"

leg2y_CI99 = data.frame(
  "Est." = leg2y_est99$par$U,
  "Lower" = leg2y_est99$par.lowCI$U,
  "Upper" = leg2y_est99$par.upCI$U)
leg2y_CI99$Parameter = rownames(leg2y_CI99)
leg2y_CI99[,1:3] = round(leg2y_CI99[,1:3], 3)
leg2y_CI99$Model = "2 year duration"

leg3y_CI99 = data.frame(
  "Est." = leg3y_est99$par$U,
  "Lower" = leg3y_est99$par.lowCI$U,
  "Upper" = leg3y_est99$par.upCI$U)
leg3y_CI99$Parameter = rownames(leg3y_CI99)
leg3y_CI99[,1:3] = round(leg3y_CI99[,1:3], 3)
leg3y_CI99$Model = "3 year duration"

leg4y_CI99 = data.frame(
  "Est." = leg4y_est99$par$U,
  "Lower" = leg4y_est99$par.lowCI$U,
  "Upper" = leg4y_est99$par.upCI$U)
leg4y_CI99$Parameter = rownames(leg4y_CI99)
leg4y_CI99[,1:3] = round(leg4y_CI99[,1:3], 3)
leg4y_CI99$Model = "4 year duration"

leg5y_CI99 = data.frame(
  "Est." = leg5y_est99$par$U,
  "Lower" = leg5y_est99$par.lowCI$U,
  "Upper" = leg5y_est99$par.upCI$U)
leg5y_CI99$Parameter = rownames(leg5y_CI99)
leg5y_CI99[,1:3] = round(leg5y_CI99[,1:3], 3)
leg5y_CI99$Model = "5 year duration"

# Then for 95th percentile credible intervals.
noleg_est95 <- MARSSparamCIs(noleg_2state, alpha = 0.05)
leg1y_est95 <- MARSSparamCIs(leg1_2state, alpha = 0.05)
leg2y_est95 <- MARSSparamCIs(leg2_2state, alpha = 0.05)
leg3y_est95 <- MARSSparamCIs(leg3_2state, alpha = 0.05)
leg4y_est95 <- MARSSparamCIs(leg4_2state, alpha = 0.05)
leg5y_est95 <- MARSSparamCIs(leg5_2state, alpha = 0.05)

# Format confidence intervals into dataframes
noleg_CI95 = data.frame(
  "Est." = noleg_est95$par$U,
  "Lower" = noleg_est95$par.lowCI$U,
  "Upper" = noleg_est95$par.upCI$U)
noleg_CI95$Parameter = rownames(noleg_CI95)
noleg_CI95[,1:3] = round(noleg_CI95[,1:3], 3)
noleg_CI95$Model = "immediate duration"

leg1y_CI95 = data.frame(
  "Est." = leg1y_est95$par$U,
  "Lower" = leg1y_est95$par.lowCI$U,
  "Upper" = leg1y_est95$par.upCI$U)
leg1y_CI95$Parameter = rownames(leg1y_CI95)
leg1y_CI95[,1:3] = round(leg1y_CI95[,1:3], 3)
leg1y_CI95$Model = "1 year duration"

leg2y_CI95 = data.frame(
  "Est." = leg2y_est95$par$U,
  "Lower" = leg2y_est95$par.lowCI$U,
  "Upper" = leg2y_est95$par.upCI$U)
leg2y_CI95$Parameter = rownames(leg2y_CI95)
leg2y_CI95[,1:3] = round(leg2y_CI95[,1:3], 3)
leg2y_CI95$Model = "2 year duration"

leg3y_CI95 = data.frame(
  "Est." = leg3y_est95$par$U,
  "Lower" = leg3y_est95$par.lowCI$U,
  "Upper" = leg3y_est95$par.upCI$U)
leg3y_CI95$Parameter = rownames(leg3y_CI95)
leg3y_CI95[,1:3] = round(leg3y_CI95[,1:3], 3)
leg3y_CI95$Model = "3 year duration"

leg4y_CI95 = data.frame(
  "Est." = leg4y_est95$par$U,
  "Lower" = leg4y_est95$par.lowCI$U,
  "Upper" = leg4y_est95$par.upCI$U)
leg4y_CI95$Parameter = rownames(leg4y_CI95)
leg4y_CI95[,1:3] = round(leg4y_CI95[,1:3], 3)
leg4y_CI95$Model = "4 year duration"

leg5y_CI95 = data.frame(
  "Est." = leg5y_est95$par$U,
  "Lower" = leg5y_est95$par.lowCI$U,
  "Upper" = leg5y_est95$par.upCI$U)
leg5y_CI95$Parameter = rownames(leg5y_CI95)
leg5y_CI95[,1:3] = round(leg5y_CI95[,1:3], 3)
leg5y_CI95$Model = "5 year duration"

# Bind all together
CIs99 = rbind(noleg_CI99, leg1y_CI99, leg2y_CI99, 
              leg3y_CI99, leg4y_CI99,leg5y_CI99)

CIs95 = rbind(noleg_CI95, leg1y_CI95, leg2y_CI95, 
              leg3y_CI95, leg4y_CI95,leg5y_CI95)

# Add column for site names
CIs99$Stream = gsub("_","",str_sub(CIs99$Parameter, start= -4))
CIs95$Stream = gsub("_","",str_sub(CIs95$Parameter, start= -4))

# Join both sets together
CIs99 <- CIs99 %>%
  rename(Lower99 = "Lower",
         Upper99 = "Upper")

CIs95 <- CIs95 %>%
  rename(Lower95 = "Lower",
         Upper95 = "Upper")

CIs <- full_join(CIs99, CIs95)

# Simplify parameter names
CIs$Parm_simple = c(rep("Ppt",7),
                    rep("Fire",7),
                    rep("Ppt x Fire",7),
                    
                    rep("Ppt",7),
                    rep("Fire",7),
                    rep("Ppt x Fire",7),
                    
                    rep("Ppt",7),
                    rep("Fire",7),
                    rep("Ppt x Fire",7),
                    
                    rep("Ppt",7),
                    rep("Fire",7),
                    rep("Ppt x Fire",7),
                    
                    rep("Ppt",7),
                    rep("Fire",7),
                    rep("Ppt x Fire",7),
                    
                    rep("Ppt",7),
                    rep("Fire",7),
                    rep("Ppt x Fire",7))

# Add column to designate those sites at which effects are significant.
CIs <- CIs %>%
  # make sure to go from most specific to least here, because the first
  # assigned grouping will stick, making categorical assignment here challenging
  mutate(sig = factor(case_when(`Est.` > 0 & 
                                  Lower95 > 0 & Lower99 > 0 & 
                                  Upper95 > 0 & Upper99 > 0 ~ "sig_pos",
                                `Est.` < 0 & 
                                  Lower95 < 0 & Lower99 < 0 &
                                  Upper95 < 0 & Upper99 < 0 ~ "sig_neg",
                                `Est.` > 0 & 
                                  Lower95 > 0 & 
                                  Upper95 > 0  ~ "weak_sig_pos",
                                `Est.` < 0 & 
                                  Lower95 < 0 &
                                  Upper95 < 0 ~ "weak_sig_neg",
                                TRUE ~ "not_sig"), 
                      levels = c("weak_sig_pos", "sig_pos", 
                                 "not_sig", 
                                 "sig_neg", "weak_sig_neg"))) %>%
  mutate(region = case_when(Stream %in% c("AB00", "GV01", 
                                          "HO00", "RS02") ~ "Mediterranean",
                            TRUE ~ "Monsoonal")) %>%
  mutate(site = factor(Stream,
                       levels = c("AB00", "GV01", "HO00", "RS02",
                                  "EFJ", "RED", "RSAW"))) %>%
  mutate(model = factor(Model, levels = c("immediate duration",
                                          "1 year duration",
                                          "2 year duration",
                                          "3 year duration",
                                          "4 year duration",
                                          "5 year duration"))) %>%
  mutate(Parm_simple_f = factor(case_when(Parm_simple == "Ppt x Fire" ~ "Ppt x Fire",
                                          Parm_simple == "Fire" ~ "Fire",
                                          Parm_simple == "Ppt" ~ "Ppt"),
                                levels = c("Ppt x Fire",
                                           "Fire",
                                           "Ppt")))

# Plot results
(SpCond_fig <- ggplot(CIs, aes(x = Parm_simple_f,
                            y = Est., fill = sig, 
                            shape = site, group = desc(site))) + 
    # coloring by both percentiles but using 99th perc. error bars to be most conservative
    geom_errorbar(aes(ymin = Lower99, ymax = Upper99),
                  position=position_dodge(width = 0.5), width = 0) +
    geom_point(position = position_dodge(width = 0.5), 
               alpha = 0.8, size = 8) + 
    scale_shape_manual(values = c(21, 22, 23, 24,
                                  21, 22, 23)) +
    scale_fill_manual(values = c("gray50", "black", "white", "black", "gray50")) +
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
    facet_grid(region~model))

# Export plot.
# ggsave(("MARSS_SpCon_binary_031124.png"),
#        path = "figures",
#        width = 65,
#        height = 24,
#        units = "cm"
# )

# End of script.
