# Conductivity MARSS models
# with fire x ppt interactions and legacy effects
# as well as multiple state structure for CA/NM sites
# Script started August 14, 2023 by Heili Lowman

# This script will run 18 Conductivity MARSS models and
# calculate summary statistics for conductivity data.
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

#### Summary Stats ####

# load data with fire x ppt interactions and legacy effects
dat <- readRDS("data_working/marss_data_sb_vc_091123.rds")

# select sites
# include these sites only (7 total - these have the longest most complete ts 
# for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# EFJ, RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "EFJ", "RED", "RSAW")
dat = dat[dat$site %in% sitez,]

# Create summary table.
dat_summary <- dat %>%
  group_by(region) %>%
  summarize(min_cond = min(mean_cond_uScm, na.rm = TRUE),
            max_cond = max(mean_cond_uScm, na.rm = TRUE),
            mean_cond = mean(mean_cond_uScm, na.rm = TRUE),
            sd_cond = sd(mean_cond_uScm, na.rm = TRUE),
            med_cond = median(mean_cond_uScm, na.rm = TRUE)) %>%
  ungroup()

#### 0y legacy, 7 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat <- readRDS("data_working/marss_data_sb_vc_091123.rds")

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
    fire_perc_ws, fire_perc_ws_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws, fire_perc_ws_ppt)) 

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
  # fire_perc_ws: fire effect in 2 m window
  "fire_perc_ws_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_RED" ,
  # fire_perc_ws_ppt: interaction of cum. ppt with fire % burn in 6 m window
  "fire_perc_ws_ppt_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_RED" ), 7, 21)

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
saveRDS(fit, file = "data_working/marss_fits/fit_092523_7state_cond_mBFGS.rds")

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

MARSSaic(fit) # AICc 1623.89 
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
# VC sites are the wonkiest

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 0y legacy, 2 states ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat <- readRDS("data_working/marss_data_sb_vc_091123.rds")

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
    fire_perc_ws, fire_perc_ws_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws, fire_perc_ws_ppt)) 

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
  # fire_perc_ws: fire effect in 2 m window
  "fire_perc_ws_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_RED" ,
  # fire_perc_ws_ppt: interaction of cum. ppt with fire % burn in 6 m window
  "fire_perc_ws_ppt_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_RED" ), 7, 21)

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
saveRDS(fit, file = "data_working/marss_fits/fit_092523_2state_cond_mBFGS.rds")

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

MARSSaic(fit) # AICc 1600.19  
MARSSaic(null.fit) # AICc 1666.072   

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight line)? Yes for the most part - RED looking a bit strange

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 0y legacy, 1 state ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat <- readRDS("data_working/marss_data_sb_vc_091123.rds")

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
    fire_perc_ws, fire_perc_ws_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws, fire_perc_ws_ppt)) 

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
  0,0,0,0,0,0,"cumulative_precip_mm_RED",
  # fire_perc_ws: fire effect in 2 m window
  "fire_perc_ws_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_RED",
  # fire_perc_ws_ppt: interaction of cum. ppt with fire % burn in 6 m window
  "fire_perc_ws_ppt_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_RED"), 7, 21)

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
saveRDS(fit, file = "data_working/marss_fits/fit_092523_1state_cond_mBFGS.rds")

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

MARSSaic(fit) # AICc 1597.105  
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

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No except some wonky VC sites

#### 1y legacy, 7 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_091123.rds")

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
    fire_perc_ws_1ylegacy, 
    fire_perc_ws_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_1ylegacy,
                    fire_perc_ws_ppt_1ylegacy)) 

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
  # fire_perc_ws_1ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_1ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 7, 21)

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
saveRDS(fit, file = "data_working/marss_fits/fit_092523_7state_cond_1ylegacy_mBFGS.rds")

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

MARSSaic(fit) # AICc 1616.927
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
dat = readRDS("data_working/marss_data_sb_vc_091123.rds")

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
    fire_perc_ws_1ylegacy, 
    fire_perc_ws_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_1ylegacy,
                    fire_perc_ws_ppt_1ylegacy)) 

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
  # fire_perc_ws_1ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_1ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 7, 21)

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
saveRDS(fit, file = "data_working/marss_fits/fit_092523_2state_cond_1ylegacy_mBFGS.rds")

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

MARSSaic(fit) # AICc 1587.286  
MARSSaic(null.fit) # AICc 1666.072 

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No
# RED still looking a bit strange - but this is likely due to
# missingness

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Fairly

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 1y legacy, 1 state ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_091123.rds")

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
    fire_perc_ws_1ylegacy, 
    fire_perc_ws_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_1ylegacy,
                    fire_perc_ws_ppt_1ylegacy)) 

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
  # fire_perc_ws_1ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_1ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 7, 21)

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
saveRDS(fit, file = "data_working/marss_fits/fit_092523_1state_cond_1ylegacy_mBFGS.rds")

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

MARSSaic(fit) # AICc 1591.56 
MARSSaic(null.fit) # AICc 1653.978

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns? No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No
# RED still looking a bit strange - but this is likely due to
# missingness

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Fairly, VC sites still a bit wonky

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No but a bit wonky

#### 2y legacy, 7 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_091123.rds")

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
    fire_perc_ws_2ylegacy, 
    fire_perc_ws_ppt_2ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_2ylegacy,
                    fire_perc_ws_ppt_2ylegacy)) 

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
  # fire_perc_ws_2ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_2ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 7, 21)

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
        file = "data_working/marss_fits/fit_092523_7state_cond_2ylegacy_mBFGS.rds")

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

MARSSaic(fit) # AICc 1613.598
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
# More so than in some diagnoses above.

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 2y legacy, 2 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_091123.rds")

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
    fire_perc_ws_2ylegacy, 
    fire_perc_ws_ppt_2ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_2ylegacy,
                    fire_perc_ws_ppt_2ylegacy)) 

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
  # fire_perc_ws_2ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_2ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 7, 21)

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
        file = "data_working/marss_fits/fit_092523_2state_cond_2ylegacy_mBFGS.rds")

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

MARSSaic(fit) # AICc 1578.748
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
# except a few wonky sites at RS02 and RSAW

#### 2y legacy, 1 state ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_091123.rds")

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
    fire_perc_ws_2ylegacy, 
    fire_perc_ws_ppt_2ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_2ylegacy,
                    fire_perc_ws_ppt_2ylegacy)) 

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
  # fire_perc_ws_2ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_2ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 7, 21)

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
        file = "data_working/marss_fits/fit_092523_1state_cond_2ylegacy_mBFGS.rds")

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

MARSSaic(fit) # AICc 1592.993
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
# RS02 still a bit wonky

#### 3y legacy, 7 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_091123.rds")

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
    fire_perc_ws_3ylegacy, 
    fire_perc_ws_ppt_3ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_3ylegacy,
                    fire_perc_ws_ppt_3ylegacy)) 

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
  # fire_perc_ws_3ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_3ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 7, 21)

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
        file = "data_working/marss_fits/fit_092523_7state_cond_3ylegacy_mBFGS.rds")

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

MARSSaic(fit) # AICc 1610.204
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
dat = readRDS("data_working/marss_data_sb_vc_091123.rds")

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
    fire_perc_ws_3ylegacy, 
    fire_perc_ws_ppt_3ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_3ylegacy,
                    fire_perc_ws_ppt_3ylegacy)) 

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
  # fire_perc_ws_3ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_3ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 7, 21)

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
        file = "data_working/marss_fits/fit_092523_2state_cond_3ylegacy_mBFGS.rds")

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

MARSSaic(fit) # AICc 1574.025
MARSSaic(null.fit) #AICc 1666.072

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? Only at RED

# Plots 6 & 7(qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Yes, except RED

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 3y legacy, 1 state ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_091123.rds")

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
    fire_perc_ws_3ylegacy, 
    fire_perc_ws_ppt_3ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_3ylegacy,
                    fire_perc_ws_ppt_3ylegacy)) 

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
  # fire_perc_ws_3ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_3ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 7, 21)

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
        file = "data_working/marss_fits/fit_092523_1state_cond_3ylegacy_mBFGS.rds")

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

MARSSaic(fit) # AICc 1584.687
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
dat = readRDS("data_working/marss_data_sb_vc_091123.rds")

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
    fire_perc_ws_4ylegacy, 
    fire_perc_ws_ppt_4ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_4ylegacy,
                    fire_perc_ws_ppt_4ylegacy)) 

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
  # fire_perc_ws_4ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_4ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 7, 21)

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
        file = "data_working/marss_fits/fit_092523_7state_cond_4ylegacy_mBFGS.rds")

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

MARSSaic(fit) # AICc 1605.075
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
dat = readRDS("data_working/marss_data_sb_vc_091123.rds")

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
    fire_perc_ws_4ylegacy, 
    fire_perc_ws_ppt_4ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_4ylegacy,
                    fire_perc_ws_ppt_4ylegacy)) 

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
  # fire_perc_ws_4ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_4ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 7, 21)

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
        file = "data_working/marss_fits/fit_092523_2state_cond_4ylegacy_mBFGS.rds")

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

MARSSaic(fit) # AICc 1566.696
MARSSaic(null.fit) # AICc 1666.072

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? RED is behaving strangely again

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Yes except RED

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### 4y legacy, 1 state ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_091123.rds")

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
    fire_perc_ws_4ylegacy, 
    fire_perc_ws_ppt_4ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_4ylegacy,
                    fire_perc_ws_ppt_4ylegacy)) 

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
  # fire_perc_ws_4ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_4ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 7, 21)

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
        file = "data_working/marss_fits/fit_092523_1state_cond_4ylegacy_mBFGS.rds")

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

MARSSaic(fit) # AICc 1572.164
MARSSaic(null.fit) # AICc 1653.978

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? RED is still behaving strangely

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? Fairly

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No except RSAW eee

#### 5y legacy, 6 states ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_091123.rds")

# select sites
# include these sites only (6 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "RED", "RSAW")
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
resp_cols = c(2:7)
cov_cols = c(8:25)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 392
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
  "cumulative_precip_mm_AB00",0,0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,0,
  0,0,0,"cumulative_precip_mm_RS02",0,0,
  0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_perc_ws_5ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RED",
  # fire_perc_ws_ppt_5ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 6, 18)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,
                  0,"v1",0,0,0,0,
                  0,0,"v1",0,0,0,
                  0,0,0,"v1",0,0,
                  0,0,0,0,"v2",0,
                  0,0,0,0,0,"v2"), 6, 6)

# Model setup for MARSS

mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", # 8 states
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
        file = "data_working/marss_fits/fit_092523_6state_cond_5ylegacy_mBFGS.rds")

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

MARSSaic(fit) # AICc 1506.619
MARSSaic(null.fit) #AICc 1698.884

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
dat = readRDS("data_working/marss_data_sb_vc_091123.rds")

# select sites
# include these sites only (6 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "RED", "RSAW")
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
resp_cols = c(2:7)
cov_cols = c(8:25)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 392
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
  "cumulative_precip_mm_AB00",0,0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,0,
  0,0,0,"cumulative_precip_mm_RS02",0,0,
  0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,"cumulative_precip_mm_RED",
  # fire_perc_ws_5ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RED",
  # fire_perc_ws_ppt_5ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 6, 18)

# make Q matrix
QQ <- matrix(list("s1","b12","b13","b14",0,    0,   
                  "b12","s2","b23","b24",0,    0,   
                  "b13","b23","s3","b34",0,    0,  
                  "b14","b24","b34","s4",0,    0,   
                  0,    0,    0,    0,"s5","b56",
                  0,    0,    0,    0,"b56","s6"), 6, 6)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,
                  0,"v1",0,0,0,0,
                  0,0,"v1",0,0,0,
                  0,0,0,"v1",0,0,
                  0,0,0,0,"v2",0,
                  0,0,0,0,0,"v2"), 6, 6)

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
        file = "data_working/marss_fits/fit_092523_2state_cond_5ylegacy_mBFGS.rds")

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

MARSSaic(fit) # AICc 1483.119
MARSSaic(null.fit) # AICc 1571.249

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
dat = readRDS("data_working/marss_data_sb_vc_091123.rds")

# select sites
# include these sites only (6 total - these have the longest most
# complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, & RS02 = SB
# RED, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "RS02",
          "RED", "RSAW")
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
resp_cols = c(2:7)
cov_cols = c(8:25)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols])) # 0
sum(is.na(dat_cond_log[,resp_cols])) # 392
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
  "cumulative_precip_mm_AB00",0,0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,0,
  0,0,0,"cumulative_precip_mm_RS02",0,0,
  0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_perc_ws_5ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_ppt_5ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 6, 18)

# make Q matrix
QQ <- matrix(list("s1","b12","b13","b14","b15","b16",
                  "b12","s2","b23","b24","b25","b26",
                  "b13","b23","s3","b34","b35","b36",
                  "b14","b24","b34","s4","b45","b46",
                  "b15","b25","b35","b45","s5","b56",
                  "b16","b26","b36","b46","b56","s6"), 6, 6)

# Make R matrix
RR <- matrix(list("v1",0,0,0,0,0,
                  0,"v1",0,0,0,0,
                  0,0,"v1",0,0,0,
                  0,0,0,"v1",0,0,
                  0,0,0,0,"v2",0,
                  0,0,0,0,0,"v2"), 6, 6)

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
        file = "data_working/marss_fits/fit_092523_1state_cond_5ylegacy_mBFGS.rds")

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

MARSSaic(fit) # AICc 1487.4
MARSSaic(null.fit) # AICc 1582.582

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 4 (std.model.resids.ytT): Do resids have temporal patterns?  No
# Do 95% of resids fall withing the CIs? Yes

# Plot 5 (std.state.resids.xtT): Any outliers? No

# Plots 6 & 7 (qqplot.std.model.resids.ytt1: Are resids normal (straight lines)? For the most part
# but RED not looking great

# Plot 8 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation? No

#### AICc Comparisons ####

# Compare all model fits for each legacy window to see which state 
# configuration was best.

# Presented here are the information criterion:
# AICc - Akaike Information Criterion adjusted for small sample sizes

# no legacy, 8 state
noleg_7state <- readRDS(file = "data_working/marss_fits/fit_092523_7state_cond_mBFGS.rds")

# no legacy, 2 state
noleg_2state <- readRDS(file = "data_working/marss_fits/fit_092523_2state_cond_mBFGS.rds")

# no legacy, 1 state
noleg_1state <- readRDS(file = "data_working/marss_fits/fit_092523_1state_cond_mBFGS.rds")

MARSSaic(noleg_7state) # AICc 1623.89
MARSSaic(noleg_2state) # AICc 1600.19
MARSSaic(noleg_1state) # AICc 1597.105

# So, based on lowest AICc values, 1-state wins.

###

# 1y legacy, 7 state
leg1_7state <- readRDS(file = "data_working/marss_fits/fit_092523_7state_cond_1ylegacy_mBFGS.rds")

# 1y legacy, 2 state
leg1_2state <- readRDS(file = "data_working/marss_fits/fit_092523_2state_cond_1ylegacy_mBFGS.rds")

# 1y legacy, 1 state
leg1_1state <- readRDS(file = "data_working/marss_fits/fit_092523_1state_cond_1ylegacy_mBFGS.rds")

MARSSaic(leg1_7state) # AICc 1616.927
MARSSaic(leg1_2state) # AICc 1587.286
MARSSaic(leg1_1state) # AICc 1591.56

# So, based on lowest AICc values, 2-state wins.

###

# 2y legacy, 7 state
leg2_7state <- readRDS(file = "data_working/marss_fits/fit_092523_7state_cond_2ylegacy_mBFGS.rds")

# 2y legacy, 2 state
leg2_2state <- readRDS(file = "data_working/marss_fits/fit_092523_2state_cond_2ylegacy_mBFGS.rds")

# 2y legacy, 1 state
leg2_1state <- readRDS(file = "data_working/marss_fits/fit_092523_1state_cond_2ylegacy_mBFGS.rds")

MARSSaic(leg2_7state) # AICc 1613.598
MARSSaic(leg2_2state) # AICc 1578.748
MARSSaic(leg2_1state) # AICc 1592.993

# So, based on lowest AICc values, 2-state wins.

###

# 3y legacy, 8 state
leg3_7state <- readRDS(file = "data_working/marss_fits/fit_092523_7state_cond_3ylegacy_mBFGS.rds")

# 3y legacy, 2 state
leg3_2state <- readRDS(file = "data_working/marss_fits/fit_092523_2state_cond_3ylegacy_mBFGS.rds")

# 3y legacy, 1 state
leg3_1state <- readRDS(file = "data_working/marss_fits/fit_092523_1state_cond_3ylegacy_mBFGS.rds")

MARSSaic(leg3_7state) # AICc 1610.204
MARSSaic(leg3_2state) # AICc 1574.025
MARSSaic(leg3_1state) # AICc 1584.687

# So, based on lowest AICc values, 2-state wins.

###

# 4y legacy, 8 state
leg4_7state <- readRDS(file = "data_working/marss_fits/fit_092523_7state_cond_4ylegacy_mBFGS.rds")

# 4y legacy, 2 state
leg4_2state <- readRDS(file = "data_working/marss_fits/fit_092523_2state_cond_4ylegacy_mBFGS.rds")

# 4y legacy, 1 state
leg4_1state <- readRDS(file = "data_working/marss_fits/fit_092523_1state_cond_4ylegacy_mBFGS.rds")

MARSSaic(leg4_7state) # AICc 1605.075
MARSSaic(leg4_2state) # AICc 1566.696
MARSSaic(leg4_1state) # AICc 1572.164

# So, based on lowest AICc values, 2-state wins.

###

# 5y legacy, 6 state
leg5_6state <- readRDS(file = "data_working/marss_fits/fit_092523_6state_cond_5ylegacy_mBFGS.rds")

# 5y legacy, 2 state
leg5_2state <- readRDS(file = "data_working/marss_fits/fit_092523_2state_cond_5ylegacy_mBFGS.rds")

# 5y legacy, 1 state
leg5_1state <- readRDS(file = "data_working/marss_fits/fit_092523_1state_cond_5ylegacy_mBFGS.rds")

MARSSaic(leg5_6state) # AICc 1506.619
MARSSaic(leg5_2state) # AICc 1483.119
MARSSaic(leg5_1state) # AICc 1487.4

# So, based on lowest AICc values, 2-state wins.

###

# So, it would seem the 2 "state" model structure wins out.

#### Results Figure ####

# For presentation consistency, I will only be creating figures with a
# single state configuration, whichever yielded the most parsimonious
# models. So, in this case, all conductivity figures will 
# represent the "2 state" scenario.

# Extract necessary confidence interval info
noleg_est <- MARSSparamCIs(noleg_2state, alpha = 0.01)
leg1y_est <- MARSSparamCIs(leg1_2state, alpha = 0.01)
leg2y_est <- MARSSparamCIs(leg2_2state, alpha = 0.01)
leg3y_est <- MARSSparamCIs(leg3_2state, alpha = 0.01)
leg4y_est <- MARSSparamCIs(leg4_2state, alpha = 0.01)
leg5y_est <- MARSSparamCIs(leg5_2state, alpha = 0.01)

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
CIs$Parm_simple = c(rep("Ppt",7),
                    rep("Perc. burn",7),
                    rep("Ppt x Perc. burn",7),
                    
                    rep("Ppt",7),
                    rep("Perc. burn",7),
                    rep("Ppt x Perc. burn",7),
                    
                    rep("Ppt",7),
                    rep("Perc. burn",7),
                    rep("Ppt x Perc. burn",7),
                    
                    rep("Ppt",7),
                    rep("Perc. burn",7),
                    rep("Ppt x Perc. burn",7),
                    
                    rep("Ppt",7),
                    rep("Perc. burn",7),
                    rep("Ppt x Perc. burn",7),
                    
                    rep("Ppt",6),
                    rep("Perc. burn",6),
                    rep("Ppt x Perc. burn",6))

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
                                  "EFJ", "RED", "RSA", "RSAW"))) %>%
  mutate(model = factor(Model, levels = c("immediate window",
                                          "1 year window",
                                          "2 year window",
                                          "3 year window",
                                          "4 year window",
                                          "5 year window")))

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
    facet_grid(region~model))

# Export plot.
# ggsave(("MARSS_SpCon_101023.png"),
#        path = "figures",
#        width = 65,
#        height = 24,
#        units = "cm"
# )

# End of script.
