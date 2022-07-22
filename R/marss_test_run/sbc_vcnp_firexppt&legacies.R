# SB & VC MARSS models with fire x ppt interactions and legacy effects
# Script started July 19, 2022
# Heili Lowman, Alex Webster

# The following script will run a MARSS analysis at SBC LTER & NM Valles Caldera sites for the CRASS project.

#### Data Notes ####

# A few notes about the process below:
# Data was tidied in sbc_marss_model.R script and exported to file:
# marss_data_sb_vc_072222.rds (as of July 22, 2022 when errors in NM ppt data were identified and fixed)
# It is imported from that file in this script.

# The first goal of this script is to demo a MARSS model with 12 states (each watershed as a unique state) and the following predictive vars:
# 1 - cum. monthly ppt
# 2 - fire_pa_6m  (1=fire ignited in prior 6 m, 0=no fire ignitions in 6 m, one col for all fires)
# 3 - fire_pa_6m_1ylegacy (1=fire ignited 6 m - 1 y ago, 0=no fire ignitions 6 m - 1 y ago)
# 4 - fire_pa_6m_2ylegacy (1=fire ignited 1.5 - 2 y ago, 0=no fire ignitions 1.5 - 2 y ago)
# 5 - fire_pa_6m_3ylegacy (1=fire ignited 2.5 - 3 y ago, 0=no fire ignitions 2.5 - 3 y ago)
# 6 - fire_pa_6m_4ylegacy (1=fire ignited 3.5 - 4 y ago, 0=no fire ignitions 3.5 - 4 y ago)
# 7 - fire_pa_6m_5ylegacy (1=fire ignited 4.5 - 5 y ago, 0=no fire ignitions 4.5 - 5 y ago)
# 8 - fire_pa_6m_ppt (fire ignited in prior 6 m * ppt this month)
# 9 - fire_pa_6m_ppt_1ylegacy (fire ignited 6 m - 1 y ago * ppt this month)
# 10 - fire_pa_6m_ppt_2ylegacy (fire ignited 1.5 - 2 y ago * ppt this month)
# 11 - fire_pa_6m_ppt_3ylegacy (fire ignited 2.5 - 3 y ago * ppt this month)
# 12 - fire_pa_6m_ppt_4ylegacy (fire ignited 3.5 - 4 y ago * ppt this month)
# 13 - fire_pa_6m_ppt_5ylegacy (fire ignited 4.5 - 5 y ago * ppt this month)

# Other data notes (copied from the sbc_marss_model.R script):
# Following a discussion with John Melack, the following sites have both enough chemistry
# and precip data, as well as fires that occurred within that timeframe to be included in the
# MARSS analysis performed below.

# Arroyo Burro - AB00
# Atascadero - AT07
# Gaviota - GV01
# Arroyo Hondo - HO00
# Mission Creek (at Rocky Nook) - MC06
# Refugio - RG01
# Rattlesnake - RS02
# San Pedro - SP02

# Here is a list of the wildfire events in the SBC:
# Fire Name (Start Date) - affected watersheds
# Gaviota (6/4/2004) - Canada De Santa Anita, Canada De La Gaviota, Tajiguas Creek
# Sherpa (6/15/2016) - Tajiguas Creek, Dos Pueblos Canyon
# Whittier (7/7/2017) - Tajiguas Creek, Dos Pueblos Canyon
# Gap (7/1/2008) - Dos Pueblos Canyon, San Pedro Creek
# Cave (11/25/2019) - Atascadero Creek
# Jesusita (5/5/2009) - Atascadero Creek, Mission Creek
# Tea (11/13/2008) - Mission Creek

# The following site may be added later, but a longer precipitation record needs to be
# identified for it - Bell Canyon (BC02).

# See "data_raw/VCNP_sonde_site_codes_names.csv" for other possible site names used across other files (E.g., sonde data, GIS, etc.)
# "Redondo Creek" ~ "RED",
# "East Fork Jemez River" ~ "EFJ",
# "San Antonio - West" ~ "RSAW",
# "San Antonio Creek - Toledo" ~ "RSA",
# "Indios Creek" ~ "IND",
# "Indios Creek - Post Fire (Below Burn)" ~ "IND_BB",
# "Indios Creek - above burn" ~ "IND_AB",
# "Sulfur Creek" ~ "SULF"


#### Load packages ####
library(tidyverse)
library(lubridate)
library(MARSS)
library(naniar) 

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))


#### Import data ####

# dat1 = readRDS("data_working/marss_data_sb_vc_060622.rds")
dat1 = readRDS("data_working/marss_data_sb_vc_072222.rds")


#### Consolidate fire effect to one col ####

dat2 = dat1

# first, change effect to just the month of ignition (in contrast to 1 year long effect)
# fires = 
  # AB00_Jesusita
  # AB00_Tea
  # AT07_Jesusita
  # GV01_Gaviota
  # HO00_Gaviota 
  # HO00_Sherpa
  # MC06_Tea
  # MC06_Jesusita
  # RG01_Gaviota 
  # RG01_Sherpa
  # RS02_Tea
  # RS02_Jesusita
  # SP02_Gap
  # RED_Thompson
  # EFJ_Thompson
  # EFJ_Conchas
  # RSAW_Thompson
  # RSAW_Conchas
  # RSA_Conchas
  # IND_Conchas
  # IND_BB_Conchas
  # SULF_Thompson

dat2$AB00_Jesusita[dat2$AB00_Jesusita==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
dat2$AB00_Tea[dat2$AB00_Tea==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
dat2$AT07_Jesusita[dat2$AT07_Jesusita==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
dat2$GV01_Gaviota[dat2$GV01_Gaviota==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
dat2$HO00_Gaviota[dat2$HO00_Gaviota==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
#dat2$HO00_Sherpa[dat2$HO00_Sherpa==1] = c(1,0,0,0,0,0,0,0,0,0,0,0) this one only had a 1 m effect because ignition occurred at end of record, so no editing needed
dat2$MC06_Tea[dat2$MC06_Tea==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
dat2$MC06_Jesusita[dat2$MC06_Jesusita==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
dat2$RG01_Gaviota[dat2$RG01_Gaviota==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
# dat2$RG01_Sherpa[dat2$RG01_Sherpa==1] = c(1,0,0,0,0,0,0,0,0,0,0,0) this one only had a 1 m effect because ignition occurred at end of record, so no editing needed
dat2$RS02_Tea[dat2$RS02_Tea==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
dat2$RS02_Jesusita[dat2$RS02_Jesusita==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
dat2$SP02_Gap[dat2$SP02_Gap==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
#
dat2$RED_Thompson[dat2$RED_Thompson==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
dat2$EFJ_Thompson[dat2$EFJ_Thompson==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
dat2$EFJ_Conchas[dat2$EFJ_Conchas==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
dat2$RSAW_Thompson[dat2$RSAW_Thompson==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
dat2$RSAW_Conchas[dat2$RSAW_Conchas==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
dat2$RSA_Conchas[dat2$RSA_Conchas==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
dat2$IND_Conchas[dat2$IND_Conchas==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
dat2$IND_BB_Conchas[dat2$IND_BB_Conchas==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)
dat2$SULF_Thompson[dat2$SULF_Thompson==1] = c(1,0,0,0,0,0,0,0,0,0,0,0)

# Second, combine fires in each watershed so there is one fire effect column

dat2$AB00_allfires = dat2$AB00_Jesusita + dat2$AB00_Tea
dat2$AT07_allfires = dat2$AT07_Jesusita
dat2$GV01_allfires = dat2$GV01_Gaviota
dat2$HO00_allfires = dat2$HO00_Gaviota + dat2$HO00_Sherpa
dat2$MC06_allfires = dat2$MC06_Tea + dat2$MC06_Jesusita
dat2$RG01_allfires = dat2$RG01_Gaviota + dat2$RG01_Sherpa
dat2$RGS02_allfires = dat2$RS02_Tea + dat2$RS02_Jesusita
dat2$SP02_allfires = dat2$SP02_Gap
#
dat2$RED_allfires = dat2$RED_Thompson
dat2$EFJ_allfires = dat2$EFJ_Thompson + dat2$EFJ_Conchas
dat2$RSAW_allfires = dat2$RSAW_Thompson + dat2$RSAW_Conchas
dat2$RSA_allfires = dat2$RSA_Conchas
dat2$IND_allfires = dat2$IND_Conchas
dat2$IND_BB_allfires = dat2$IND_BB_Conchas
dat2$SULF_allfires = dat2$SULF_Thompson

# reorganize
names(dat2)
dat3 = dat2[,c(1,2,3,
               37,38,
               57:71, 
               4,
               34,
               35,36)]

# add all fire effects to 1 col
dat3$fires_pa<-rowSums(dat3[,6:20])

# check that there are no 2s in fire effect cols
max(dat3[,6:20])
max(dat3$fires_pa)
# dat3[,6:20][dat3[,6:20]==2] = 1 replace 2s with 1s if needed

qplot(index, fires_pa, data=dat3, colour=site, geom="path")

#### Add 6 m window to fire effect ####

# fire_pa_6m: 1=fire ignited in prior 6 m, 0=no fire ignitions in 6 m

dat4 = dat3

dat4$date = as.Date(paste(dat4$year, dat4$month, "15", sep="-"))

firedates_6m = c(dat4$date[dat4$fires_pa == 1],
                 dat4$date[dat4$fires_pa == 1] + 30*1,
                 dat4$date[dat4$fires_pa == 1] + 30*2,
                 dat4$date[dat4$fires_pa == 1] + 30*3,
                 dat4$date[dat4$fires_pa == 1] + 30*4,
                 dat4$date[dat4$fires_pa == 1] + 30*5)

firedates_6m = data.frame(year = year(firedates_6m),
                          month = month(firedates_6m),
                          site = rep(dat4$site[dat4$fires_pa == 1], 6),
                          fires_pa_6m = 1)
dat5 = left_join(dat4, firedates_6m, by=c("year","month","site"))
dat5$fires_pa_6m[is.na(dat5$fires_pa_6m)] = 0

qplot(index, fires_pa_6m, data=dat5, colour=site, geom="path")

#### Add fire p/a legacy effects ####

# 1 year legacy
firedates_6m_1ylegacy = data.frame(year = firedates_6m$year+1,
                                   month = firedates_6m$month,
                                   site = firedates_6m$site,
                                   fire_pa_6m_1ylegacy = 1)
dat6 = left_join(dat5, firedates_6m_1ylegacy, by=c("year","month","site"))
dat6$fire_pa_6m_1ylegacy[is.na(dat6$fire_pa_6m_1ylegacy)] = 0

qplot(index, fire_pa_6m_1ylegacy, data=dat6, colour=site, geom="path")

# 2 year legacy
firedates_6m_2ylegacy = data.frame(year = firedates_6m$year+2,
                                   month = firedates_6m$month,
                                   site = firedates_6m$site,
                                   fire_pa_6m_2ylegacy = 1)
dat7 = left_join(dat6, firedates_6m_2ylegacy, by=c("year","month","site"))

# 3 year legacy
firedates_6m_3ylegacy = data.frame(year = firedates_6m$year+3,
                                   month = firedates_6m$month,
                                   site = firedates_6m$site,
                                   fire_pa_6m_3ylegacy = 1)
dat8 = left_join(dat7, firedates_6m_3ylegacy, by=c("year","month","site"))

# 4 year legacy
firedates_6m_4ylegacy = data.frame(year = firedates_6m$year+4,
                                   month = firedates_6m$month,
                                   site = firedates_6m$site,
                                   fire_pa_6m_4ylegacy = 1)
dat9 = left_join(dat8, firedates_6m_4ylegacy, by=c("year","month","site"))

# 5 year legacy
firedates_6m_5ylegacy = data.frame(year = firedates_6m$year+5,
                                   month = firedates_6m$month,
                                   site = firedates_6m$site,
                                   fire_pa_6m_5ylegacy = 1)
dat10 = left_join(dat9, firedates_6m_5ylegacy, by=c("year","month","site"))

# replace NAs with 0
dat10[,28:32][is.na(dat10[,28:32])] = 0

#### Add fire p/a x ppt legacy effects ####

dat10$fire_pa_6m_ppt = dat10$fires_pa_6m*dat10$cumulative_precip_mm
dat10$fire_pa_6m_ppt_1ylegacy = dat10$fire_pa_6m_1ylegacy*dat10$cumulative_precip_mm
dat10$fire_pa_6m_ppt_2ylegacy = dat10$fire_pa_6m_2ylegacy*dat10$cumulative_precip_mm
dat10$fire_pa_6m_ppt_3ylegacy = dat10$fire_pa_6m_3ylegacy*dat10$cumulative_precip_mm
dat10$fire_pa_6m_ppt_4ylegacy = dat10$fire_pa_6m_4ylegacy*dat10$cumulative_precip_mm
dat10$fire_pa_6m_ppt_5ylegacy = dat10$fire_pa_6m_5ylegacy*dat10$cumulative_precip_mm

p1 = qplot(index, cumulative_precip_mm, data=dat10, colour=site, geom="path")
p2 = qplot(index, fire_pa_6m_ppt, data=dat10, colour=site, geom="path")
p3 = qplot(index, fires_pa_6m, data=dat10, colour=site, geom="path")
gridExtra::grid.arrange(p1,p3,p2, ncol=1)

#### Export data with fire x ppt interactions and legacy effects ####

saveRDS(dat10, "data_working/marss_data_sb_vc_072222.rds")

#### MARSS: ppt, fire_6m x ppt, no legacy effects ####
#### Set up data for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_072222.rds")

# order in time and reset index to assure an ordered ts for MARSS
dat = dat[with(dat, order(region, site, date)), ]
table(dat$site)
dat$index = rep(seq(1,166), 15)

# select sites
# include these sites only (11 total - these have the longest most complete ts for SpC):
# AB00, AT07, GV01, HO00, MC06, RG01, RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "AT07", "GV01", "HO00", "MC06", "RG01", "RS02", 
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, mean_cond_uScm, cumulative_precip_mm, 
    fires_pa_6m,
    fire_pa_6m_ppt) %>% 
  pivot_wider(
    names_from = site, values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                                       fires_pa_6m, 
                                       fire_pa_6m_ppt)) 

dat_cond[is.nan(dat_cond)] <- NA

# log and scale transform response var
names(dat_cond)
dat_cond_log = dat_cond
dat_cond_log[,2:12] = log10(dat_cond_log[,2:12])
dat_cond_log[,2:12] = scale(dat_cond_log[,2:12])
sum(is.nan(dat_cond_log[,2:12]))
sum(is.na(dat_cond_log[,2:12]))
range(dat_cond_log[,2:12], na.rm = T)

# Pull out only response var
names(dat_cond_log)
dat_dep <- t(dat_cond_log[,c(2:12)])
row.names(dat_dep)

# check covars for nas, nans, or infs
sum(is.nan(dat_cond_log[,13:45]))
sum(is.na(dat_cond_log[,13:45]))
sum(is.infinite(dat_cond_log[,13:45]))

# Make covariate inputs
dat_cov <- dat_cond_log[,c(13:45)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov))
sum(is.na(dat_cov))
sum(is.infinite(dat_cov))

#### make C matrix  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CC <- matrix(list( 
  # precip by site
  "cumulative_precip_mm_AB00",0,0,0,0,0,0,0,0,0,0,
  0,"cumulative_precip_mm_AT07",0,0,0,0,0,0,0,0,0,
  0,0,"cumulative_precip_mm_GV01",0,0,0,0,0,0,0,0,
  0,0,0,"cumulative_precip_mm_HO00",0,0,0,0,0,0,0,
  0,0,0,0,"cumulative_precip_mm_MC06",0,0,0,0,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RG01",0,0,0,0,0,
  0,0,0,0,0,0,"cumulative_precip_mm_RS02",0,0,0,0,
  0,0,0,0,0,0,0,"cumulative_precip_mm_EFJ", 0,0,0,
  0,0,0,0,0,0,0,0,"cumulative_precip_mm_RED", 0,0,
  0,0,0,0,0,0,0,0,0,"cumulative_precip_mm_RSA", 0,
  0,0,0,0,0,0,0,0,0,0,"cumulative_precip_mm_RSAW",
  # fires p/a in 6 m window
  "fires_pa_6m_AB00",0,0,0,0,0,0,0,0,0,0,
  0,"fires_pa_6m_AT07",0,0,0,0,0,0,0,0,0,
  0,0,"fires_pa_6m_GV01",0,0,0,0,0,0,0,0,
  0,0,0,"fires_pa_6m_HO00",0,0,0,0,0,0,0,
  0,0,0,0,"fires_pa_6m_MC06",0,0,0,0,0,0,
  0,0,0,0,0,"fires_pa_6m_RG01",0,0,0,0,0,
  0,0,0,0,0,0,"fires_pa_6m_RS02",0,0,0,0,
  0,0,0,0,0,0,0,"fires_pa_6m_EFJ", 0,0,0,
  0,0,0,0,0,0,0,0,"fires_pa_6m_RED", 0,0,
  0,0,0,0,0,0,0,0,0,"fires_pa_6m_RSA", 0,
  0,0,0,0,0,0,0,0,0,0,"fires_pa_6m_RSAW",
  # interaction of cum. ppt with fire p/a in 6 m window
  "fire_pa_6m_ppt_AB00",0,0,0,0,0,0,0,0,0,0,
  0,"fire_pa_6m_ppt_AT07",0,0,0,0,0,0,0,0,0,
  0,0,"fire_pa_6m_ppt_GV01",0,0,0,0,0,0,0,0,
  0,0,0,"fire_pa_6m_ppt_HO00",0,0,0,0,0,0,0,
  0,0,0,0,"fire_pa_6m_ppt_MC06",0,0,0,0,0,0,
  0,0,0,0,0,"fire_pa_6m_ppt_RG01",0,0,0,0,0,
  0,0,0,0,0,0,"fire_pa_6m_ppt_RS02",0,0,0,0,
  0,0,0,0,0,0,0,"fire_pa_6m_ppt_EFJ", 0,0,0,
  0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_RED", 0,0,
  0,0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_RSA", 0,
  0,0,0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_RSAW"
  ),
  11,33)

#### Model setup for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++

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
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

#### Fit MARSS model +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_07232022_11state_cond_firexppt_mBFGS.rds")

#### DIAGNOSES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Plot coef and coef estimates ###
## estimates
# hessian method is fast but not ideal for final results - bootstrap for final results
est_fit <- MARSSparamCIs(fit)

# formatting confidence intervals into dataframe
CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)

### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00","AT07","GV01","HO00","MC06","RG01","RS02","EFJ","RED","RSA","RSAW")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  dplyr::rename(Site = i) %>%
  rename(Parameter = parm) # rename site column

# plot results
(RESULTS_ALL_d <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 07/23/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(.~Site, scales = "free"))

#### MARSS: ppt, fire_6m x ppt, no legacy effects, no HO00 ####

#### Set up data for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_072222.rds")

# # order in time and reset index to assure an ordered ts for MARSS
# dat = dat[with(dat, order(region, site, date)), ]
# table(dat$site)
# dat$index = rep(seq(1,166), 15)

# select sites
# include these sites only (10 total - these have the longest most complete ts for SpC and I have also removed HO00 because it was causing issues with missing data):
# AB00, AT07, GV01, , MC06, RG01, RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "AT07", "GV01", "MC06", "RG01", "RS02", 
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, mean_cond_uScm, cumulative_precip_mm, fires_pa_6m, fire_pa_6m_ppt) %>% 
  pivot_wider(
    names_from = site, values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                                       fires_pa_6m, 
                                       fire_pa_6m_ppt)) 

dat_cond[is.nan(dat_cond)] <- NA

# log and scale transform response var
names(dat_cond)
dat_cond_log = dat_cond
#dat_cond_log[,2:11] = log10(dat_cond_log[,2:11])
dat_cond_log[,2:11] = scale(dat_cond_log[,2:11])
sum(is.nan(dat_cond_log[,2:11]))
sum(is.na(dat_cond_log[,2:11]))
range(dat_cond_log[,2:11], na.rm = T)

# Pull out only response var
names(dat_cond_log)
dat_dep <- t(dat_cond_log[,c(2:11)])
row.names(dat_dep)

# check covars for nas, nans, or infs
sum(is.nan(dat_cond_log[,12:41]))
sum(is.na(dat_cond_log[,12:41]))
sum(is.infinite(dat_cond_log[,12:41]))

# Make covariate inputs
dat_cov <- dat_cond_log[,c(12:41)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov))
sum(is.na(dat_cov))
sum(is.infinite(dat_cov))
# # fire_pa_6m_ppt_1ylegacy_HO00 has no values to scale, so must replace with zeros
# dat_cov["fire_pa_6m_ppt_1ylegacy_HO00",]= 0.001

#### make C matrix  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# "XXXX_AB00",0,0,0,0,0,0,0,0,0,
# 0,"XXXX_AT07",0,0,0,0,0,0,0,0,
# 0,0,"XXXX_GV01",0,0,0,0,0,0,0,
# 0,0,0,"XXXX_MC06",0,0,0,0,0,0,
# 0,0,0,0,"XXXX_RG01",0,0,0,0,0,
# 0,0,0,0,0,"XXXX_RS02",0,0,0,0,
# 0,0,0,0,0,0,"XXXX_EFJ", 0,0,0,
# 0,0,0,0,0,0,0,"XXXX_RED", 0,0,
# 0,0,0,0,0,0,0,0,"XXXX_RSA", 0,
# 0,0,0,0,0,0,0,0,0,"XXXX_RSAW",

CC <- matrix(list( 
  # precip by site
  "cumulative_precip_mm_AB00",0,0,0,0,0,0,0,0,0,
  0,"cumulative_precip_mm_AT07",0,0,0,0,0,0,0,0,
  0,0,"cumulative_precip_mm_GV01",0,0,0,0,0,0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,0,0,0,0,0,
  0,0,0,0,"cumulative_precip_mm_RG01",0,0,0,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RS02",0,0,0,0,
  0,0,0,0,0,0,"cumulative_precip_mm_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"cumulative_precip_mm_RED", 0,0,
  0,0,0,0,0,0,0,0,"cumulative_precip_mm_RSA", 0,
  0,0,0,0,0,0,0,0,0,"cumulative_precip_mm_RSAW",
  # fires p/a in 6 m window
  "fires_pa_6m_AB00",0,0,0,0,0,0,0,0,0,
  0,"fires_pa_6m_AT07",0,0,0,0,0,0,0,0,
  0,0,"fires_pa_6m_GV01",0,0,0,0,0,0,0,
  0,0,0,"fires_pa_6m_MC06",0,0,0,0,0,0,
  0,0,0,0,"fires_pa_6m_RG01",0,0,0,0,0,
  0,0,0,0,0,"fires_pa_6m_RS02",0,0,0,0,
  0,0,0,0,0,0,"fires_pa_6m_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"fires_pa_6m_RED", 0,0,
  0,0,0,0,0,0,0,0,"fires_pa_6m_RSA", 0,
  0,0,0,0,0,0,0,0,0,"fires_pa_6m_RSAW",
  # interaction of cum. ppt with fire p/a in 6 m window
  "fire_pa_6m_ppt_AB00",0,0,0,0,0,0,0,0,0,
  0,"fire_pa_6m_ppt_AT07",0,0,0,0,0,0,0,0,
  0,0,"fire_pa_6m_ppt_GV01",0,0,0,0,0,0,0,
  0,0,0,"fire_pa_6m_ppt_MC06",0,0,0,0,0,0,
  0,0,0,0,"fire_pa_6m_ppt_RG01",0,0,0,0,0,
  0,0,0,0,0,"fire_pa_6m_ppt_RS02",0,0,0,0,
  0,0,0,0,0,0,"fire_pa_6m_ppt_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"fire_pa_6m_ppt_RED", 0,0,
  0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_RSA", 0,
  0,0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_RSAW"
),
10,30)

#### Model setup for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++

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
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

#### Fit MARSS model +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_07232022_10state_cond_firexppt_mBFGS.rds")

#### DIAGNOSES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)

# formatting confidence intervals into dataframe
CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)

### Plot Results for All Sites ###

# First, create dataset of all outputs
# RSA and RSAW get mixed up in plotting, so replacing RSAW with SAW
CIs_fit$parm =gsub("RSAW","SAW",CIs_fit$parm)
rownames(CIs_fit) =gsub("RSAW","SAW",rownames(CIs_fit))

# This works for AB00 alone
CIs_AB00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("AB00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00","AT07","GV01","MC06","RG01","RS02","EFJ","RED","RSA","SAW")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  dplyr::rename(Site = i) %>%
  rename(Parameter = parm) # rename site column
CIs_fit_ed$Parameter = rep(c("Cum. Ppt", "Fire p/a (6m window)","Cum. Ppt * Fire p/a"), 10)
CIs_fit_ed$Region = c(rep(c("SB"),6*3), rep(c("VC"),4*3))

# plot results
(RESULTS_ALL_d <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 07/23/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

### plots of input data ###

EFJ = dat[dat$site=="EFJ",]

par(mfrow=c(2,1))
plot(EFJ$cumulative_precip_mm~EFJ$date , ylim = rev(range(EFJ$cumulative_precip_mm)), type="b")
plot(EFJ$mean_cond_uScm~EFJ$date)
#
EFJ_2011 = EFJ[EFJ$year==2011,]
plot(EFJ_2011$cumulative_precip_mm~EFJ_2011$date , ylim = rev(range(EFJ_2011$cumulative_precip_mm)), type="b")
plot(EFJ_2011$mean_cond_uScm~EFJ_2011$date)
#
EFJ_2013 = EFJ[EFJ$year==2013,]
plot(EFJ_2013$cumulative_precip_mm~EFJ_2013$index , ylim = rev(range(EFJ_2013$cumulative_precip_mm)), type="b")
plot(EFJ_2013$mean_cond_uScm~EFJ_2013$index)

#### MARSS: ppt, fire x ppt, no legacy effects, no HO00 - BEST CURRENT MODEL ####

#### Set up data for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_072222.rds")

# # order in time and reset index to assure an ordered ts for MARSS
# dat = dat[with(dat, order(region, site, date)), ]
# table(dat$site)
# dat$index = rep(seq(1,166), 15)

# select sites
# include these sites only (10 total - these have the longest most complete ts for SpC and I have also removed HO00 because it was causing issues with missing data):
# AB00, AT07, GV01, , MC06, RG01, RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "AT07", "GV01", "MC06", "RG01", "RS02", 
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, mean_cond_uScm, cumulative_precip_mm, fires_pa, fire_pa_6m_ppt) %>% 
  pivot_wider(
    names_from = site, values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                                       fires_pa, 
                                       fire_pa_6m_ppt)) 

dat_cond[is.nan(dat_cond)] <- NA

# log and scale transform response var
names(dat_cond)
dat_cond_log = dat_cond
#dat_cond_log[,2:11] = log10(dat_cond_log[,2:11])
dat_cond_log[,2:11] = scale(dat_cond_log[,2:11])
sum(is.nan(dat_cond_log[,2:11]))
sum(is.na(dat_cond_log[,2:11]))
range(dat_cond_log[,2:11], na.rm = T)

# Pull out only response var
names(dat_cond_log)
dat_dep <- t(dat_cond_log[,c(2:11)])
row.names(dat_dep)

# check covars for nas, nans, or infs
sum(is.nan(dat_cond_log[,12:41]))
sum(is.na(dat_cond_log[,12:41]))
sum(is.infinite(dat_cond_log[,12:41]))

# Make covariate inputs
dat_cov <- dat_cond_log[,c(12:41)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov))
sum(is.na(dat_cov))
sum(is.infinite(dat_cov))
# # fire_pa_6m_ppt_1ylegacy_HO00 has no values to scale, so must replace with zeros
# dat_cov["fire_pa_6m_ppt_1ylegacy_HO00",]= 0.001

#### make C matrix  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# "XXXX_AB00",0,0,0,0,0,0,0,0,0,
# 0,"XXXX_AT07",0,0,0,0,0,0,0,0,
# 0,0,"XXXX_GV01",0,0,0,0,0,0,0,
# 0,0,0,"XXXX_MC06",0,0,0,0,0,0,
# 0,0,0,0,"XXXX_RG01",0,0,0,0,0,
# 0,0,0,0,0,"XXXX_RS02",0,0,0,0,
# 0,0,0,0,0,0,"XXXX_EFJ", 0,0,0,
# 0,0,0,0,0,0,0,"XXXX_RED", 0,0,
# 0,0,0,0,0,0,0,0,"XXXX_RSA", 0,
# 0,0,0,0,0,0,0,0,0,"XXXX_RSAW",

CC <- matrix(list( 
  # precip by site
  "cumulative_precip_mm_AB00",0,0,0,0,0,0,0,0,0,
  0,"cumulative_precip_mm_AT07",0,0,0,0,0,0,0,0,
  0,0,"cumulative_precip_mm_GV01",0,0,0,0,0,0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,0,0,0,0,0,
  0,0,0,0,"cumulative_precip_mm_RG01",0,0,0,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RS02",0,0,0,0,
  0,0,0,0,0,0,"cumulative_precip_mm_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"cumulative_precip_mm_RED", 0,0,
  0,0,0,0,0,0,0,0,"cumulative_precip_mm_RSA", 0,
  0,0,0,0,0,0,0,0,0,"cumulative_precip_mm_RSAW",
  # fires p/a
  "fires_pa_AB00",0,0,0,0,0,0,0,0,0,
  0,"fires_pa_AT07",0,0,0,0,0,0,0,0,
  0,0,"fires_pa_GV01",0,0,0,0,0,0,0,
  0,0,0,"fires_pa_MC06",0,0,0,0,0,0,
  0,0,0,0,"fires_pa_RG01",0,0,0,0,0,
  0,0,0,0,0,"fires_pa_RS02",0,0,0,0,
  0,0,0,0,0,0,"fires_pa_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"fires_pa_RED", 0,0,
  0,0,0,0,0,0,0,0,"fires_pa_RSA", 0,
  0,0,0,0,0,0,0,0,0,"fires_pa_RSAW",
  # interaction of cum. ppt with fire p/a in 6 m window
  "fire_pa_6m_ppt_AB00",0,0,0,0,0,0,0,0,0,
  0,"fire_pa_6m_ppt_AT07",0,0,0,0,0,0,0,0,
  0,0,"fire_pa_6m_ppt_GV01",0,0,0,0,0,0,0,
  0,0,0,"fire_pa_6m_ppt_MC06",0,0,0,0,0,0,
  0,0,0,0,"fire_pa_6m_ppt_RG01",0,0,0,0,0,
  0,0,0,0,0,"fire_pa_6m_ppt_RS02",0,0,0,0,
  0,0,0,0,0,0,"fire_pa_6m_ppt_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"fire_pa_6m_ppt_RED", 0,0,
  0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_RSA", 0,
  0,0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_RSAW"
),
10,30)

#### Model setup for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++

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
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

#### Fit MARSS model +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

#### DIAGNOSES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)

# formatting confidence intervals into dataframe
CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)

### Plot Results for All Sites ###

# First, create dataset of all outputs
# RSA and RSAW get mixed up in plotting, so replacing RSAW with SAW
CIs_fit$parm =gsub("RSAW","SAW",CIs_fit$parm)
rownames(CIs_fit) =gsub("RSAW","SAW",rownames(CIs_fit))

# This works for AB00 alone
CIs_AB00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("AB00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00","AT07","GV01","MC06","RG01","RS02","EFJ","RED","RSA","SAW")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  dplyr::rename(Site = i) %>%
  rename(Parameter = parm) # rename site column
CIs_fit_ed$Parameter = rep(c("Cum. Ppt", "Fire p/a","Cum. Ppt * Fire p/a (6m)"), 10)
CIs_fit_ed$Region = c(rep(c("SB"),6*3), rep(c("VC"),4*3))

# plot results
(RESULTS_ALL_d <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 07/23/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))


#
#### MARSS: ppt, fire x ppt, 5 years of legacy effects ####
# notes: trouble with convergence

#### Set up data for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_072222.rds")

# order in time and reset index to assure an ordered ts for MARSS
dat = dat[with(dat, order(region, site, date)), ]
table(dat$site)
dat$index = rep(seq(1,166), 15)

# select sites
# include these sites only (10 total - these have the longest most complete ts for SpC and I have also removed HO00 because it was causing issues with missing data):
# AB00, AT07, GV01, , MC06, RG01, RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "AT07", "GV01", "MC06", "RG01", "RS02", 
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, mean_cond_uScm, cumulative_precip_mm, 
    fires_pa_6m, fire_pa_6m_1ylegacy, fire_pa_6m_2ylegacy, fire_pa_6m_3ylegacy, fire_pa_6m_4ylegacy, fire_pa_6m_5ylegacy,
    fire_pa_6m_ppt, fire_pa_6m_ppt_1ylegacy, fire_pa_6m_ppt_2ylegacy, fire_pa_6m_ppt_3ylegacy, fire_pa_6m_ppt_4ylegacy, fire_pa_6m_ppt_5ylegacy) %>% 
  pivot_wider(
    names_from = site, values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                                       fires_pa_6m, fire_pa_6m_1ylegacy, fire_pa_6m_2ylegacy, fire_pa_6m_3ylegacy, fire_pa_6m_4ylegacy, fire_pa_6m_5ylegacy,
                                       fire_pa_6m_ppt, fire_pa_6m_ppt_1ylegacy, fire_pa_6m_ppt_2ylegacy, fire_pa_6m_ppt_3ylegacy, fire_pa_6m_ppt_4ylegacy, fire_pa_6m_ppt_5ylegacy)) 

dat_cond[is.nan(dat_cond)] <- NA

# log and scale transform response var
names(dat_cond)
dat_cond_log = dat_cond
dat_cond_log[,2:11] = log10(dat_cond_log[,2:11])
dat_cond_log[,2:11] = scale(dat_cond_log[,2:11])
sum(is.nan(dat_cond_log[,2:11]))
sum(is.na(dat_cond_log[,2:11]))
range(dat_cond_log[,2:11], na.rm = T)

# Pull out only response var
names(dat_cond_log)
dat_dep <- t(dat_cond_log[,c(2:11)])
row.names(dat_dep)

# check covars for nas, nans, or infs
sum(is.nan(dat_cond_log[,12:141]))
sum(is.na(dat_cond_log[,12:141]))
sum(is.infinite(dat_cond_log[,12:141]))

# Make covariate inputs
dat_cov <- dat_cond_log[,c(12:141)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov))
sum(is.na(dat_cov))
sum(is.infinite(dat_cov))
# # fire_pa_6m_ppt_1ylegacy_HO00 has no values to scale, so must replace with zeros
# dat_cov["fire_pa_6m_ppt_1ylegacy_HO00",]= 0.001

#### make C matrix  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# "XXXX_AB00",0,0,0,0,0,0,0,0,0,
# 0,"XXXX_AT07",0,0,0,0,0,0,0,0,
# 0,0,"XXXX_GV01",0,0,0,0,0,0,0,
# 0,0,0,"XXXX_MC06",0,0,0,0,0,0,
# 0,0,0,0,"XXXX_RG01",0,0,0,0,0,
# 0,0,0,0,0,"XXXX_RS02",0,0,0,0,
# 0,0,0,0,0,0,"XXXX_EFJ", 0,0,0,
# 0,0,0,0,0,0,0,"XXXX_RED", 0,0,
# 0,0,0,0,0,0,0,0,"XXXX_RSA", 0,
# 0,0,0,0,0,0,0,0,0,"XXXX_RSAW",

CC <- matrix(list( 
  # precip by site
  "cumulative_precip_mm_AB00",0,0,0,0,0,0,0,0,0,
  0,"cumulative_precip_mm_AT07",0,0,0,0,0,0,0,0,
  0,0,"cumulative_precip_mm_GV01",0,0,0,0,0,0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,0,0,0,0,0,
  0,0,0,0,"cumulative_precip_mm_RG01",0,0,0,0,0,
  0,0,0,0,0,"cumulative_precip_mm_RS02",0,0,0,0,
  0,0,0,0,0,0,"cumulative_precip_mm_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"cumulative_precip_mm_RED", 0,0,
  0,0,0,0,0,0,0,0,"cumulative_precip_mm_RSA", 0,
  0,0,0,0,0,0,0,0,0,"cumulative_precip_mm_RSAW",
  # fires p/a in 6 m window
  "fires_pa_6m_AB00",0,0,0,0,0,0,0,0,0,
  0,"fires_pa_6m_AT07",0,0,0,0,0,0,0,0,
  0,0,"fires_pa_6m_GV01",0,0,0,0,0,0,0,
  0,0,0,"fires_pa_6m_MC06",0,0,0,0,0,0,
  0,0,0,0,"fires_pa_6m_RG01",0,0,0,0,0,
  0,0,0,0,0,"fires_pa_6m_RS02",0,0,0,0,
  0,0,0,0,0,0,"fires_pa_6m_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"fires_pa_6m_RED", 0,0,
  0,0,0,0,0,0,0,0,"fires_pa_6m_RSA", 0,
  0,0,0,0,0,0,0,0,0,"fires_pa_6m_RSAW",
  # 1 y legacy of fires p/a in 6 m window: fires_pa_6m_1ylegacy
  "fires_pa_6m_1ylegacy_AB00",0,0,0,0,0,0,0,0,0,
  0,"fires_pa_6m_1ylegacy_AT07",0,0,0,0,0,0,0,0,
  0,0,"fires_pa_6m_1ylegacy_GV01",0,0,0,0,0,0,0,
  0,0,0,"fires_pa_6m_1ylegacy_MC06",0,0,0,0,0,0,
  0,0,0,0,"fires_pa_6m_1ylegacy_RG01",0,0,0,0,0,
  0,0,0,0,0,"fires_pa_6m_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,0,0,"fires_pa_6m_1ylegacy_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"fires_pa_6m_1ylegacy_RED", 0,0,
  0,0,0,0,0,0,0,0,"fires_pa_6m_1ylegacy_RSA", 0,
  0,0,0,0,0,0,0,0,0,"fires_pa_6m_1ylegacy_RSAW",
  # 2 y legacy of fires p/a in 6 m window: fires_pa_6m_2ylegacy
  "fires_pa_6m_2ylegacy_AB00",0,0,0,0,0,0,0,0,0,
  0,"fires_pa_6m_2ylegacy_AT07",0,0,0,0,0,0,0,0,
  0,0,"fires_pa_6m_2ylegacy_GV01",0,0,0,0,0,0,0,
  0,0,0,"fires_pa_6m_2ylegacy_MC06",0,0,0,0,0,0,
  0,0,0,0,"fires_pa_6m_2ylegacy_RG01",0,0,0,0,0,
  0,0,0,0,0,"fires_pa_6m_2ylegacy_RS02",0,0,0,0,
  0,0,0,0,0,0,"fires_pa_6m_2ylegacy_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"fires_pa_6m_2ylegacy_RED", 0,0,
  0,0,0,0,0,0,0,0,"fires_pa_6m_2ylegacy_RSA", 0,
  0,0,0,0,0,0,0,0,0,"fires_pa_6m_2ylegacy_RSAW",
  # 3y legacy of fires p/a in 6 m window
  "fires_pa_6m_3ylegacy_AB00",0,0,0,0,0,0,0,0,0,
  0,"fires_pa_6m_3ylegacy_AT07",0,0,0,0,0,0,0,0,
  0,0,"fires_pa_6m_3ylegacy_GV01",0,0,0,0,0,0,0,
  0,0,0,"fires_pa_6m_3ylegacy_MC06",0,0,0,0,0,0,
  0,0,0,0,"fires_pa_6m_3ylegacy_RG01",0,0,0,0,0,
  0,0,0,0,0,"fires_pa_6m_3ylegacy_RS02",0,0,0,0,
  0,0,0,0,0,0,"fires_pa_6m_3ylegacy_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"fires_pa_6m_3ylegacy_RED", 0,0,
  0,0,0,0,0,0,0,0,"fires_pa_6m_3ylegacy_RSA", 0,
  0,0,0,0,0,0,0,0,0,"fires_pa_6m_3ylegacy_RSAW",
  # 4y legacy of fires p/a in 6 m window
  "fires_pa_6m_3ylegacy_AB00",0,0,0,0,0,0,0,0,0,
  0,"fires_pa_6m_3ylegacy_AT07",0,0,0,0,0,0,0,0,
  0,0,"fires_pa_6m_3ylegacy_GV01",0,0,0,0,0,0,0,
  0,0,0,"fires_pa_6m_3ylegacy_MC06",0,0,0,0,0,0,
  0,0,0,0,"fires_pa_6m_3ylegacy_RG01",0,0,0,0,0,
  0,0,0,0,0,"fires_pa_6m_3ylegacy_RS02",0,0,0,0,
  0,0,0,0,0,0,"fires_pa_6m_3ylegacy_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"fires_pa_6m_3ylegacy_RED", 0,0,
  0,0,0,0,0,0,0,0,"fires_pa_6m_3ylegacy_RSA", 0,
  0,0,0,0,0,0,0,0,0,"fires_pa_6m_3ylegacy_RSAW",
  # 5y legacy of fires p/a in 6 m window
  "fires_pa_6m_5ylegacy_AB00",0,0,0,0,0,0,0,0,0,
  0,"fires_pa_6m_5ylegacy_AT07",0,0,0,0,0,0,0,0,
  0,0,"fires_pa_6m_5ylegacy_GV01",0,0,0,0,0,0,0,
  0,0,0,"fires_pa_6m_5ylegacy_MC06",0,0,0,0,0,0,
  0,0,0,0,"fires_pa_6m_5ylegacy_RG01",0,0,0,0,0,
  0,0,0,0,0,"fires_pa_6m_5ylegacy_RS02",0,0,0,0,
  0,0,0,0,0,0,"fires_pa_6m_5ylegacy_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"fires_pa_6m_5ylegacy_RED", 0,0,
  0,0,0,0,0,0,0,0,"fires_pa_6m_5ylegacy_RSA", 0,
  0,0,0,0,0,0,0,0,0,"fires_pa_6m_5ylegacy_RSAW",
  # interaction of cum. ppt with fire p/a in 6 m window
  "fire_pa_6m_ppt_AB00",0,0,0,0,0,0,0,0,0,
  0,"fire_pa_6m_ppt_AT07",0,0,0,0,0,0,0,0,
  0,0,"fire_pa_6m_ppt_GV01",0,0,0,0,0,0,0,
  0,0,0,"fire_pa_6m_ppt_MC06",0,0,0,0,0,0,
  0,0,0,0,"fire_pa_6m_ppt_RG01",0,0,0,0,0,
  0,0,0,0,0,"fire_pa_6m_ppt_RS02",0,0,0,0,
  0,0,0,0,0,0,"fire_pa_6m_ppt_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"fire_pa_6m_ppt_RED", 0,0,
  0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_RSA", 0,
  0,0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_RSAW",
  # 1y legacy of interaction of cum. ppt with fire p/a in 6 m window
  "fire_pa_6m_ppt_1ylegacy_AB00",0,0,0,0,0,0,0,0,0,
  0,"fire_pa_6m_ppt_1ylegacy_AT07",0,0,0,0,0,0,0,0,
  0,0,"fire_pa_6m_ppt_1ylegacy_GV01",0,0,0,0,0,0,0,
  0,0,0,"fire_pa_6m_ppt_1ylegacy_MC06",0,0,0,0,0,0,
  0,0,0,0,"fire_pa_6m_ppt_1ylegacy_RG01",0,0,0,0,0,
  0,0,0,0,0,"fire_pa_6m_ppt_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,0,0,"fire_pa_6m_ppt_1ylegacy_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"fire_pa_6m_ppt_1ylegacy_RED", 0,0,
  0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_1ylegacy_RSA", 0,
  0,0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_1ylegacy_RSAW",
  # 2y legacy of interaction of cum. ppt with fire p/a in 6 m window
  "fire_pa_6m_ppt_2ylegacy_AB00",0,0,0,0,0,0,0,0,0,
  0,"fire_pa_6m_ppt_2ylegacy_AT07",0,0,0,0,0,0,0,0,
  0,0,"fire_pa_6m_ppt_2ylegacy_GV01",0,0,0,0,0,0,0,
  0,0,0,"fire_pa_6m_ppt_2ylegacy_MC06",0,0,0,0,0,0,
  0,0,0,0,"fire_pa_6m_ppt_2ylegacy_RG01",0,0,0,0,0,
  0,0,0,0,0,"fire_pa_6m_ppt_2ylegacy_RS02",0,0,0,0,
  0,0,0,0,0,0,"fire_pa_6m_ppt_2ylegacy_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"fire_pa_6m_ppt_2ylegacy_RED", 0,0,
  0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_2ylegacy_RSA", 0,
  0,0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_2ylegacy_RSAW",
  # 3y legacy of interaction of cum. ppt with fire p/a in 6 m window
  "fire_pa_6m_ppt_3ylegacy_AB00",0,0,0,0,0,0,0,0,0,
  0,"fire_pa_6m_ppt_3ylegacy_AT07",0,0,0,0,0,0,0,0,
  0,0,"fire_pa_6m_ppt_3ylegacy_GV01",0,0,0,0,0,0,0,
  0,0,0,"fire_pa_6m_ppt_3ylegacy_MC06",0,0,0,0,0,0,
  0,0,0,0,"fire_pa_6m_ppt_3ylegacy_RG01",0,0,0,0,0,
  0,0,0,0,0,"fire_pa_6m_ppt_3ylegacy_RS02",0,0,0,0,
  0,0,0,0,0,0,"fire_pa_6m_ppt_3ylegacy_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"fire_pa_6m_ppt_3ylegacy_RED", 0,0,
  0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_3ylegacy_RSA", 0,
  0,0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_3ylegacy_RSAW",
  # 4y legacy of interaction of cum. ppt with fire p/a in 6 m window
  "fire_pa_6m_ppt_4ylegacy_AB00",0,0,0,0,0,0,0,0,0,
  0,"fire_pa_6m_ppt_4ylegacy_AT07",0,0,0,0,0,0,0,0,
  0,0,"fire_pa_6m_ppt_4ylegacy_GV01",0,0,0,0,0,0,0,
  0,0,0,"fire_pa_6m_ppt_4ylegacy_MC06",0,0,0,0,0,0,
  0,0,0,0,"fire_pa_6m_ppt_4ylegacy_RG01",0,0,0,0,0,
  0,0,0,0,0,"fire_pa_6m_ppt_4ylegacy_RS02",0,0,0,0,
  0,0,0,0,0,0,"fire_pa_6m_ppt_4ylegacy_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"fire_pa_6m_ppt_4ylegacy_RED", 0,0,
  0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_4ylegacy_RSA", 0,
  0,0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_4ylegacy_RSAW",
  # 5y legacy of interaction of cum. ppt with fire p/a in 6 m window
  "fire_pa_6m_ppt_5ylegacy_AB00",0,0,0,0,0,0,0,0,0,
  0,"fire_pa_6m_ppt_5ylegacy_AT07",0,0,0,0,0,0,0,0,
  0,0,"fire_pa_6m_ppt_5ylegacy_GV01",0,0,0,0,0,0,0,
  0,0,0,"fire_pa_6m_ppt_5ylegacy_MC06",0,0,0,0,0,0,
  0,0,0,0,"fire_pa_6m_ppt_5ylegacy_RG01",0,0,0,0,0,
  0,0,0,0,0,"fire_pa_6m_ppt_5ylegacy_RS02",0,0,0,0,
  0,0,0,0,0,0,"fire_pa_6m_ppt_5ylegacy_EFJ", 0,0,0,
  0,0,0,0,0,0,0,"fire_pa_6m_ppt_5ylegacy_RED", 0,0,
  0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_5ylegacy_RSA", 0,
  0,0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_5ylegacy_RSAW"
),
10,130)

#### Model setup for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++

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
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

#### Fit MARSS model +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_07232022_10state_cond_firexppt&legacies_mBFGS.rds")

#### DIAGNOSES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
# Warning messages:
# 1: In MARSShessian(MLEobj, method = hessian.fun) :
#   MARSShessian: Hessian could not be inverted to compute the parameter var-cov matrix. parSigma set to NULL.  See MARSSinfo("HessianNA").
# 
# 2: In MARSSparamCIs(fit) :
#   MARSSparamCIs: No parSigma element returned by Hessian function.  See marssMLE object errors (MLEobj$errors)
# > MARSSinfo("HessianNA")
# The variance-covariance matrix can be estimated (large sample estimator) from the
# inverse of the Hessian of the log-likelihood function at the MLE parameter values.
# The Hessian is the second partial derivative of a matrix function. The Hessian of
# the log-likelihood function at the MLEs is the observed Fisher information. The
# observed Fisher information is an estimator of large-sample variance-covariance
# matrix of the estimated parameters.
# 
# The MARSS package provides 3 ways to compute the Hessian: 1. The recursive algorithm
# by Harvey (1989) 2. a numerical estimate using the dfHess() function from the nmle
# package, and 3. a numerical estimate from the optim() function.
# 
# The calculation of the Hessian associated with the variance terms (Q & R) is prone
# to numerical errors. When this happens, an NA is put on the diagonal of the Hessian
# for that parameter value. No standard errors or CIs can be computed for that value.
# 
# A Hessian with many NAs is probably a sign that you have a poor model (meaning your model is not a good description of the data) or you do not have enough data given the complexity of your model. <--****************************

# better to do parametric/non-parametric bootstrapping once model is decided upon
# Maybe increase to over 100 boots, 100 is standard
# est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)


# formatting confidence intervals into dataframe
CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)

### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00","AT07","GV01","HO00","MC06","RG01","RS02","EFJ","RED","RSA","RSAW")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  dplyr::rename(Site = i) %>%
  rename(Parameter = parm) # rename site column

# plot results
(RESULTS_ALL_d <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 07/23/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(.~Site, scales = "free"))
