# SB & VC MARSS models with fire x ppt interactions and legacy effects
# Script started July 19, 2022
# Heili Lowman, Alex Webster

# CURRENT TO DO (as of July 26, 2022):
# [X] Heili BY July 29: check that a 2 month fire effect is defensible for SB data
    # yes. See email btw heili & alex on July 27, 2022
# [X] Heili BY July 29: check that a 8 month window for fire to interact with ppt is defensible for SB data - we want to allow for rains starting the winter after a fire season to produce an interaction effect
    # 6 or 8 mo are both defensible. See email btw heili & alex on July 27, 2022
# [X] Heili BY July 29: find within-watersheds fire area data for SB
    # need to get from spatial analysis. Have asked Stevan for help. 
# [] Alex BY Aug 1: find within-watersheds fire area data for VC
# [] Alex BY Aug 8 : Once all above is done, create demo models with and w/o legacy effects and w and w/o fire area. These should be ready to iterate in structure to a) different numbers of legacy effects for model comparisions, b) other solutes in SB, and c) to modified Q matrices to test hypotheses about different numbers of sate processes. 

#### READ ME ####

# The following script will prepare data and run a MARSS analysis at SBC LTER & NM Valles Caldera sites for the CRASS project. Modeling sections each contain one model or one model comparision and restart the script and import data anew each time. 

# A few notes about the process below:
# Data was tidied in sbc_marss_model.R script and exported to file:
# marss_data_sb_vc_072222.rds (as of July 22, 2022 when errors in NM ppt data were identified and fixed)
# It is imported from that file into this script and futher prepared for modeling here.

# The first goal of this script is to demo a MARSS model with each watershed as a unique state and the following predictive vars:
# 1 - cum. monthly ppt
# 2 - fire_pa_6m  (1=fire ignited in prior 6 m, 0=no fire ignitions in 6 m, one col for all fires) OR fire_pa (1-fire ignited in this month, 0=no fire ignition this month). fire_pa is an immediate effect of the fire by itself, whereas fire_pa_6m is an effect of the fire by itself that is potentially lagged by up to 6 m. 
# 3 - fire_pa_6m_1ylegacy (1=fire ignited 6 m - 1 y ago, 0=no fire ignitions 6 m - 1 y ago) OR fire_pa_1ylegacy (1=fire ignited 1 y ago)
# 4 - fire_pa_6m_2ylegacy (1=fire ignited 1.5 - 2 y ago, 0=no fire ignitions 1.5 - 2 y ago) OR fire_pa_2ylegacy (1=fire ignited 2 y ago)
# 5 - fire_pa_6m_3ylegacy (1=fire ignited 2.5 - 3 y ago, 0=no fire ignitions 2.5 - 3 y ago) OR fire_pa_3ylegacy (1=fire ignited 3 y ago)
# 6 - fire_pa_6m_4ylegacy (1=fire ignited 3.5 - 4 y ago, 0=no fire ignitions 3.5 - 4 y ago) OR fire_pa_4ylegacy (1=fire ignited 4 y ago)
# 7 - fire_pa_6m_5ylegacy (1=fire ignited 4.5 - 5 y ago, 0=no fire ignitions 4.5 - 5 y ago) OR fire_pa_5ylegacy (1=fire ignited 5 y ago)
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


#### Import and edit data to include fire x ppt interactions and legacy effects ####

# dat1 = readRDS("data_working/marss_data_sb_vc_060622.rds")
dat1 = readRDS("data_working/marss_data_sb_vc_072222.rds")


#### Consolidate fire effect to one col ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

#### Add 6 m window to fire effect ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

#### Add fire p/a legacy effects ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

#### Add fire p/a x ppt legacy effects ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

#### Check for known fire x ppt effect in EFJ ####

# Sherson et al 2015 describes the effect of the 2011 Las Conchas fire on SpC in the East Fork Jemez River using high frequency data: almost no SpC response to storms before the fire, then observations of SpC flushing after. Those observations should be apparent in this dataset as well as it is the same as site VC EFJ here. This serves as a good tests that the data used here is correct (i.e. hasn't gotten messed up in editing processes), so I am checking for that observation here before proceeding. 
# For model results, this means we should also expect to see a ppt x fire effect at least in thr EFJ stream, in association with at least the 2011 fire. If we do not, it is a red flag that the model might be poorly specified, overfit, etc. and we should double check everything. It's possible the effect isn't strong enough to produce a significant coefficent, but it is something to look out for as a 'gut check'.

# pull data from EFJ 
EFJ = dat10[dat10$site=="EFJ",]
# plot
par(mfrow=c(2,1))
plot(EFJ$cumulative_precip_mm~EFJ$index, type="b",
     xlab="Date", ylab="Cum. Precip.")
abline(v=74, col="red")
abline(v=97, col="red")
plot(EFJ$mean_cond_uScm~EFJ$index, type="b",
     xlab="Date", ylab="SpC")
abline(v=74, col="red")
abline(v=97, col="red")
# Notes:
# I wish I knew whether that high SpC value at index 14 (in 2006) was an outlier or not. It does not look like a real C-Q response because there is no falling limb despite the elevated precip occuring for more than 1 month, unlike later in the dataset. However, precip. IS high that month so it could be real. Looking in other sites, there seems to be a similar event in IND but not in RSA or RSAW; data is not available on those dates in RED, IND_BB, or SULF. 
# Other than that point, this data does seem to corroborate Sherson 2015's observations. Similar precip elsewhere in the timeseries does not have a SpC flush until after the fires. 

# pull data from other sites:
# "RSAW"  "RSA" "IND"  "IND_BB"  "RED"  "SULF"  
dat_site = dat10[dat10$site=="IND",]
# plot
par(mfrow=c(2,1))
plot(dat_site$cumulative_precip_mm~dat_site$index, type="b",
     xlab="Date", ylab="Cum. Precip.")
abline(v=14, col="blue")
abline(v=74, col="red")
abline(v=97, col="red")
plot(dat_site$mean_cond_uScm~dat_site$index, type="b",
     xlab="Date", ylab="SpC")
abline(v=14, col="blue")
abline(v=74, col="red")
abline(v=97, col="red")

# My gut tells me the 2006 high point is not real. I am going to treat it as an outlier for now, but must discuss this with group.

# Replace outlier with previous month's value
dat10$mean_cond_uScm[dat10$site=="EFJ" & dat10$index==14] = dat10$mean_cond_uScm[dat10$site=="EFJ" & dat10$index==13] 


# other data checks
#Stevan's fire mapping suggests that the Gaviota fire impacted HO00
dat_site = dat10[dat10$site=="HO00",]
# plot
par(mfrow=c(2,1))
plot(dat_site$HO00_allfires~dat_site$date, type="b",
     xlab="Date", ylab="")
abline(v=14, col="blue")
abline(v=74, col="red")
abline(v=97, col="red")
plot(dat_site$mean_cond_uScm~dat_site$index, type="b",
     xlab="Date", ylab="SpC")
abline(v=14, col="blue")
abline(v=74, col="red")
abline(v=97, col="red")


#
#### Export data with fire x ppt interactions and legacy effects ####

saveRDS(dat10, "data_working/marss_data_sb_vc_072222_2.rds")

#### MARSS: ppt, fire pa 1m, fire pa 6m x ppt, no legacy effects ####

#### Set up data for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_072222_2.rds")

# select sites
# include these sites only (10 total - these have the longest most complete ts for SpC and I have also removed HO00 because it was causing issues with missing fire effects data):
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
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_07262022_10state_cond_fire1mpa_fire6mpaxPPT_nolegacies_mBFGS.rds")

#### DIAGNOSES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_07262022_10state_cond_fire1mpa_fire6mpaxPPT_nolegacies_mBFGS.rds")

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
# don't use residuals() - this will pull an entirely different df
resids <- MARSSresiduals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

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
#fit        0.0 60
#null.fit 221.6 30
# RESULT: covar model is better than null

### Plot response vars ###
par(mfrow=c(5,2),oma = c(0, 0, 0, 0))
plot(dat_dep[1,], type="o")
plot(dat_dep[2,], type="o")
plot(dat_dep[3,], type="o")
plot(dat_dep[4,], type="o")
plot(dat_dep[5,], type="o")
plot(dat_dep[6,], type="o")
plot(dat_dep[7,], type="o")
plot(dat_dep[8,], type="o")
plot(dat_dep[9,], type="o")
plot(dat_dep[10,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(5,2),oma = c(0, 0, 1, 0))
for(i in c(1:10)){
  #forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1)
}
# These look good except for the 4th state, MC06. This state has a significant sinusidal pattern. Need to consider a seasonal effect.
# What you don't want is a consistent lag at 1, 6, or 12.
# Patterns are bad (esp. sinusoidal), random is good.
# Should definitely examine these without seasonal effect to see how necessary this is.

### Are resids normal? ###
par(mfrow=c(5,2),oma = c(0, 0, 1.5, 0))
for(i in c(1:10)){
  # qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"),
  #        pch=16,
  #        xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  # qqline(resids$model.residuals[i,])
  qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16, 
          xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
   qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1)
}
# state residuals - not looking great
# they are qq plots that should look like a straight line
# shapiro test scores should be closer to 1
# Press back arrow to see all 12 states
# flat lines likely due to low variation in some sites

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

#  PLOT RESULTS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# only do this if diagnoses look acceptable! 

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results - should bootstrap for final
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
CIs_fit_ed$Region = c(rep(c("SB"),6*3), rep(c("VC"),4*3)) # *****CHECK ORDER OF SITES FOR THIS!!*****

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
         title = "Sp. Conductivity MARSS modeling results - 07/26/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))


#
#### MARSS: ppt, fire pa 1m, fire pa 6m x ppt, 5y legacy effects - IN PROGRESS ####

#### Set up data for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_072222_2.rds")

# select sites
# include these sites only (10 total - these have the longest most complete ts for SpC and I have also removed HO00 because it was causing issues with missing fire effects data):
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
    fires_pa, fire_pa_6m_1ylegacy, fire_pa_6m_2ylegacy, fire_pa_6m_3ylegacy, fire_pa_6m_4ylegacy, fire_pa_6m_5ylegacy,
    fire_pa_6m_ppt, fire_pa_6m_ppt_1ylegacy, fire_pa_6m_ppt_2ylegacy, fire_pa_6m_ppt_3ylegacy, fire_pa_6m_ppt_4ylegacy, fire_pa_6m_ppt_5ylegacy) %>% 
  pivot_wider(
    names_from = site, values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                                       fires_pa, fire_pa_6m_1ylegacy, fire_pa_6m_2ylegacy, fire_pa_6m_3ylegacy, fire_pa_6m_4ylegacy, fire_pa_6m_5ylegacy,
                                       fire_pa_6m_ppt, fire_pa_6m_ppt_1ylegacy, fire_pa_6m_ppt_2ylegacy, fire_pa_6m_ppt_3ylegacy, fire_pa_6m_ppt_4ylegacy, fire_pa_6m_ppt_5ylegacy)) 

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

#### make C matrix  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_07262022_10state_cond_fire1mpa_fire6mpaxPPT_5ylegacies_mBFGS.rds")

#### DIAGNOSES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# IF YOU START HERE, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_07262022_10state_cond_fire1mpa_fire6mpaxPPT_5ylegacies_mBFGS.rds")

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# 

## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
# don't use residuals() - this will pull an entirely different df
resids <- MARSSresiduals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

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
# you want the null model to have a HIGHER AIC score than the covariate model
# at least 2 pts different
# if null is lower or equal, suggests that adding covariates does not help the model explain variation in the data and there is a trade off with the added complexity.
#           dAIC df
# fit        0.0 150
# null.fit 269.1 30 
# RESULT: BAD. Null should not be lower.

### Plot response vars ###
par(mfrow=c(5,2),oma = c(0, 0, 0, 0))
plot(dat_dep[1,], type="o")
plot(dat_dep[2,], type="o")
plot(dat_dep[3,], type="o")
plot(dat_dep[4,], type="o")
plot(dat_dep[5,], type="o")
plot(dat_dep[6,], type="o")
plot(dat_dep[7,], type="o")
plot(dat_dep[8,], type="o")
plot(dat_dep[9,], type="o")
plot(dat_dep[10,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(5,2),oma = c(0, 0, 1, 0))
for(i in c(1:10)){
  #forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1)
}
# What you don't want is a consistent lag at 1, 6, or 12.
# Patterns are bad (esp. sinusoidal), random is good.
# RESULTS: Looks good EXCEPT for state 4, which looks like it has a lingering seasonal pattern

### Are resids normal? ###
par(mfrow=c(5,2),oma = c(0, 0, 1.5, 0))
for(i in c(1:10)){
  # qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"),
  #        pch=16,
  #        xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  # qqline(resids$model.residuals[i,])
  qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16, 
         xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1)
}
# These are qq plots that should look like a straight line
# shapiro test scores should be close to 1
# flat lines likely due to low variation in some sites, so may be unavoidable in some cases
# RESULTS: looks decent except in low variation sites

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

#  PLOT RESULTS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# only do this if diagnoses look acceptable! 

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results - should bootstrap for final
est_fit <- MARSSparamCIs(fit)
# DID NOT CALCULATE
# A Hessian with many NAs is probably a sign that you have a poor model (meaning your model is not a good description of the data) or you do not have enough data given the complexity of your model.

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
CIs_fit_ed$Parameter = rep(c("Cum. Ppt", 
                             "Fire p/a","1y Legacy Fire p/a","2y Legacy Fire p/a","3y Legacy Fire p/a","4y Legacy Fire p/a","5y Legacy Fire p/a",
                             "Cum. Ppt * Fire p/a (6m)","1y Legacy Cum. Ppt * Fire p/a (6m)","2y Legacy Cum. Ppt * Fire p/a (6m)","3y Legacy Cum. Ppt * Fire p/a (6m)","4y Legacy Cum. Ppt * Fire p/a (6m)","5y Legacy Cum. Ppt * Fire p/a (6m)"), 
                           10)
CIs_fit_ed$Region = c(rep(c("SB"),6*3), rep(c("VC"),4*3)) # *****CHECK ORDER OF SITES FOR THIS!!*****

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
         title = "Sp. Conductivity MARSS modeling results - 07/26/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))


#