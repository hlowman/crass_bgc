# SB & VC MARSS models with fire x ppt interactions and legacy effects
# Script started July 19, 2022
# Heili Lowman, Alex Webster

# The following script will run a MARSS analysis at SBC LTER & NM Valles Caldera sites for the CRASS project.

#### Data Notes ####

# A few notes about the process below:
# Data was tidied in sbc_marss_model.R script and exported to file:
# marss_data_sb_vc_060622.rds (as of July 19, 2022)
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
library(here)
library(MARSS)
library(naniar) 

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))


#### Import data ####

dat1 = readRDS("data_working/marss_data_sb_vc_060622.rds")

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