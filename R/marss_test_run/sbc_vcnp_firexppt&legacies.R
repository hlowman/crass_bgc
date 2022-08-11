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
# [X] Alex BY Aug 1: find within-watersheds fire area data for VC
# [] Alex BY Aug 8 : Once all above is done, create demo models with and w/o legacy effects and w and w/o fire area. These should be ready to iterate in structure to a) different numbers of legacy effects for model comparisions, b) other solutes in SB, and c) to modified Q matrices to test hypotheses about different numbers of sate processes. 

#### READ ME ####

# The following script will prepare data and run a MARSS analysis at SBC LTER & NM Valles Caldera sites for the CRASS project. Modeling sections each contain one model or one model comparision and restart the script and import data anew each time. 

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

#
#### Compile data ####
# Load datasets
# Stream Chemistry - all sites
chem <- readRDS("data_working/SBchem_edited_120321.rds")
chem_nm <- readRDS("data_working/VCNPchem_edited_120521.rds")
# Precipitation - all sites
precip <- readRDS("data_working/SBprecip_edited_120121.rds")
precip_nm <- readRDS("data_working/VCNPprecip_m_cum_edited_20220721.rds")
# Fire Events - all sites
fire <- readRDS("data_working/SBfire_edited_072922.rds")
fire_nm <- readRDS("data_working/VCNPfire_edited_072922.rds")
# Site Location information
location <- read_csv("data_raw/sbc_sites_stream_hydro.csv")
location_nm <- read_csv("data_raw/VCNP_sonde_site_codes_names.csv")

# I first need to identify the matching stream sites for the precip data
precip_ed <- precip %>%
  mutate(sitecode_match = factor(case_when(
    sitecode == "GV202" ~ "GV01",
    sitecode == "BARA" ~ "HO00", # HO201 doesn't start until 2002
    sitecode == "RG202" ~ "RG01",
    sitecode == "SMPA" ~ "SP02",
    sitecode == "GORY" ~ "AT07",
    sitecode == "CAWTP" ~ "AB00",
    sitecode == "STFS" ~ "MC06", # BOGA doesn't start until 2005
    sitecode == "ELDE" ~ "RS02"))) %>%
  dplyr::rename(cumulative_precip_mm = c_precip_mm,
                sitecode_precip = sitecode) %>%
  mutate(Day = 1) %>% # new column of "days"
  mutate(Date = make_date(Year, Month, Day))

# check edited fire data
sitez = c("AB00", "GV01", "MC06", "RG01", "RS02", "HO00",
          "EFJ", "RED", "RSA", "RSAW")
ggplot(fire %>% filter(site %in% sitez), 
                   aes(x = date, y = fire_pa)) +
    geom_point() +
    labs(x = "Date") +
    facet_wrap(.~site, scales = "free",
               ncol = 1) +
    theme_bw()
ggplot(fire %>% filter(site %in% sitez), 
       aes(x = date, y = fire_perc_ws)) +
  geom_point() +
  labs(x = "Date") +
  facet_wrap(.~site, scales = "free",
             ncol = 1) +
  theme_bw()
ggplot(fire_nm %>% filter(site %in% sitez), 
       aes(x = date, y = fire_pa)) +
  geom_point() +
  labs(x = "Date") +
  facet_wrap(.~site, scales = "free",
             ncol = 1) +
  theme_bw()
ggplot(fire_nm %>% filter(site %in% sitez), 
       aes(x = date, y = fire_perc_ws)) +
  geom_point() +
  labs(x = "Date") +
  facet_wrap(.~site, scales = "free",
             ncol = 1) +
  theme_bw()

## Timeframe Selection

# Examine precip data coverage for MARSS timescale delineation
precip_ed %>%
  drop_na(sitecode_match) %>%
  ggplot(aes(x = Year, y = sitecode_match, color = sitecode_match)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none")

# So, covariate data, which cannot be missing, can run from a maximum of 9/2002 to 7/2016.

# To avoid strange gaps in data, I'm going to start with the fire data,
# since I know it extends from 9/1/2002 to 7/1/2016. This will ensure all other datasets are
# joined to these dates in full (which was causing problems earlier).
fire_precip <- left_join(fire, precip_ed, by = c("site" = "sitecode_match", "date" = "Date"))

fire_precip <- fire_precip %>%
  mutate(year = year(date),
         month = month(date)) # and the Year/Month didn't populate, so adding in new columns

# Then, left join with chemistry so as not to lose any data.
dat <- left_join(fire_precip, chem, by = c("site", "year" = "Year", "month" = "Month"))

# Adding in dummy covariates by season
n_months <- dat$date %>%
  unique() %>%
  length()

seas_1 <- sin(2 * pi * seq(n_months) / 12)
seas_2 <- cos(2 * pi * seq(n_months) / 12)

dat <- dat %>%
  mutate(Season1 = rep(seas_1, 8),
         Season2 = rep(seas_2, 8),
         index = rep(seq(1,166), 8))

# AJW: replace NaNs with NAs
dat[is.nan(dat)] = NA

# And, inspect dataset for missing covariate data.
sum(is.na(dat$fire)) # 0
sum(is.na(dat$cumulative_precip_mm)) # 0
# Great!

#### Valles Caldera ##

# Now to do the same for NM data
precip_nm_ed <- precip_nm %>%
  dplyr::rename(sitecode_precip = ID) %>%
  mutate(Day = 1) %>% # new column of "days"
  mutate(Date = as.Date(paste(year, month, "01", sep="-")))%>%
  mutate(Year = year,
         Month = month)

## Timeframe Selection
# Examine precip data coverage for MARSS timescale delineation
precip_nm_ed %>%
  ggplot(aes(x = Year, y = sitecode_precip, color = sitecode_precip)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none")

# So, covariate data, which cannot be missing, can run from a maximum of 2005 to 2021.
# To avoid strange gaps in data, I'm going to start with the fire data,
# since I know it extends from 6/1/2005 to 3/1/2019. This will ensure all other datasets are 
# joined to these dates in full (which was causing problems earlier).
fire_precip_nm <- left_join(fire_nm, precip_nm_ed, by = c("site" = "sitecode_precip", "date" = "Date"))

# Edit chemistry dataset to permit joining
namelist <- c("San Antonio Creek - Toledo", "San Antonio Creek- Toledo", "San Antonio Creek -Toledo")

# Make abbreviations appear at appropriate sites
chem_nm_ed <- chem_nm %>%
  mutate(site = factor(case_when(
    site_code == "Redondo Creek" ~ "RED",
    site_code == "East Fork Jemez River" ~ "EFJ",
    site_code == "San Antonio - West" ~ "RSAW",
    site_code %in% namelist ~ "RSA",
    site_code == "Indios Creek" ~ "IND",
    site_code == "Indios Creek - Post Fire (Below Burn)" ~ "IND_BB",
    site_code == "Indios Creek - above burn" ~ "IND_AB",
    site_code == "Sulfur Creek" ~ "SULF",
    TRUE ~ NA_character_))) %>%
  mutate(NH4_mgL = gsub("<", "", nh4_mgL),
         NO3_mgL = gsub("<", "", nO2_nO3_mgL),
         PO4_mgL = gsub("<", "", po4_mgL)) # remove all "<" symbols

# need to also add in LODs and convert analytes appropriately
mylist <- c("Contaminated")

chem_nm_ed2 <- chem_nm_ed %>%
  replace_with_na_at(.vars = c("NH4_mgL", "NO3_mgL", "PO4_mgL"),
                     condition = ~.x %in% mylist) %>% # Remove "contaminated" and replace with NA
  mutate(NH4_mgL_lod = as.numeric(ifelse(NH4_mgL <= 0.1, 0.05, NH4_mgL)),
         NO3_mgL_lod = as.numeric(ifelse(NO3_mgL <= 0.1, 0.05, NO3_mgL)),
         PO4_mgL_lod = as.numeric(ifelse(PO4_mgL <= 0.1, 0.05, PO4_mgL))) # Report low values at 1/2 LOD

# convert to uM and summarize by month
chem_nm_ed3 <- chem_nm_ed2 %>%
  mutate(nh4_uM = (NH4_mgL_lod/80.043)*1000, # convert to uM
         no3_uM = (NO3_mgL_lod/62.0049)*1000,
         po4_uM = (PO4_mgL_lod/94.9714)*1000) 

chem_nm_monthly <- chem_nm_ed3 %>%
  dplyr::group_by(site, Year, Month) %>%
  dplyr::summarize(mean_nh4_uM = mean(nh4_uM, na.rm = TRUE),
                   mean_no3_uM = mean(no3_uM, na.rm = TRUE),
                   mean_po4_uM = mean(po4_uM, na.rm = TRUE),
                   mean_cond_uScm = mean(mean_cond_uScm, na.rm = TRUE)) %>%
  dplyr::ungroup()

# Then, left join with chemistry so as not to lose any data.
dat_nm <- left_join(fire_precip_nm, chem_nm_monthly, by = c("site", "Year", "Month"))

# Adding a plot to examine analyte availability
chem_nm_monthly[is.nan(chem_nm_monthly)] = NA

# plot of chem data availability
chem_nm_monthly %>%
  ggplot(aes(x = Year, y = site, color = site)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none")

dat_nm_trim <- dat_nm %>%
  filter(date >= "2005-06-01") 
# should be 166 records pre site to match "dat" above
# 166 records x 7 sites = 1162 records
# dates between SB and NM datasets don't need to match
# but lengths of time series do
table(dat_nm_trim$site)

# Adding in dummy covariates by season
n_months_nm <- dat_nm_trim$date %>%
  unique() %>%
  length()

seas_1_nm <- sin(2 * pi * seq(n_months_nm) / 12)
seas_2_nm <- cos(2 * pi * seq(n_months_nm) / 12)

dat_nm_trim <- dat_nm_trim %>%
  mutate(Season1 = rep(seas_1_nm, 7),
         Season2 = rep(seas_2_nm, 7),
         index = rep(seq(1,166), 7))

# AJW: replace NaNs with NAs
dat_nm_trim[is.nan(dat_nm_trim)] = NA

# Note for future me - be VERY careful with the joining above. Something weird was happening previously where precip data that IS present was simply dropping off.

# Also, because the Valles Caldera data isn't always monthly, some fire 1s were dropping off when they shouldn't, so Heili made edits 2/24/22 to rectify the situation.

#### Join Full Dataset

# And finally, join the Santa Barbara and New Mexico datasets

# Created NM dataset for joining with SB data
dat_nm_select <- dat_nm_trim %>%
  # dplyr::rename(year = Year,
  #        month = Month) %>%
  select(year, month, index, site, ws_area_m2,
         cumulative_precip_mm, 
         fire_ID, ig_date, fire_pa, ws_fire_area_m2, fire_perc_ws,
         mean_nh4_uM, mean_no3_uM, mean_po4_uM, mean_cond_uScm, 
         Season1, Season2) %>%
  mutate(region = "VC") 

dat_select <- dat %>%
  select(year, month, index, site, ws_area_m2,
         cumulative_precip_mm, 
         fire_ID, ig_date, fire_pa, ws_fire_area_m2, fire_perc_ws,
         mean_nh4_uM, mean_no3_uM, mean_po4_uM, mean_cond_uScm, 
         Season1, Season2) %>%
  mutate(region = "SB")

#dat_agu <- rbind(dat_select, dat_nm_select)
dat_new22 <- rbind(dat_select, dat_nm_select)

# And export to save progress
#saveRDS(dat_new22, "data_working/marss_data_sb_vc_022422.rds")
#saveRDS(dat_new22, "data_working/marss_data_sb_vc_060622.rds")
# with fixed NM ppt data (AJW):
#saveRDS(dat_new22, "data_working/marss_data_sb_vc_072222.rds")
# removed decay terms and fixed errors in sb fire data:
#saveRDS(dat_new22, "data_working/marss_data_sb_vc_072822.rds")
# add fire areas and make all fire data tidy/long:
saveRDS(dat_new22, "data_working/marss_data_sb_vc_072922.rds")

#### Import compiled data and add fire x ppt interactions and legacy effects ####

# dat1 = readRDS("data_working/marss_data_sb_vc_060622.rds")
# dat1 = readRDS("data_working/marss_data_sb_vc_072222.rds")
# dat1 = readRDS("data_working/marss_data_sb_vc_072822.rds")
dat1 = readRDS("data_working/marss_data_sb_vc_072922.rds")

# reorganize
names(dat1)
dat2 = dat1[,c(1,2,3,18,4,5, #"year" "month" "index" "region" "site" "ws_area_m2"
               6, # cumulative_precip_mm
               7:11, # fire info
               16,17, # seasonal effects
               12:15)] # solutes

qplot(index, fire_pa, data=dat2, colour=site, geom="path", facets = "region")

# create df of sites x fire info
firez = unique(dat1[,c(4,7,8,10,11)])
firez = firez[complete.cases(firez),]
firez$ig_date = gsub("/","-",firez$ig_date)
firez$year = year(as.Date(firez$ig_date))
firez$month = month(as.Date(firez$ig_date))
firez$date = as.Date(paste(firez$year, firez$month, "01", sep="-"))

#### Add 6 m window to fire effect ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
dat4 = dat2

# fire_pa_6m: 1=fire ignited in prior 6 m, 0=no fire ignitions in 6 m

#dat4$date = as.Date(paste(dat4$year, dat4$month, "01", sep="-"))

firedates_6m = c(firez$date,
                 firez$date + 31*1,
                 firez$date + 31*2,
                 firez$date + 31*3,
                 firez$date + 31*4,
                 firez$date + 31*5)

firedates_6m = data.frame(year = year(firedates_6m),
                          month = month(firedates_6m),
                          site = rep(firez$site, 6),
                          fire_pa_6m = 1,
                          ws_fire_area_m2_6m = rep(firez$ws_fire_area_m2, 6),
                          fire_perc_ws_6m = rep(firez$fire_perc_ws, 6))
dat5 = left_join(dat4, firedates_6m, by=c("year","month","site"))
dat5[,19:21][is.na(dat5[,19:21])] = 0

qplot(index, fire_pa_6m, data=dat5, colour=site, geom="path", facets = "region")
qplot(index, fire_perc_ws_6m, data=dat5, colour=site, geom="point", facets = "region")

#### Add fire 6 m window legacy effects ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 1 year legacy
firedates_6m_1ylegacy = data.frame(year = firedates_6m$year+1,
                                   month = firedates_6m$month,
                                   site = firedates_6m$site,
                                   fire_pa_6m_1ylegacy = 1,
                                   ws_fire_area_m2_6m_1ylegacy = rep(firez$ws_fire_area_m2, 6),
                                   fire_perc_ws_6m_1ylegacy = rep(firez$fire_perc_ws, 6))
dat6 = left_join(dat5, firedates_6m_1ylegacy, by=c("year","month","site"))
dat6[,22:24][is.na(dat6[,22:24])] = 0

# 2 year legacy
firedates_6m_2ylegacy = data.frame(year = firedates_6m$year+2,
                                   month = firedates_6m$month,
                                   site = firedates_6m$site,
                                   fire_pa_6m_2ylegacy = 1,
                                   ws_fire_area_m2_6m_2ylegacy = rep(firez$ws_fire_area_m2, 6),
                                   fire_perc_ws_6m_2ylegacy = rep(firez$fire_perc_ws, 6))
dat7 = left_join(dat6, firedates_6m_2ylegacy, by=c("year","month","site"))
dat7[,25:27][is.na(dat7[,25:27])] = 0

# 3 year legacy
firedates_6m_3ylegacy = data.frame(year = firedates_6m$year+3,
                                   month = firedates_6m$month,
                                   site = firedates_6m$site,
                                   fire_pa_6m_3ylegacy = 1,
                                   ws_fire_area_m2_6m_3ylegacy = rep(firez$ws_fire_area_m2, 6),
                                   fire_perc_ws_6m_3ylegacy = rep(firez$fire_perc_ws, 6))
dat8 = left_join(dat7, firedates_6m_3ylegacy, by=c("year","month","site"))
dat8[,28:30][is.na(dat8[,28:30])] = 0

# 4 year legacy
firedates_6m_4ylegacy = data.frame(year = firedates_6m$year+4,
                                   month = firedates_6m$month,
                                   site = firedates_6m$site,
                                   fire_pa_6m_4ylegacy = 1,
                                   ws_fire_area_m2_6m_4ylegacy = rep(firez$ws_fire_area_m2, 6),
                                   fire_perc_ws_6m_4ylegacy = rep(firez$fire_perc_ws, 6))
dat9 = left_join(dat8, firedates_6m_4ylegacy, by=c("year","month","site"))
dat9[,31:33][is.na(dat9[,31:33])] = 0

# 5 year legacy
firedates_6m_5ylegacy = data.frame(year = firedates_6m$year+5,
                                   month = firedates_6m$month,
                                   site = firedates_6m$site,
                                   fire_pa_6m_5ylegacy = 1,
                                   ws_fire_area_m2_6m_5ylegacy = rep(firez$ws_fire_area_m2, 6),
                                   fire_perc_ws_6m_5ylegacy = rep(firez$fire_perc_ws, 6))
dat10 = left_join(dat9, firedates_6m_5ylegacy, by=c("year","month","site"))
dat10[,34:36][is.na(dat10[,34:36])] = 0


#### Add fire x ppt legacy effects ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# for fire pa
dat10$fire_pa_6m_ppt = dat10$fire_pa_6m*dat10$cumulative_precip_mm
dat10$fire_pa_6m_ppt_1ylegacy = dat10$fire_pa_6m_1ylegacy*dat10$cumulative_precip_mm
dat10$fire_pa_6m_ppt_2ylegacy = dat10$fire_pa_6m_2ylegacy*dat10$cumulative_precip_mm
dat10$fire_pa_6m_ppt_3ylegacy = dat10$fire_pa_6m_3ylegacy*dat10$cumulative_precip_mm
dat10$fire_pa_6m_ppt_4ylegacy = dat10$fire_pa_6m_4ylegacy*dat10$cumulative_precip_mm
dat10$fire_pa_6m_ppt_5ylegacy = dat10$fire_pa_6m_5ylegacy*dat10$cumulative_precip_mm

# for fire area
dat10$ws_fire_area_m2_6m_ppt = dat10$ws_fire_area_m2_6m*dat10$cumulative_precip_mm
dat10$ws_fire_area_m2_6m_ppt_1ylegacy = dat10$ws_fire_area_m2_6m_1ylegacy*dat10$cumulative_precip_mm
dat10$ws_fire_area_m2_6m_ppt_2ylegacy = dat10$ws_fire_area_m2_6m_2ylegacy*dat10$cumulative_precip_mm
dat10$ws_fire_area_m2_6m_ppt_3ylegacy = dat10$ws_fire_area_m2_6m_3ylegacy*dat10$cumulative_precip_mm
dat10$ws_fire_area_m2_6m_ppt_4ylegacy = dat10$ws_fire_area_m2_6m_4ylegacy*dat10$cumulative_precip_mm
dat10$ws_fire_area_m2_6m_ppt_5ylegacy = dat10$ws_fire_area_m2_6m_5ylegacy*dat10$cumulative_precip_mm

# for fire perc burn ws
dat10$fire_perc_ws_6m_ppt = dat10$fire_perc_ws_6m*dat10$cumulative_precip_mm
dat10$fire_perc_ws_6m_ppt_1ylegacy = dat10$fire_perc_ws_6m_1ylegacy*dat10$cumulative_precip_mm
dat10$fire_perc_ws_6m_ppt_2ylegacy = dat10$fire_perc_ws_6m_2ylegacy*dat10$cumulative_precip_mm
dat10$fire_perc_ws_6m_ppt_3ylegacy = dat10$fire_perc_ws_6m_3ylegacy*dat10$cumulative_precip_mm
dat10$fire_perc_ws_6m_ppt_4ylegacy = dat10$fire_perc_ws_6m_4ylegacy*dat10$cumulative_precip_mm
dat10$fire_perc_ws_6m_ppt_5ylegacy = dat10$fire_perc_ws_6m_5ylegacy*dat10$cumulative_precip_mm


#### Add fire 2 m window legacy effects ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

firedates_2m = c(firez$date,
                 firez$date + 31*1)

firedates_2m = data.frame(year = year(firedates_2m),
                          month = month(firedates_2m),
                          site = rep(firez$site, 2),
                          fire_pa_2m = 1,
                          ws_fire_area_m2_2m = rep(firez$ws_fire_area_m2, 2),
                          fire_perc_ws_2m = rep(firez$fire_perc_ws, 2))

# 1 year legacy
firedates_2m_1ylegacy = data.frame(year = firedates_2m$year+1,
                                   month = firedates_2m$month,
                                   site = firedates_2m$site,
                                   fire_pa_2m_1ylegacy = 1,
                                   ws_fire_area_m2_2m_1ylegacy = rep(firez$ws_fire_area_m2, 2),
                                   fire_perc_ws_2m_1ylegacy = rep(firez$fire_perc_ws, 2))
dat11 = left_join(dat10, firedates_2m_1ylegacy, by=c("year","month","site"))

# 2 year legacy
firedates_2m_2ylegacy = data.frame(year = firedates_2m$year+2,
                                   month = firedates_2m$month,
                                   site = firedates_2m$site,
                                   fire_pa_2m_2ylegacy = 1,
                                   ws_fire_area_m2_2m_2ylegacy = rep(firez$ws_fire_area_m2, 2),
                                   fire_perc_ws_2m_2ylegacy = rep(firez$fire_perc_ws, 2))
dat11 = left_join(dat11, firedates_2m_2ylegacy, by=c("year","month","site"))

# 3 year legacy
firedates_2m_3ylegacy = data.frame(year = firedates_2m$year+3,
                                   month = firedates_2m$month,
                                   site = firedates_2m$site,
                                   fire_pa_2m_3ylegacy = 1,
                                   ws_fire_area_m2_2m_3ylegacy = rep(firez$ws_fire_area_m2, 2),
                                   fire_perc_ws_2m_3ylegacy = rep(firez$fire_perc_ws, 2))
dat11 = left_join(dat11, firedates_2m_3ylegacy, by=c("year","month","site"))

# 4 year legacy
firedates_2m_4ylegacy = data.frame(year = firedates_2m$year+4,
                                   month = firedates_2m$month,
                                   site = firedates_2m$site,
                                   fire_pa_2m_4ylegacy = 1,
                                   ws_fire_area_m2_2m_4ylegacy = rep(firez$ws_fire_area_m2, 2),
                                   fire_perc_ws_2m_4ylegacy = rep(firez$fire_perc_ws, 2))
dat11 = left_join(dat11, firedates_2m_4ylegacy, by=c("year","month","site"))

# 5 year legacy
firedates_2m_5ylegacy = data.frame(year = firedates_2m$year+5,
                                   month = firedates_2m$month,
                                   site = firedates_2m$site,
                                   fire_pa_2m_5ylegacy = 1,
                                   ws_fire_area_m2_2m_5ylegacy = rep(firez$ws_fire_area_m2, 2),
                                   fire_perc_ws_2m_5ylegacy = rep(firez$fire_perc_ws, 2))
dat11 = left_join(dat11, firedates_2m_5ylegacy, by=c("year","month","site"))

dat11[,55:69][is.na(dat11[,55:69])] = 0

#
#### Check for known fire x ppt effect in EFJ ####

# Sherson et al 2015 describes the effect of the 2011 Las Conchas fire on SpC in the East Fork Jemez River using high frequency data: almost no SpC response to storms before the fire, then observations of SpC flushing after. Those observations should be apparent in this dataset as well as it is the same as site VC EFJ here. This serves as a good tests that the data used here is correct (i.e. hasn't gotten messed up in editing processes), so I am checking for that observation here before proceeding. 
# For model results, this means we should also expect to see a ppt x fire effect at least in thr EFJ stream, in association with at least the 2011 fire. If we do not, it is a red flag that the model might be poorly specified, overfit, etc. and we should double check everything. It's possible the effect isn't strong enough to produce a significant coefficent, but it is something to look out for as a 'gut check'.

# pull data from EFJ 
EFJ = dat11[dat11$site=="EFJ",]
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
dat_site = dat11[dat11$site=="IND",]
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




#
#### Remove outliers ####

# See notes in previous section
# Replace outlier with previous month's value
dat11$mean_cond_uScm[dat11$site=="EFJ" & dat11$index==14] = dat11$mean_cond_uScm[dat11$site=="EFJ" & dat11$index==13] 

# winter-time unrealistic ppt values in RSAW:
# Replace outlier with previous month's value
dat11$cumulative_precip_mm[dat11$site=="RSAW" & dat11$year==2016 & dat11$month==12] = dat11$cumulative_precip_mm[dat11$site=="RSAW" & dat11$year==2016 & dat11$month==11]
dat11$cumulative_precip_mm[dat11$site=="RSAW" & dat11$year==2017 & dat11$month==1] = dat11$cumulative_precip_mm[dat11$site=="RSAW" & dat11$year==2017 & dat11$month==2]


#### Export data with fire x ppt interactions and legacy effects ####

#saveRDS(dat10, "data_working/marss_data_sb_vc_072222_2.rds")
#saveRDS(dat10, "data_working/marss_data_sb_vc_072822_2.rds")
saveRDS(dat11, "data_working/marss_data_sb_vc_072922_2.rds")

#
#### Plot data ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_072922_2.rds")

# select sites
# include these sites only (9 total - these have the longest most complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, MC06, RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "MC06", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]

# add date
dat$date = as.Date(paste(dat$year, dat$month, "01", sep="-"))

# SpC in all sites
ggplot(dat, aes(x = date, y = mean_cond_uScm)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="AT07" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RG01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="SP02" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="EFJ" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="IND" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="IND_BB" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RED" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RSA" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RSAW" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="SULF" & fire_pa==1), aes(xintercept=date), 
             colour="red")

# Fire area in all sites
ggplot(dat, aes(x = date, y = ws_fire_area_m2)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="AT07" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RG01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="SP02" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="EFJ" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="IND" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="IND_BB" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RED" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RSA" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RSAW" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="SULF" & fire_pa==1), aes(xintercept=date), 
             colour="red")

# Fire % of watershed burned in all sites
ggplot(dat, aes(x = date, y = fire_perc_ws)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="AT07" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RG01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="SP02" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="EFJ" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="IND" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="IND_BB" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RED" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RSA" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RSAW" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="SULF" & fire_pa==1), aes(xintercept=date), 
             colour="red")

# ppt in all sites
ggplot(dat, aes(x = date, y = cumulative_precip_mm)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="AT07" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RG01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="SP02" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="EFJ" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="IND" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="IND_BB" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RED" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RSA" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RSAW" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="SULF" & fire_pa==1), aes(xintercept=date), 
             colour="red")

# Fire p/a X ppt in all sites
ggplot(dat, aes(x = date, y = fire_pa_6m_ppt)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="AT07" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RG01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="SP02" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="EFJ" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="IND" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="IND_BB" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RED" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RSA" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RSAW" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="SULF" & fire_pa==1), aes(xintercept=date), 
             colour="red")

# Fire area X ppt in all sites
ggplot(dat, aes(x = date, y = ws_fire_area_m2_6m_ppt)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="AT07" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RG01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="SP02" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="EFJ" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="IND" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="IND_BB" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RED" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RSA" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RSAW" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="SULF" & fire_pa==1), aes(xintercept=date), 
             colour="red")

# Fire % watershed burned X ppt in all sites
ggplot(dat, aes(x = date, y = fire_perc_ws_6m_ppt)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="AT07" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RG01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="SP02" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="EFJ" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="IND" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="IND_BB" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RED" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RSA" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RSAW" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="SULF" & fire_pa==1), aes(xintercept=date), 
             colour="red")


#### MARSS: ppt, fire pa 2m, fire pa 6m x ppt, no legacy effects ####

#### Set up data for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_072922_2.rds")

# select sites
# include these sites only (9 total - these have the longest most complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, MC06, RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "MC06", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)
table(dat$site,dat$fire_pa)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_pa, fire_pa_6m_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_pa, fire_pa_6m_ppt)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:10)
cov_cols = c(11:37)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols]))
sum(is.na(dat_cond_log[,resp_cols]))
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_cond_log[,cov_cols]))
sum(is.na(dat_cond_log[,cov_cols]))
sum(is.infinite(dat_cond_log[,cov_cols]))

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0)
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov))
sum(is.na(dat_cov))
sum(is.infinite(dat_cov))
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
#yes

#### make C matrix  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# "XXXX_AB00",0,0,0,0,0,0,0,0,
# 0,"XXXX_GV01",0,0,0,0,0,0,0,
# 0,0,"XXXX_HO00",0,0,0,0,0,0,
# 0,0,0,"XXXX_MC06",0,0,0,0,0,
# 0,0,0,0,"XXXX_RS02",0,0,0,0,
# 0,0,0,0,0,"XXXX_EFJ" ,0,0,0,
# 0,0,0,0,0,0,"XXXX_RSAW",0,0,
# 0,0,0,0,0,0,0,"XXXX_RSA" ,0,
# 0,0,0,0,0,0,0,0,"XXXX_RED" ,

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,0,0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,0,0,0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,0,0,0,0,
  0,0,0,0,"cumulative_precip_mm_RS02",0,0,0,0,
  0,0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,0,
  0,0,0,0,0,0,"cumulative_precip_mm_RSAW",0,0,
  0,0,0,0,0,0,0,"cumulative_precip_mm_RSA" ,0,
  0,0,0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_pa: fire effect in 2 m window
  "fire_pa_AB00",0,0,0,0,0,0,0,0,
  0,"fire_pa_GV01",0,0,0,0,0,0,0,
  0,0,"fire_pa_HO00",0,0,0,0,0,0,
  0,0,0,"fire_pa_MC06",0,0,0,0,0,
  0,0,0,0,"fire_pa_RS02",0,0,0,0,
  0,0,0,0,0,"fire_pa_EFJ" ,0,0,0,
  0,0,0,0,0,0,"fire_pa_RSAW",0,0,
  0,0,0,0,0,0,0,"fire_pa_RSA" ,0,
  0,0,0,0,0,0,0,0,"fire_pa_RED" ,
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_pa_6m_ppt_AB00",0,0,0,0,0,0,0,0,
  0,"fire_pa_6m_ppt_GV01",0,0,0,0,0,0,0,
  0,0,"fire_pa_6m_ppt_HO00",0,0,0,0,0,0,
  0,0,0,"fire_pa_6m_ppt_MC06",0,0,0,0,0,
  0,0,0,0,"fire_pa_6m_ppt_RS02",0,0,0,0,
  0,0,0,0,0,"fire_pa_6m_ppt_EFJ" ,0,0,0,
  0,0,0,0,0,0,"fire_pa_6m_ppt_RSAW",0,0,
  0,0,0,0,0,0,0,"fire_pa_6m_ppt_RSA" ,0,
  0,0,0,0,0,0,0,0,"fire_pa_6m_ppt_RED" ), 9,27)

#### Model setup for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# x0_fixed = c(dat_dep[,1])
# x0_fixed[4] = mean(dat_dep[4,],na.rm=T)

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

#### Fit MARSS model +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)
# Convergence warnings
# 102 warnings. First 10 shown.  Type cat(object$errors) to see the full list.
# Warning: the  logLik  parameter value has not converged.
# Type MARSSinfo("convergence") for more info on this warning.
# 
# MARSSkem warnings. Type MARSSinfo() for help.
# iter=2,t=1 B update is outside the unit circle......

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_08102022_9state_cond_fire2mpa_fire6mpaxPPT_nolegacies_mBFGS.rds")

# the EM algorithm did not converge but still provided something to pass on to the BFGS fit so I was able to get a model fit. Unclear to me at this point whoe the BFGS model did but I should be warry. see MARSSinfo('denominv') for possibilities as to why the EM algorithm struggled. One possibility is that there are several fire and firexppt covar rows that hold the same values. Ideally that would be coded as shared covar data among states in the C matrix, but that's a pain. However I didn't have this issue until now, so could be something else. 

#### DIAGNOSES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_08102022_9state_cond_fire2mpa_fire6mpaxPPT_nolegacies_mBFGS.rds")

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
# fit        0.0 54
# null.fit 281.1 27
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

### Do resids have temporal autocorrelation? ###
par(mfrow=c(5,2),oma = c(0, 0, 1, 0))
for(i in c(1:9)){
  #forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1)
}
# These look ok.
# What you don't want is a consistent lag at 1, 6, or 12.
# Patterns are bad (esp. sinusoidal), random is good.
# Should definitely examine these without seasonal effect to see how necessary this is.

### Are resids normal? ###
par(mfrow=c(5,2),oma = c(0, 0, 1.5, 0))
for(i in c(1:9)){
  # qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"),
  #        pch=16,
  #        xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  # qqline(resids$model.residuals[i,])
  qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16, 
          xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
   qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1)
}
# state residuals - looking pretty good
# they are qq plots that should look like a straight line
# shapiro test scores should be closer to 1
# flat lines likely due to low variation in some sites

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

# None of these diagnoses look prohibitively bad

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
my_list <- c("AB00","GV01","HO00","MC06","RS02","EFJ","RED","RSA","SAW")

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
CIs_fit_ed$Parameter = rep(c("Cum. Ppt", "Fire p/a","Cum. Ppt * Fire p/a (6m)"), 9)
CIs_fit_ed$Region = c(rep(c("SB"),5*3), rep(c("VC"),4*3)) # *****CHECK ORDER OF SITES FOR THIS!!*****

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
#### MARSS: ppt, % burn 2m, % burn 6m x ppt, no legacy effects ####

#### Set up data for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_072922_2.rds")

# select sites
# include these sites only (9 total - these have the longest most complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, MC06, RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "MC06", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws, fire_perc_ws_6m_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws, fire_perc_ws_6m_ppt)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:10)
cov_cols = c(11:37)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols]))
sum(is.na(dat_cond_log[,resp_cols]))
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_cond_log[,cov_cols]))
sum(is.na(dat_cond_log[,cov_cols]))
sum(is.infinite(dat_cond_log[,cov_cols]))

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0)
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov))
sum(is.na(dat_cov))
sum(is.infinite(dat_cov))
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
#no

#### make C matrix  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# "XXXX_AB00",0,0,0,0,0,0,0,0,
# 0,"XXXX_GV01",0,0,0,0,0,0,0,
# 0,0,"XXXX_HO00",0,0,0,0,0,0,
# 0,0,0,"XXXX_MC06",0,0,0,0,0,
# 0,0,0,0,"XXXX_RS02",0,0,0,0,
# 0,0,0,0,0,"XXXX_EFJ" ,0,0,0,
# 0,0,0,0,0,0,"XXXX_RSAW",0,0,
# 0,0,0,0,0,0,0,"XXXX_RSA" ,0,
# 0,0,0,0,0,0,0,0,"XXXX_RED" ,

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,0,0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,0,0,0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,0,0,0,0,
  0,0,0,0,"cumulative_precip_mm_RS02",0,0,0,0,
  0,0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,0,
  0,0,0,0,0,0,"cumulative_precip_mm_RSAW",0,0,
  0,0,0,0,0,0,0,"cumulative_precip_mm_RSA" ,0,
  0,0,0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_perc_ws: fire effect in 2 m window
  "fire_perc_ws_AB00",0,0,0,0,0,0,0,0,
  0,"fire_perc_ws_GV01",0,0,0,0,0,0,0,
  0,0,"fire_perc_ws_HO00",0,0,0,0,0,0,
  0,0,0,"fire_perc_ws_MC06",0,0,0,0,0,
  0,0,0,0,"fire_perc_ws_RS02",0,0,0,0,
  0,0,0,0,0,"fire_perc_ws_EFJ" ,0,0,0,
  0,0,0,0,0,0,"fire_perc_ws_RSAW",0,0,
  0,0,0,0,0,0,0,"fire_perc_ws_RSA" ,0,
  0,0,0,0,0,0,0,0,"fire_perc_ws_RED" ,
  # fire_perc_ws_6m_ppt: interaction of cum. ppt with fire % burn in 6 m window
  "fire_perc_ws_6m_ppt_AB00",0,0,0,0,0,0,0,0,
  0,"fire_perc_ws_6m_ppt_GV01",0,0,0,0,0,0,0,
  0,0,"fire_perc_ws_6m_ppt_HO00",0,0,0,0,0,0,
  0,0,0,"fire_perc_ws_6m_ppt_MC06",0,0,0,0,0,
  0,0,0,0,"fire_perc_ws_6m_ppt_RS02",0,0,0,0,
  0,0,0,0,0,"fire_perc_ws_6m_ppt_EFJ" ,0,0,0,
  0,0,0,0,0,0,"fire_perc_ws_6m_ppt_RSAW",0,0,
  0,0,0,0,0,0,0,"fire_perc_ws_6m_ppt_RSA" ,0,
  0,0,0,0,0,0,0,0,"fire_perc_ws_6m_ppt_RED" ), 9,27)

#### Model setup for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# x0_fixed = c(dat_dep[,1])
# x0_fixed[4] = mean(dat_dep[4,],na.rm=T)

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

#### Fit MARSS model +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)
# Convergence warnings
# 102 warnings. First 10 shown.  Type cat(object$errors) to see the full list.
# Warning: the  logLik  parameter value has not converged.
# Type MARSSinfo("convergence") for more info on this warning.
# 
# MARSSkem warnings. Type MARSSinfo() for help.
# iter=2,t=1 B update is outside the unit circle......

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_08102022_9state_cond_percburn2m_percburn6mxppt_nolegacies_mBFGS.rds")

# the EM algorithm had warnings but still provided something to pass on to the BFGS fit so I was able to get a model fit. 

#### DIAGNOSES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# NOTE: If you start here, make sure you run the parts of the script above that prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_08102022_9state_cond_percburn2m_percburn6mxppt_nolegacies_mBFGS.rds")

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
# fit        0.0 54
# null.fit 271.5 27
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

### Do resids have temporal autocorrelation? ###
par(mfrow=c(5,2),oma = c(0, 0, 1, 0))
for(i in c(1:9)){
  #forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1)
}
# These look ok.
# What you don't want is a consistent lag at 1, 6, or 12.
# Patterns are bad (esp. sinusoidal), random is good.
# Should definitely examine these without seasonal effect to see how necessary this is.

### Are resids normal? ###
par(mfrow=c(5,2),oma = c(0, 0, 1.5, 0))
for(i in c(1:9)){
  # qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"),
  #        pch=16,
  #        xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  # qqline(resids$model.residuals[i,])
  qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16, 
         xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1)
}
# state residuals - looking pretty good
# they are qq plots that should look like a straight line
# shapiro test scores should be closer to 1
# flat lines likely due to low variation in some sites

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

# None of these diagnoses look prohibitively bad

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
my_list <- c("AB00","GV01","HO00","MC06","RS02","EFJ","RED","RSA","SAW")

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
CIs_fit_ed$Parameter = rep(c("Cum. Ppt", "% Burn","Cum. Ppt * % Burn (6m)"), 9)
CIs_fit_ed$Region = c(rep(c("SB"),5*3), rep(c("VC"),4*3)) # *****CHECK ORDER OF SITES FOR THIS!!*****

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
         title = "Sp. Conductivity MARSS modeling results - 08/10/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))


#

#### MARSS: ppt, % burn 2m, % burn 6m x ppt, 1y legacy effects ####

#### Set up data for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_072922_2.rds")

# select sites
# include these sites only (9 total - these have the longest most complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, MC06, RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "MC06", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
dat = dat[dat$site %in% sitez,]
table(dat$site)

# pivot wider for MARSS format
dat_cond <- dat %>%
  select(
    site, index, 
    mean_cond_uScm, 
    cumulative_precip_mm, 
    fire_perc_ws, fire_perc_ws_6m_1ylegacy, 
    fire_perc_ws_6m_ppt, fire_perc_ws_6m_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws, fire_perc_ws_6m_1ylegacy,
                    fire_perc_ws_6m_ppt, fire_perc_ws_6m_ppt_1ylegacy)) 

# indicate column #s of response and predictor vars
names(dat_cond)
resp_cols = c(2:10)
cov_cols = c(11:55)

# log and scale transform response var
dat_cond_log = dat_cond
dat_cond_log[,resp_cols] = log10(dat_cond_log[,resp_cols])
dat_cond_log[,resp_cols] = scale(dat_cond_log[,resp_cols])
# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_cond_log[,resp_cols]))
sum(is.na(dat_cond_log[,resp_cols]))
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_cond_log[,cov_cols]))
sum(is.na(dat_cond_log[,cov_cols]))
sum(is.infinite(dat_cond_log[,cov_cols]))

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0)
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov))

# fire_perc_ws_6m_ppt_1ylegacy_HO00 is all NaNs because there was no precip. in the months where the 1y fire legacy falls. I will replace with all zeros and ignore results from this effect.
dat_cov["fire_perc_ws_6m_ppt_1ylegacy_HO00",] = 0

sum(is.na(dat_cov))
sum(is.infinite(dat_cov))
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
#no

#### make C matrix  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# "XXXX_AB00",0,0,0,0,0,0,0,0,
# 0,"XXXX_GV01",0,0,0,0,0,0,0,
# 0,0,"XXXX_HO00",0,0,0,0,0,0,
# 0,0,0,"XXXX_MC06",0,0,0,0,0,
# 0,0,0,0,"XXXX_RS02",0,0,0,0,
# 0,0,0,0,0,"XXXX_EFJ" ,0,0,0,
# 0,0,0,0,0,0,"XXXX_RSAW",0,0,
# 0,0,0,0,0,0,0,"XXXX_RSA" ,0,
# 0,0,0,0,0,0,0,0,"XXXX_RED" ,

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,0,0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,0,0,0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,0,0,0,0,
  0,0,0,0,"cumulative_precip_mm_RS02",0,0,0,0,
  0,0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,0,
  0,0,0,0,0,0,"cumulative_precip_mm_RSAW",0,0,
  0,0,0,0,0,0,0,"cumulative_precip_mm_RSA" ,0,
  0,0,0,0,0,0,0,0,"cumulative_precip_mm_RED" ,
  # fire_perc_ws: fire effect in 2 m window
  "fire_perc_ws_AB00",0,0,0,0,0,0,0,0,
  0,"fire_perc_ws_GV01",0,0,0,0,0,0,0,
  0,0,"fire_perc_ws_HO00",0,0,0,0,0,0,
  0,0,0,"fire_perc_ws_MC06",0,0,0,0,0,
  0,0,0,0,"fire_perc_ws_RS02",0,0,0,0,
  0,0,0,0,0,"fire_perc_ws_EFJ" ,0,0,0,
  0,0,0,0,0,0,"fire_perc_ws_RSAW",0,0,
  0,0,0,0,0,0,0,"fire_perc_ws_RSA" ,0,
  0,0,0,0,0,0,0,0,"fire_perc_ws_RED" ,
  # fire_perc_ws_6m_1ylegacy
  "fire_perc_ws_6m_1ylegacy_AB00",0,0,0,0,0,0,0,0,
  0,"fire_perc_ws_6m_1ylegacy_GV01",0,0,0,0,0,0,0,
  0,0,"fire_perc_ws_6m_1ylegacy_HO00",0,0,0,0,0,0,
  0,0,0,"fire_perc_ws_6m_1ylegacy_MC06",0,0,0,0,0,
  0,0,0,0,"fire_perc_ws_6m_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,0,"fire_perc_ws_6m_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,0,"fire_perc_ws_6m_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,0,"fire_perc_ws_6m_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,0,"fire_perc_ws_6m_1ylegacy_RED" ,
  # fire_perc_ws_6m_ppt: interaction of cum. ppt with fire % burn in 6 m window
  "fire_perc_ws_6m_ppt_AB00",0,0,0,0,0,0,0,0,
  0,"fire_perc_ws_6m_ppt_GV01",0,0,0,0,0,0,0,
  0,0,"fire_perc_ws_6m_ppt_HO00",0,0,0,0,0,0,
  0,0,0,"fire_perc_ws_6m_ppt_MC06",0,0,0,0,0,
  0,0,0,0,"fire_perc_ws_6m_ppt_RS02",0,0,0,0,
  0,0,0,0,0,"fire_perc_ws_6m_ppt_EFJ" ,0,0,0,
  0,0,0,0,0,0,"fire_perc_ws_6m_ppt_RSAW",0,0,
  0,0,0,0,0,0,0,"fire_perc_ws_6m_ppt_RSA" ,0,
  0,0,0,0,0,0,0,0,"fire_perc_ws_6m_ppt_RED" ,
  # fire_perc_ws_6m_ppt_1ylegacy
  "fire_perc_ws_6m_ppt_1ylegacy_AB00",0,0,0,0,0,0,0,0,
  0,"fire_perc_ws_6m_ppt_1ylegacy_GV01",0,0,0,0,0,0,0,
  0,0,"fire_perc_ws_6m_ppt_1ylegacy_HO00",0,0,0,0,0,0,
  0,0,0,"fire_perc_ws_6m_ppt_1ylegacy_MC06",0,0,0,0,0,
  0,0,0,0,"fire_perc_ws_6m_ppt_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,0,"fire_perc_ws_6m_ppt_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,0,"fire_perc_ws_6m_ppt_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,0,"fire_perc_ws_6m_ppt_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,0,"fire_perc_ws_6m_ppt_1ylegacy_RED" ), 9,45)

#### Model setup for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# x0_fixed = c(dat_dep[,1])
# x0_fixed[4] = mean(dat_dep[4,],na.rm=T)

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

#### Fit MARSS model +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)
# Stopped at iter=1 in MARSSkem at U update. denom is not invertible.
This means some of the U (+ C) terms cannot be estimated.
Type MARSSinfo('denominv') for more info. 
par, kf, states, iter, loglike are the last values before the error.

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_08102022_9state_cond_percburn2m_percburn6mxppt_1ylegacies_mBFGS.rds")

# the EM algorithm did not converge but still provided something to pass on to the BFGS fit so I was able to get a model fit. 

#### DIAGNOSES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# NOTE: If you start here, make sure you run the parts of the script above that prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_08102022_9state_cond_percburn2m_percburn6mxppt_1ylegacies_mBFGS.rds")

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
# fit        0.0 54
# null.fit 291.9 27
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

### Do resids have temporal autocorrelation? ###
par(mfrow=c(5,2),oma = c(0, 0, 1, 0))
for(i in c(1:9)){
  #forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1)
}
# These look ok.
# What you don't want is a consistent lag at 1, 6, or 12.
# Patterns are bad (esp. sinusoidal), random is good.
# Should definitely examine these without seasonal effect to see how necessary this is.

### Are resids normal? ###
par(mfrow=c(5,2),oma = c(0, 0, 1.5, 0))
for(i in c(1:9)){
  # qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"),
  #        pch=16,
  #        xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  # qqline(resids$model.residuals[i,])
  qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16, 
         xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1)
}
# state residuals - looking pretty good
# they are qq plots that should look like a straight line
# shapiro test scores should be closer to 1
# flat lines likely due to low variation in some sites

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

# None of these diagnoses look prohibitively bad

#  PLOT RESULTS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# only do this if diagnoses look acceptable! 

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results - should bootstrap for final
est_fit <- MARSSparamCIs(fit)
# warnings!!! Means model is bad!

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
my_list <- c("AB00","GV01","HO00","MC06","RS02","EFJ","RED","RSA","SAW")

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
                             "% Burn", "% Burn 1y",
                             "Cum. Ppt * % Burn (6m)","Cum. Ppt * % Burn (6m) 1y"), 9)
CIs_fit_ed$Region = c(rep(c("SB"),5*5), rep(c("VC"),4*5)) # *****CHECK ORDER OF SITES FOR THIS!!*****

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
         title = "Sp. Conductivity MARSS modeling results - 08/10/2022") +
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