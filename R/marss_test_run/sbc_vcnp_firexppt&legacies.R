# SB & VC MARSS models with fire x ppt interactions and legacy effects
# Script started July 19, 2022
# Heili Lowman, Alex Webster

#### READ ME ####

# The following script will prepare data and run MARSS analyses at SBC LTER & NM Valles Caldera sites for the CRASS project. 
# MARSS modeling sections each contain 1 model or 1 model comparision and restart the script and import data anew each time. 
# Be sure to load libraries before starting modeling sections. 

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

# TO DO (as of July 26, 2022):
# [X] Heili BY July 29: check that a 2 month fire effect is defensible for SB data
# yes. See email btw heili & alex on July 27, 2022
# [X] Heili BY July 29: check that a 8 month window for fire to interact with ppt is defensible for SB data - we want to allow for rains starting the winter after a fire season to produce an interaction effect
# 6 or 8 mo are both defensible. See email btw heili & alex on July 27, 2022
# [X] Heili BY July 29: find within-watersheds fire area data for SB
# need to get from spatial analysis. Have asked Stevan for help. 
# [X] Alex BY Aug 1: find within-watersheds fire area data for VC
# [X] Alex BY Aug 8 : Once all above is done, create demo models with and w/o legacy effects and w and w/o fire area. These should be ready to iterate in structure to a) different numbers of legacy effects for model comparisions, b) other solutes in SB, and c) to modified Q matrices to test hypotheses about different numbers of state processes. 
# Done...ish. Models with only immediate fire and fireXppt effects are converging well. Including immediate + any legacy effects does not provide acceptable model convergence. Replacing immediate effects with legacy effects works, but you have to be careful with sites included. HO00 cannot be included with 1 y legcy effects for the firexppt interaction because there is no rain where the 1y legacy fire effect falls and this produces a all-zeros covariate row that MARSS does not like. RED cannot be included in 2y or greater legacy effects because there are too few data points after the fire, much less when you push the fire effect into the future.
# [X] Alex: replace legacy effects that ask if a fire effect *developed* years after the fire with legacy effects that ask if a fire effect *persisted* years after the fire. I think we'll get fewer spurious results with this approach. 
# [X] Alex: summary figure of results


#### Load packages - ***do this first for all sections!*** ####
library(tidyverse)
library(lubridate)
library(MARSS)
library(naniar) 
library(beepr)
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

#### Add fire persistance legacy effects ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

firez$effect_date = firez$date

# 0.5 year legacy (this is to allow a window for the fire x ppt interaction only)
firedates_0.5ylegacy = rbind(firez, 
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(1)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(2)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(3)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(4)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(5))
)
firedates_0.5ylegacy = data.frame(year = year(firedates_0.5ylegacy$effect_date),
                                month = month(firedates_0.5ylegacy$effect_date),
                                site = firedates_0.5ylegacy$site,
                                fire_pa_0.5ylegacy = 1,
                                ws_fire_area_m2_0.5ylegacy = firedates_0.5ylegacy$ws_fire_area_m2,
                                fire_perc_ws_0.5ylegacy = firedates_0.5ylegacy$fire_perc_ws)
firedates_0.5ylegacy = 
  firedates_0.5ylegacy %>% 
  group_by(year, month, site, fire_pa_0.5ylegacy) %>% 
  summarise(ws_fire_area_m2_0.5ylegacy = sum(ws_fire_area_m2_0.5ylegacy),
            fire_perc_ws_0.5ylegacy = sum(fire_perc_ws_0.5ylegacy))

# 1 year legacy
firedates_1ylegacy = rbind(firez, 
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(1)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(2)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(3)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(4)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(5)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(6)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(7)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(8)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(9)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(10)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(11))
                           )
firedates_1ylegacy = data.frame(year = year(firedates_1ylegacy$effect_date),
                                month = month(firedates_1ylegacy$effect_date),
                                site = firedates_1ylegacy$site,
                                fire_pa_1ylegacy = 1,
                                ws_fire_area_m2_1ylegacy = firedates_1ylegacy$ws_fire_area_m2,
                                fire_perc_ws_1ylegacy = firedates_1ylegacy$fire_perc_ws)
firedates_1ylegacy = 
  firedates_1ylegacy %>% 
  group_by(year, month, site, fire_pa_1ylegacy) %>% 
  summarise(ws_fire_area_m2_1ylegacy = sum(ws_fire_area_m2_1ylegacy),
            fire_perc_ws_1ylegacy = sum(fire_perc_ws_1ylegacy))


# 2 year legacy
firedates_2ylegacy = rbind(firez, 
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(1)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(2)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(3)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(4)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(5)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(6)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(7)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(8)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(9)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(10)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(11)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(12)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(13)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(14)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(15)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(16)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(17)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(18)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(19)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(20)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(21)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(22)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(23))
)
firedates_2ylegacy = data.frame(year = year(firedates_2ylegacy$effect_date),
                                month = month(firedates_2ylegacy$effect_date),
                                site = firedates_2ylegacy$site,
                                fire_pa_2ylegacy = 1,
                                ws_fire_area_m2_2ylegacy = firedates_2ylegacy$ws_fire_area_m2,
                                fire_perc_ws_2ylegacy = firedates_2ylegacy$fire_perc_ws)
firedates_2ylegacy = 
  firedates_2ylegacy %>% 
  group_by(year, month, site, fire_pa_2ylegacy) %>% 
  summarise(ws_fire_area_m2_2ylegacy = sum(ws_fire_area_m2_2ylegacy),
            fire_perc_ws_2ylegacy = sum(fire_perc_ws_2ylegacy))


# 3 year legacy
firedates_3ylegacy = rbind(firez, 
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(1)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(2)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(3)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(4)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(5)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(6)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(7)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(8)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(9)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(10)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(11)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(12)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(13)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(14)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(15)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(16)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(17)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(18)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(19)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(20)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(21)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(22)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(23)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(24)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(25)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(26)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(27)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(28)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(29)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(30)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(31)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(32)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(33)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(34)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(35))
)
firedates_3ylegacy = data.frame(year = year(firedates_3ylegacy$effect_date),
                                month = month(firedates_3ylegacy$effect_date),
                                site = firedates_3ylegacy$site,
                                fire_pa_3ylegacy = 1,
                                ws_fire_area_m2_3ylegacy = firedates_3ylegacy$ws_fire_area_m2,
                                fire_perc_ws_3ylegacy = firedates_3ylegacy$fire_perc_ws)
firedates_3ylegacy = 
  firedates_3ylegacy %>% 
  group_by(year, month, site, fire_pa_3ylegacy) %>% 
  summarise(ws_fire_area_m2_3ylegacy = sum(ws_fire_area_m2_3ylegacy),
            fire_perc_ws_3ylegacy = sum(fire_perc_ws_3ylegacy))

# 4 year legacy
firedates_4ylegacy = rbind(firez, 
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(1)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(2)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(3)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(4)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(5)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(6)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(7)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(8)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(9)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(10)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(11)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(12)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(13)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(14)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(15)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(16)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(17)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(18)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(19)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(20)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(21)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(22)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(23)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(24)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(25)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(26)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(27)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(28)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(29)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(30)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(31)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(32)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(33)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(34)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(35)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(36)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(37)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(38)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(39)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(40)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(41)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(42)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(43)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(44)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(45)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(46)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(47))
)
firedates_4ylegacy = data.frame(year = year(firedates_4ylegacy$effect_date),
                                month = month(firedates_4ylegacy$effect_date),
                                site = firedates_4ylegacy$site,
                                fire_pa_4ylegacy = 1,
                                ws_fire_area_m2_4ylegacy = firedates_4ylegacy$ws_fire_area_m2,
                                fire_perc_ws_4ylegacy = firedates_4ylegacy$fire_perc_ws)
firedates_4ylegacy = 
  firedates_4ylegacy %>% 
  group_by(year, month, site, fire_pa_4ylegacy) %>% 
  summarise(ws_fire_area_m2_4ylegacy = sum(ws_fire_area_m2_4ylegacy),
            fire_perc_ws_4ylegacy = sum(fire_perc_ws_4ylegacy))

# 5 year legacy
firedates_5ylegacy = rbind(firez, 
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(1)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(2)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(3)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(4)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(5)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(6)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(7)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(8)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(9)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(10)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(11)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(12)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(13)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(14)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(15)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(16)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(17)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(18)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(19)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(20)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(21)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(22)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(23)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(24)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(25)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(26)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(27)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(28)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(29)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(30)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(31)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(32)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(33)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(34)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(35)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(36)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(37)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(38)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(39)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(40)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(41)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(42)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(43)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(44)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(45)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(46)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(47)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(48)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(49)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(50)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(51)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(52)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(53)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(54)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(55)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(56)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(57)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(58)),
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(59))
)
firedates_5ylegacy = data.frame(year = year(firedates_5ylegacy$effect_date),
                                month = month(firedates_5ylegacy$effect_date),
                                site = firedates_5ylegacy$site,
                                fire_pa_5ylegacy = 1,
                                ws_fire_area_m2_5ylegacy = firedates_5ylegacy$ws_fire_area_m2,
                                fire_perc_ws_5ylegacy = firedates_5ylegacy$fire_perc_ws)
firedates_5ylegacy = 
  firedates_5ylegacy %>% 
  group_by(year, month, site, fire_pa_5ylegacy) %>% 
  summarise(ws_fire_area_m2_5ylegacy = sum(ws_fire_area_m2_5ylegacy),
            fire_perc_ws_5ylegacy = sum(fire_perc_ws_5ylegacy))

# join legacy dates to data
dat9 = left_join(dat2, firedates_0.5ylegacy, by=c("year","month","site"))
dat10 = left_join(dat9, firedates_1ylegacy, by=c("year","month","site"))
dat11 = left_join(dat10, firedates_2ylegacy, by=c("year","month","site"))
dat12 = left_join(dat11, firedates_3ylegacy, by=c("year","month","site"))
dat13 = left_join(dat12, firedates_4ylegacy, by=c("year","month","site"))
dat14 = left_join(dat13, firedates_5ylegacy, by=c("year","month","site"))

dat14[,19:34][is.na(dat14[,19:34])] = 0

# plot legacy effects
qplot(index, fire_pa_0.5ylegacy, data=dat14, colour=site, geom="path", facets = "region")
qplot(index, ws_fire_area_m2_1ylegacy, data=dat14, colour=site, geom="path", facets = "region")
qplot(index, fire_perc_ws_1ylegacy, data=dat14, colour=site, geom="path", facets = "region")

# plot legacy effects
qplot(index, fire_pa_5ylegacy, data=dat14, colour=site, geom="path", facets = "region")
qplot(index, ws_fire_area_m2_5ylegacy, data=dat14, colour=site, geom="path", facets = "region")
qplot(index, fire_perc_ws_5ylegacy, data=dat14, colour=site, geom="path", facets = "region")


#### Add fire x ppt legacy effects ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# for fire pa
dat14$fire_pa_ppt = dat14$fire_pa_0.5ylegacy*dat14$cumulative_precip_mm
dat14$fire_pa_ppt_1ylegacy = dat14$fire_pa_1ylegacy*dat14$cumulative_precip_mm
dat14$fire_pa_ppt_2ylegacy = dat14$fire_pa_2ylegacy*dat14$cumulative_precip_mm
dat14$fire_pa_ppt_3ylegacy = dat14$fire_pa_3ylegacy*dat14$cumulative_precip_mm
dat14$fire_pa_ppt_4ylegacy = dat14$fire_pa_4ylegacy*dat14$cumulative_precip_mm
dat14$fire_pa_ppt_5ylegacy = dat14$fire_pa_5ylegacy*dat14$cumulative_precip_mm

# for fire area
dat14$ws_fire_area_m2_ppt = dat14$ws_fire_area_m2_0.5ylegacy*dat14$cumulative_precip_mm
dat14$ws_fire_area_m2_ppt_1ylegacy = dat14$ws_fire_area_m2_1ylegacy*dat14$cumulative_precip_mm
dat14$ws_fire_area_m2_ppt_2ylegacy = dat14$ws_fire_area_m2_2ylegacy*dat14$cumulative_precip_mm
dat14$ws_fire_area_m2_ppt_3ylegacy = dat14$ws_fire_area_m2_3ylegacy*dat14$cumulative_precip_mm
dat14$ws_fire_area_m2_ppt_4ylegacy = dat14$ws_fire_area_m2_4ylegacy*dat14$cumulative_precip_mm
dat14$ws_fire_area_m2_ppt_5ylegacy = dat14$ws_fire_area_m2_5ylegacy*dat14$cumulative_precip_mm

# for fire perc burn ws
dat14$fire_perc_ws_ppt = dat14$fire_perc_ws_0.5ylegacy*dat14$cumulative_precip_mm
dat14$fire_perc_ws_ppt_1ylegacy = dat14$fire_perc_ws_1ylegacy*dat14$cumulative_precip_mm
dat14$fire_perc_ws_ppt_2ylegacy = dat14$fire_perc_ws_2ylegacy*dat14$cumulative_precip_mm
dat14$fire_perc_ws_ppt_3ylegacy = dat14$fire_perc_ws_3ylegacy*dat14$cumulative_precip_mm
dat14$fire_perc_ws_ppt_4ylegacy = dat14$fire_perc_ws_4ylegacy*dat14$cumulative_precip_mm
dat14$fire_perc_ws_ppt_5ylegacy = dat14$fire_perc_ws_5ylegacy*dat14$cumulative_precip_mm

# plot interaction effects
qplot(index, fire_pa_ppt, data=dat14, colour=site, geom="path", facets = "region")
qplot(index, ws_fire_area_m2_ppt, data=dat14, colour=site, geom="path", facets = "region")
qplot(index, fire_perc_ws_ppt, data=dat14, colour=site, geom="path", facets = "region")

qplot(index, fire_pa_ppt_1ylegacy, data=dat14, colour=site, geom="path", facets = "region")
qplot(index, ws_fire_area_m2_ppt_1ylegacy, data=dat14, colour=site, geom="path", facets = "region")
qplot(index, fire_perc_ws_ppt_1ylegacy, data=dat14, colour=site, geom="path", facets = "region")


#
#### Examine data closely ####

### Check for known fire x ppt effect in EFJ ###

# Sherson et al 2015 describes the effect of the 2011 Las Conchas fire on SpC in the East Fork Jemez River using high frequency data: almost no SpC response to storms before the fire, then observations of SpC flushing after. Those observations should be apparent in this dataset as well as it is the same as site VC EFJ here. This serves as a good tests that the data used here is correct (i.e. hasn't gotten messed up in editing processes), so I am checking for that observation here before proceeding. 
# For model results, this means we should also expect to see a ppt x fire effect at least in thr EFJ stream, in association with at least the 2011 fire. If we do not, it is a red flag that the model might be poorly specified, overfit, etc. and we should double check everything. It's possible the effect isn't strong enough to produce a significant coefficent, but it is something to look out for as a 'gut check'.

# pull data from EFJ 
EFJ = dat14[dat14$site=="EFJ",]
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

# My gut tells me the 2006 high point is not real. I am going to treat it as an outlier for now, but must discuss this with group.


# pull data from other sites:
# "RSAW"  "RSA" "IND"  "IND_BB"  "RED"  "SULF"  
dat_site = dat14[dat14$site=="RSAW",]
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

# In RED (Redondo Creek), there is an unusually low value the same month of the Thompson Ridge fire. Same thing on same date in RSAW. Both values are much lower than any other SpC values in the records of any other site. They also occured in May, which is typically a dry month in this area, and precip. records confirm that this couldn't be dilution. No other snowmelt periods in either record have low values at this time either. I think these low values are malfunction associated with the fire or lack of mantinance or something. I'm going to remove them. 

#
#### Remove outliers ####

# See notes in previous section
# Replace outlier with previous month's value
dat14$mean_cond_uScm[dat14$site=="EFJ" & dat14$index==14] = dat14$mean_cond_uScm[dat14$site=="EFJ" & dat14$index==13] 

# Replace with NAs
dat14$mean_cond_uScm[dat14$site=="RED" & dat14$index==96] =NA
dat14$mean_cond_uScm[dat14$site=="RSAW" & dat14$index==96] =NA

# winter-time unrealistic ppt values in RSAW:
# Replace outlier with previous month's value
dat14$cumulative_precip_mm[dat14$site=="RSAW" & dat14$year==2016 & dat14$month==12] = dat14$cumulative_precip_mm[dat14$site=="RSAW" & dat14$year==2016 & dat14$month==11]
dat14$cumulative_precip_mm[dat14$site=="RSAW" & dat14$year==2017 & dat14$month==1] = dat14$cumulative_precip_mm[dat14$site=="RSAW" & dat14$year==2017 & dat14$month==2]


#### Export data ready for modeling ####

#saveRDS(dat10, "data_working/marss_data_sb_vc_072222_2.rds")
#saveRDS(dat10, "data_working/marss_data_sb_vc_072822_2.rds")
#saveRDS(dat11, "data_working/marss_data_sb_vc_072922_2.rds")
saveRDS(dat14, "data_working/marss_data_sb_vc_082622.rds")

#
#### Plot data ready for modeling ####

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_082622.rds")

# select sites that are currently included in modeling, or skip this to plot all sites
# these have the longest most complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, MC06, RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "MC06", "RS02",
          "EFJ", "RED", "RSA", "RSAW")
sitez = "RSA"
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
ggplot(dat, aes(x = date, y = fire_pa_ppt)) +
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
ggplot(dat, aes(x = date, y = ws_fire_area_m2_ppt)) +
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
ggplot(dat, aes(x = date, y = fire_perc_ws_ppt)) +
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

# Fire % watershed burned X ppt 2y legacy
ggplot(dat, aes(x = date, y = fire_perc_ws_2ylegacy)) +
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

#
#### MARSS: ppt, % burn (2m win), % burn x ppt, no legacy ####

#### Set up data for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_082622.rds")

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
    fire_perc_ws, fire_perc_ws_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws, fire_perc_ws_ppt)) 

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
#which(colSums(dat_cov)==0)
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
  # fire_perc_ws_ppt: interaction of cum. ppt with fire % burn in 6 m window
  "fire_perc_ws_ppt_AB00",0,0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_GV01",0,0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_HO00",0,0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_MC06",0,0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_RS02",0,0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_EFJ" ,0,0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_RSAW",0,0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_RSA" ,0,
  0,0,0,0,0,0,0,0,"fire_perc_ws_ppt_RED" ), 9,27)

#### Model setup for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#x0_fixed = c(dat_dep[,1, drop = FALSE])
#x0_fixed[4] = mean(dat_dep[4,],na.rm=T)

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
        file = "data_working/marss_test_run/fit_08262022_9state_cond_percburn2m_percburn6mxppt_nolegacies_mBFGS.rds")

#### DIAGNOSES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# NOTE: If you start here, make sure you run the parts of the script above that prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_08262022_9state_cond_percburn2m_percburn6mxppt_nolegacies_mBFGS.rds")

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

#           dAIC df
# fit        0.0 54
# null.fit 251.5 27
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable?
# yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs?
# looks good

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off").

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# looks good, possibly with the exception of HO00

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# these look ok - the violations look random, possible with the exception of RS02, which may have a sinusoidal pattern

### Overall ###
# None of these diagnoses look prohibitively bad.

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
CIs_fit_ed$Parameter = rep(c("Cum. Ppt", "% Burn","Cum. Ppt * % Burn"), 9)
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
         title = "Sp. Conductivity MARSS modeling results - 08/26/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))


#

#### MARSS: ppt, % burn 1y legacy, % burn x ppt 1y legacy ####

#### Set up data for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_082622.rds")

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
    fire_perc_ws_1ylegacy, 
    fire_perc_ws_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_cond_uScm, cumulative_precip_mm, 
                    fire_perc_ws_1ylegacy,
                    fire_perc_ws_ppt_1ylegacy)) 

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

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols]))
sum(is.na(dat_cond_log[,cov_cols]))
sum(is.infinite(dat_cond_log[,cov_cols]))

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov))
sum(is.na(dat_cov))
sum(is.infinite(dat_cov))
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# check for cols with all zeros. this can cause model convergence issues
any(colSums(dat_cov)==0)

# check RSA data
par(mfrow=c(2,1))
dat_RSA = dat_cov[c(8,17,26),]
dat_RSA = rbind(dat_RSA, dat_dep[8,])
dat_RSA = as.data.frame(t(dat_RSA))
plot(dat_RSA$fire_perc_ws_ppt_1ylegacy_RSA, type="l", col="purple", ylim=c(-2,10))
lines(dat_RSA$cumulative_precip_mm_RSA, type="l", col="blue")
lines(dat_RSA$fire_perc_ws_1ylegacy_RSA, type="l", col="red")
lines(dat_RSA$V4, type="o", col="black")
dat_SAW = dat_cov[c(7,16,25),]
dat_SAW = rbind(dat_SAW, dat_dep[7,])
dat_SAW = as.data.frame(t(dat_SAW))
plot(dat_SAW$fire_perc_ws_ppt_1ylegacy_RSAW, type="l", col="purple", ylim=c(-2,10))
lines(dat_SAW$cumulative_precip_mm_RSAW, type="l", col="blue")
lines(dat_SAW$fire_perc_ws_1ylegacy_RSAW, type="l", col="red")
lines(dat_SAW$V4, type="o", col="black")

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
  # fire_perc_ws_1ylegacy
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,0,0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_MC06",0,0,0,0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,0,"fire_perc_ws_1ylegacy_RED" ,
  # fire_perc_ws_6m_ppt_1ylegacy
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_MC06",0,0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02",0,0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_EFJ" ,0,0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSAW",0,0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RSA" ,0,
  0,0,0,0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RED" ), 9,27)

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

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par); beep(2)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_08262022_9state_cond_percburn1ylegacy_percburn6mxppt1ylegacy_mBFGS.rds")

#### DIAGNOSES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# NOTE: If you start here, make sure you run the parts of the script above that prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_08262022_9state_cond_percburn1ylegacy_percburn6mxppt1ylegacy_mBFGS.rds")

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

#            dAIC df
# fit        0.0 54
# null.fit 279.2 27
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable?
# yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs?
# looks good

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off").
# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# looks ok

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# looks ok

### Overall ###
# None of these diagnoses look prohibitively bad.

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
CIs_fit_ed$Parameter = rep(c("Cum. Ppt", 
                             "% Burn 1y",
                             "Cum. Ppt * % Burn 1y"), 9)
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
         title = "Sp. Conductivity MARSS modeling results - 08/26/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))


#
#### MARSS: ppt, % burn 2y legacy, % burn x ppt 2y legacy NO RED ####

#### Set up data for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_082622.rds")

# select sites
# include these sites only (9 total - these have the longest most complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, MC06, RS02 = SB
# EFJ, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "MC06", "RS02",
          "EFJ", "RSA", "RSAW")
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
sum(is.nan(dat_cond_log[,resp_cols]))
sum(is.na(dat_cond_log[,resp_cols]))
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols]))
sum(is.na(dat_cond_log[,cov_cols]))
sum(is.infinite(dat_cond_log[,cov_cols]))

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov))
sum(is.na(dat_cov))
sum(is.infinite(dat_cov))
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# check for cols with all zeros. this can cause model convergence issues
any(colSums(dat_cov)==0)

#### make C matrix  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# "XXXX_AB00",0,0,0,0,0,0,0,
# 0,"XXXX_GV01",0,0,0,0,0,0,
# 0,0,"XXXX_HO00",0,0,0,0,0,
# 0,0,0,"XXXX_MC06",0,0,0,0,
# 0,0,0,0,"XXXX_RS02",0,0,0,
# 0,0,0,0,0,"XXXX_EFJ" ,0,0,
# 0,0,0,0,0,0,"XXXX_RSAW",0,
# 0,0,0,0,0,0,0,"XXXX_RSA" ,

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,0,0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,0,0,0,
  0,0,0,0,"cumulative_precip_mm_RS02",0,0,0,
  0,0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,0,"cumulative_precip_mm_RSA" ,
  # fire_perc_ws_2ylegacy
  "fire_perc_ws_2ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_2ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_2ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_2ylegacy_MC06",0,0,0,0,
  0,0,0,0,"fire_perc_ws_2ylegacy_RS02",0,0,0,
  0,0,0,0,0,"fire_perc_ws_2ylegacy_EFJ" ,0,0,
  0,0,0,0,0,0,"fire_perc_ws_2ylegacy_RSAW",0,
  0,0,0,0,0,0,0,"fire_perc_ws_2ylegacy_RSA" ,
  # fire_perc_ws_ppt_2ylegacy
  "fire_perc_ws_ppt_2ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_2ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_2ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_2ylegacy_MC06",0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_2ylegacy_RS02",0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_2ylegacy_EFJ" ,0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_2ylegacy_RSAW",0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_2ylegacy_RSA"), 8,24)

#### Model setup for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#x0_fixed = c(dat_dep[,1])
#x0_fixed[4] = mean(dat_dep[4,],na.rm=T)

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

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par); beep(2)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_08262022_8state_cond_percburn2ylegacy_percburn6mxppt2ylegacy_mBFGS.rds")

#### DIAGNOSES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# NOTE: If you start here, make sure you run the parts of the script above that prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_08262022_8state_cond_percburn2ylegacy_percburn6mxppt2ylegacy_mBFGS.rds")

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

#            dAIC df
# fit        0.0 48
# null.fit 272.8 24
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable?

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs?

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off").

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.

### Overall ###
# None of these diagnoses look prohibitively bad.

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
my_list <- c("AB00","GV01","HO00","MC06","RS02","EFJ","RSA","SAW")

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
                             "% Burn 2y",
                             "Cum. Ppt * % Burn 2y"), 8)
CIs_fit_ed$Region = c(rep(c("SB"),5*3), rep(c("VC"),3*3)) # *****CHECK ORDER OF SITES FOR THIS!!*****

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
         title = "Sp. Conductivity MARSS modeling results - 08/26/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))


#
#### MARSS: ppt, % burn 3y legacy, % burn x ppt 3y legacy NO RED ####

#### Set up data for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# remove data
rm(list=ls())
# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))
# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_vc_082622.rds")

# select sites
# include these sites only (9 total - these have the longest most complete ts for SpC and have SpC data coverage before and after fires):
# AB00, GV01, HO00, MC06, RS02 = SB
# EFJ, RSA, & RSAW = VC
sitez = c("AB00", "GV01", "HO00", "MC06", "RS02",
          "EFJ", "RSA", "RSAW")
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
sum(is.nan(dat_cond_log[,resp_cols]))
sum(is.na(dat_cond_log[,resp_cols]))
range(dat_cond_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_cond_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs b4 scaling (none allowed)
sum(is.nan(dat_cond_log[,cov_cols]))
sum(is.na(dat_cond_log[,cov_cols]))
sum(is.infinite(dat_cond_log[,cov_cols]))

# Make covariate inputs
dat_cov <- dat_cond_log[,c(cov_cols)]
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs after scaling (none allowed)
sum(is.nan(dat_cov))
sum(is.na(dat_cov))
sum(is.infinite(dat_cov))
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# check for cols with all zeros. this can cause model convergence issues
any(colSums(dat_cov)==0)

#### make C matrix  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# "XXXX_AB00",0,0,0,0,0,0,0,
# 0,"XXXX_GV01",0,0,0,0,0,0,
# 0,0,"XXXX_HO00",0,0,0,0,0,
# 0,0,0,"XXXX_MC06",0,0,0,0,
# 0,0,0,0,"XXXX_RS02",0,0,0,
# 0,0,0,0,0,"XXXX_EFJ" ,0,0,
# 0,0,0,0,0,0,"XXXX_RSAW",0,
# 0,0,0,0,0,0,0,"XXXX_RSA" ,

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,0,0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,0,0,0,
  0,0,0,0,"cumulative_precip_mm_RS02",0,0,0,
  0,0,0,0,0,"cumulative_precip_mm_EFJ" ,0,0,
  0,0,0,0,0,0,"cumulative_precip_mm_RSAW",0,
  0,0,0,0,0,0,0,"cumulative_precip_mm_RSA" ,
  # fire_perc_ws_3ylegacy
  "fire_perc_ws_3ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_3ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_3ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_3ylegacy_MC06",0,0,0,0,
  0,0,0,0,"fire_perc_ws_3ylegacy_RS02",0,0,0,
  0,0,0,0,0,"fire_perc_ws_3ylegacy_EFJ" ,0,0,
  0,0,0,0,0,0,"fire_perc_ws_3ylegacy_RSAW",0,
  0,0,0,0,0,0,0,"fire_perc_ws_3ylegacy_RSA" ,
  # fire_perc_ws_ppt_3ylegacy
  "fire_perc_ws_ppt_3ylegacy_AB00",0,0,0,0,0,0,0,
  0,"fire_perc_ws_ppt_3ylegacy_GV01",0,0,0,0,0,0,
  0,0,"fire_perc_ws_ppt_3ylegacy_HO00",0,0,0,0,0,
  0,0,0,"fire_perc_ws_ppt_3ylegacy_MC06",0,0,0,0,
  0,0,0,0,"fire_perc_ws_ppt_3ylegacy_RS02",0,0,0,
  0,0,0,0,0,"fire_perc_ws_ppt_3ylegacy_EFJ" ,0,0,
  0,0,0,0,0,0,"fire_perc_ws_ppt_3ylegacy_RSAW",0,
  0,0,0,0,0,0,0,"fire_perc_ws_ppt_3ylegacy_RSA"), 8,24)

#### Model setup for MARSS +++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#x0_fixed = c(dat_dep[,1])
#x0_fixed[4] = mean(dat_dep[4,],na.rm=T)

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

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par); beep(2)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_08262022_8state_cond_percburn3ylegacy_percburn6mxppt3ylegacy_mBFGS.rds")

#### DIAGNOSES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# NOTE: If you start here, make sure you run the parts of the script above that prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_08262022_8state_cond_percburn3ylegacy_percburn6mxppt3ylegacy_mBFGS.rds")

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

#            dAIC df
# fit        0.0 48
# null.fit 342 24
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable?
# missing data forcasts seem really low variation in VC

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs?
# EFJ resids have a bit of a temporal pattern that matches the trend in the data

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off").

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Flat lines likely due to low variation in some sites
# HO00 abd RS02 are a little concerning but overall these look acceptable.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# RS02 has a bit of a sinusoidal decaying pattern but these generally look pretty random. 

### Overall ###
# None of these diagnoses look prohibitively bad.

#  PLOT RESULTS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# only do this if diagnoses look acceptable! 

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results - should bootstrap for final
est_fit <- MARSSparamCIs(fit); beep(2)

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
my_list <- c("AB00","GV01","HO00","MC06","RS02","EFJ","RSA","SAW")

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
                             "% Burn 3y",
                             "Cum. Ppt * % Burn 3y"), 8)
CIs_fit_ed$Region = c(rep(c("SB"),5*3), rep(c("VC"),3*3)) # *****CHECK ORDER OF SITES FOR THIS!!*****

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
         title = "Sp. Conductivity MARSS modeling results - 08/16/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))


#
#### Summary figure of SpC results ####

# import model fits

noleg_fit = readRDS(file = "data_working/marss_test_run/fit_08262022_9state_cond_percburn2m_percburn6mxppt_nolegacies_mBFGS.rds")
leg1y_fit = readRDS(file = "data_working/marss_test_run/fit_08262022_9state_cond_percburn1ylegacy_percburn6mxppt1ylegacy_mBFGS.rds")
leg2y_fit = readRDS(file = "data_working/marss_test_run/fit_08262022_8state_cond_percburn2ylegacy_percburn6mxppt2ylegacy_mBFGS.rds")
leg3y_fit = readRDS(file = "data_working/marss_test_run/fit_08262022_8state_cond_percburn3ylegacy_percburn6mxppt3ylegacy_mBFGS.rds")

noleg_est <- MARSSparamCIs(noleg_fit); beep(2)
leg1y_est <- MARSSparamCIs(leg1y_fit); beep(2)
leg2y_est <- MARSSparamCIs(leg2y_fit); beep(2)
leg3y_est <- MARSSparamCIs(leg3y_fit); beep(2)


# formatting confidence intervals into dataframe
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

CIs = rbind(
  noleg_CI, leg1y_CI, leg2y_CI, leg3y_CI
)

# RSA and RSAW get mixed up in plotting, so replacing RSAW with SAW
CIs$Parameter =gsub("RSAW","SAW",CIs$Parameter)
rownames(CIs) =gsub("RSAW","SAW",rownames(CIs))

# add col for site names
CIs$Stream = gsub("_","",str_sub(CIs$Parameter, start= -4))
# add col for region
CIs$Region = ifelse(CIs$Stream %in% c("EFJ","SAW","RSA","RED"),"VC","SB")
# simplify parm names
CIs$Parm_simple = c(rep("Ppt",9),
                    rep("Perc. burn",9),
                    rep("Ppt x Perc. burn",9),
                    
                    rep("Ppt",9),
                    rep("Perc. burn",9),
                    rep("Ppt x Perc. burn",9),
                    
                    rep("Ppt",8),
                    rep("Perc. burn",8),
                    rep("Ppt x Perc. burn",8),
                    
                    rep("Ppt",8),
                    rep("Perc. burn",8),
                    rep("Ppt x Perc. burn",8))

# Add column to designate those sites at which effects are significant.
CIs_ALL <- CIs %>%
  mutate(sig = factor(case_when(`Est.` > 0 & 
                                  Lower > 0 & 
                                  Upper > 0 ~ "sig_pos",
                                `Est.` < 0 & 
                                  Lower < 0 & 
                                  Upper < 0 ~ "sig_neg",
                                TRUE ~ "not_sig"), 
                      levels = c("sig_pos", "not_sig", "sig_neg"))) %>%
  # and column to denote the site/year where model did not
  # converge properly
  mutate(converged = factor(case_when(Stream == c("MC06","RSA") & Model == "0 year window" ~ FALSE,
                                      TRUE ~ TRUE))) 

# factorize levels
CIs_ALL$Region = factor(CIs_ALL$Region, levels=c("VC","SB"))
CIs_ALL$Parm_simple = factor(CIs_ALL$Parm_simple, levels=c("Ppt x Perc. burn","Perc. burn","Ppt"))
CIs_ALL$Stream = factor(CIs_ALL$Stream, levels=c("EFJ","RED","RSA","SAW",
                                                     "AB00","GV01","HO00","MC06","RS02"))

my_palette <- c("black", "white", "black")

# plot results
RESULTS_ALL_d <- ggplot(CIs_ALL %>%
                          # filter for models that converged
                          filter(converged == TRUE), 
                        aes(x = factor(Parm_simple, 
                                       levels = c("Ppt x Perc. burn",
                                                  "Perc. burn",
                                                  "Ppt")),
                            Est., fill=sig, shape = Stream)) + 
  geom_errorbar(aes(ymin=Lower, ymax=Upper),
                position=position_dodge(width=0.5), width=0) +
  geom_point(position=position_dodge(width=0.5), 
             alpha = 0.8, size=8) + 
  scale_shape_manual(values = c(21, 22, 23, 24,
                                21, 22, 23, 24, 25)) +
  scale_fill_manual(values = my_palette) +
  theme_bw()+
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 20),
        strip.text.x = element_text(size = 24),
        strip.text.y = element_text(size = 24),
        legend.title=element_text(size = 20), 
        legend.text=element_text(size = 20)) +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  coord_flip(ylim = c(-1, 1)) + 
  ylim(-1, 1) +
  labs(y = "Effect Size", 
       x = "Covariates",
       fill = "Significance") +
  theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
  guides(shape=guide_legend("Stream"), fill = "none") +
  facet_grid(Region~Model)

# Export plot.
ggsave(("Fig_MARSS_SpC_SB_VC.png"),
       path = "figures",
       width = 50,
       height = 24,
       units = "cm"
)

#### Example time series of predictors figure ####

dat = readRDS("data_working/marss_data_sb_vc_082622.rds")
sitez = c("EFJ")
dat = dat[dat$site %in% sitez,]
par(mfrow=c(3,4), mai = c(0.05, .6, 0.05, 0.05))
plot(dat$index, dat$cumulative_precip_mm, type="l", col="blue", xlab="", xaxt='n', ylab="Ppt",axes=F)
plot(dat$index, dat$cumulative_precip_mm, type="l", col="blue", xlab="", xaxt='n', ylab="",axes=F)
plot(dat$index, dat$cumulative_precip_mm, type="l", col="blue", xlab="", xaxt='n', ylab="",axes=F)
plot(dat$index, dat$cumulative_precip_mm, type="l", col="blue", xlab="", xaxt='n', ylab="",axes=F)
plot(dat$index, dat$fire_perc_ws, type="l", col="red", xlab="", xaxt='n', ylab="Perc. burn",axes=F, ylim = c(0,60))
plot(dat$index, dat$fire_perc_ws_1ylegacy, type="l", col="red", xlab="", xaxt='n', ylab="",axes=F, ylim = c(0,60))
plot(dat$index, dat$fire_perc_ws_2ylegacy, type="l", col="red", xlab="", xaxt='n', ylab="",axes=F, ylim = c(0,60))
plot(dat$index, dat$fire_perc_ws_3ylegacy, type="l", col="red", xlab="", xaxt='n', ylab="",axes=F, ylim = c(0,60))
plot(dat$index, dat$fire_perc_ws_ppt, type="l", col="purple", xaxt='n', ylab="Ppt x Perc. burn",axes=F, ylim = c(0,14000))
plot(dat$index, dat$fire_perc_ws_ppt_1ylegacy, type="l", col="purple", xaxt='n', ylab="",axes=F, ylim = c(0,14000))
plot(dat$index, dat$fire_perc_ws_ppt_2ylegacy, type="l", col="purple", xaxt='n', ylab="",axes=F, ylim = c(0,14000))
plot(dat$index, dat$fire_perc_ws_ppt_3ylegacy, type="l", col="purple", xaxt='n', ylab="",axes=F, ylim = c(0,14000))
