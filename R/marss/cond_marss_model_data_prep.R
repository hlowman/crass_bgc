# Conductivity MARSS models
# with fire x ppt interactions and legacy effects
# as well as 8 state/2 state structures for NM + CA sites
# Script started August 11, 2023 by Heili Lowman

# This script will prep data for running in the conductivity MARSS models.

#### Setup ####

# Load packages.
library(tidyverse)
library(lubridate)
library(MARSS)
library(naniar) 

# Load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))

# Load data.
# Stream Chemistry
chem_sb <- readRDS("data_working/SBchem_edited_120321.rds")
chem_vc <- readRDS("data_working/VCNP_monthly_conductivity_sonde_grab_081123.rds")

# Precipitation
precip_sb <- readRDS("data_working/SBprecip_edited_120121.rds")
precip_vc <- readRDS("data_working/VCNPprecip_m_cum_edited_20220721.rds")

# Fire Events
fire_sb <- readRDS("data_working/SBfire_edited_011323.rds")
fire_vc <- readRDS("data_working/VCNPfire_edited_081123.rds")

#### Compile data ####

##### Santa Barbara, CA #####

# These are the sites we are fitting models for in Santa Barbara.
sites <- c("AB00", "GV01", "HO00", "RS02")

# I first need to identify the matching stream sites for the precip data
precip_sb_ed <- precip_sb %>%
  mutate(sitecode_match = factor(case_when(
    sitecode == "GV202" ~ "GV01",
    sitecode == "HO202" ~ "HO00", # switched back from "BARA" since Q data doesn't
    # start until 2002 either, using HO202 instead of 201 because record is longer
    sitecode == "CAWTP" ~ "AB00",
    sitecode == "ELDE" ~ "RS02"))) %>%
  dplyr::rename(cumulative_precip_mm = c_precip_mm,
                sitecode_precip = sitecode) %>%
  filter(sitecode_match %in% sites) %>%
  mutate(Day = 1) %>% # new column of "days"
  mutate(Date = make_date(Year, Month, Day))

# And filter the fire dataset also.
fire_sb_ed <- fire_sb %>%
  filter(site %in% sites)

ggplot(fire_sb_ed, 
       aes(x = date, y = fire_perc_ws)) +
  geom_point() +
  labs(x = "Date") +
  facet_wrap(.~site, scales = "free",
             ncol = 1) +
  theme_bw() # looks good based on fire dates in "sbc_fire_compilation.R" script

# Examine precip data coverage for MARSS timescale delineation.
precip_sb_ed %>%
  ggplot(aes(x = Year, y = sitecode_match, color = sitecode_match)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") # 2002-2017

# Covariate data, which cannot be missing, can run from a maximum of 9/2002 to 10/2017.
# All datasets MUST be the same length, and therefore I must always go with the shortest
# common timeseries, erring on the side of post-fire data inclusion.

# To avoid strange gaps in data, I'm going to start with the fire data,
# since I know it extends from 9/1/2002 to 10/1/2017. This will ensure all other 
# data are joined to these dates in full (which was causing problems earlier).

# I've already filtered the precip & fire data down to just the 4 sites I want above.
fire_precip_sb <- left_join(fire_sb_ed, precip_sb_ed, 
                            by = c("site" = "sitecode_match", "date" = "Date"))

fire_precip_sb <- fire_precip_sb %>%
  mutate(year = year(date),
         month = month(date)) 
# the Year/Month didn't populate, 
# so adding in new columns

# Select only conductivity chem data
cond_sb <- chem_sb %>%
  select(site, Year, Month, mean_cond_uScm)

# Then, left join with chemistry so as not to lose any data.
fire_precip_chem_sb <- left_join(fire_precip_sb, cond_sb, by = c("site", "Year", "Month"))

# Adding index values in for future use.
dat_sb <- fire_precip_chem_sb %>%
  mutate(index = rep(seq(1,182), 4))

# Replace NaNs with NAs
dat_sb[is.nan(dat_sb)] <- NA

# And, inspect dataset for missing covariate data.
sum(is.na(dat_sb$fire_perc_ws)) # 0
sum(is.na(dat_sb$cumulative_precip_mm)) # 2
# Need to trim out October 2017. Ok, because it's at the end 
# rather than in the middle of time series.

dat_sb <- dat_sb %>%
  filter(date != as_date(as.character("2017-10-01")))

# Trim down to columns of interest
dat_select_sb <- dat_sb %>%
  select(year, month, index, site, ws_area_m2, # date/site info
         cumulative_precip_mm,                 # precip data
         fire_ID, ig_date, fire_pa, ws_fire_area_m2, fire_perc_ws, # fire data
         mean_cond_uScm) %>%                   # chem data
  mutate(region = "SB")

##### Valles Caldera, NM #####

# These are the sites we are fitting models for in Valles Caldera.
sites2 <- c("EFJ", "RED", "RSA", "RSAW")

# So, unlike SB, VC precip data is already identified by site.
# I am simply going to filter by the sites we want.
precip_vc_ed <- precip_vc %>%
  filter(ID %in% sites2) %>%
  mutate(day = 1) %>% # new column of "days"
  mutate(Date = make_date(year, month, day))

# And filter the fire dataset also.
fire_vc_ed <- fire_vc %>%
  filter(site %in% sites2)

ggplot(fire_vc_ed, 
       aes(x = date, y = fire_perc_ws)) +
  geom_point() +
  labs(x = "Date") +
  facet_wrap(.~site, scales = "free",
             ncol = 1) +
  theme_bw() # looks good based on fire dates in "vc_fire_compilation.R" script

# Examine precip data coverage for MARSS timescale delineation.
precip_vc_ed %>%
  ggplot(aes(x = year, y = ID, color = ID)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") # 2004-2021

# Covariate data, which cannot be missing, can run from a maximum of 8/2004 to 10/2021.

# To avoid strange gaps in data, I'm going to start with the fire data.
# This will ensure all other data are joined to these dates in full.

# I've already filtered the precip & fire data down to just the 4 sites I want above.
fire_precip_vc <- left_join(fire_vc_ed, precip_vc_ed, 
                            by = c("site" = "ID", 
                                   "date" = "Date"))

# Then, left join with chemistry so as not to lose any data.
fire_precip_chem_vc <- left_join(fire_precip_vc, chem_vc, by = c("site" = "site_name", 
                                                                 "year" = "Y",
                                                                 "month" = "M"))

# Adding index values in for future use.
dat_vc <- fire_precip_chem_vc %>%
  mutate(index = rep(seq(1, 207), 4))

# Replace NaNs with NAs
dat_vc[is.nan(dat_vc)] <- NA

# And, inspect dataset for missing covariate data.
sum(is.na(dat_vc$fire_perc_ws)) # 0
sum(is.na(dat_vc$cumulative_precip_mm)) # 0

# Trim down to columns of interest
dat_select_vc <- dat_vc %>%
  select(year, month, index, site, ws_area_m2, # date/site info
         cumulative_precip_mm,                 # precip data
         fire_ID, ig_date, fire_pa, ws_fire_area_m2, fire_perc_ws, # fire data
         mean_cond_uScm) %>%                   # chem data
  mutate(region = "VC")

# Because this dataset is longer than the SB dataset, I must trim it down.
# Each timeseries must be the exact same length, although not necessarily during the
# same time period, so I need to trim down to 181 samples.

# Using the same reasoning as above, I am going to trim off the earlier data.
dat_select_vc_181 <- dat_select_vc %>%
  filter(index > 26) %>%
  select(-index) %>%
  mutate(index = rep(seq(1, 181), 4))

# Now, both datasets are the same size (724 observations).

# Join and export to save progress.
dat_select <- rbind(dat_select_sb, dat_select_vc_181)
#saveRDS(dat_select, "data_working/marss_data_sb_vc_nolegacies_081423.rds")

#### Add fire legacy effects ####

dat1 <- readRDS("data_working/marss_data_sb_vc_nolegacies_081423.rds")

# reorganize
names(dat1)
dat2 = dat1[,c(1,2,3,13,4,5, #"year" "month" "index" "region" "site" "ws_area_m2"
               6, # cumulative_precip_mm
               7:11, # fire info
               12)] # conductivity

# create df of sites x fire info
firez = unique(dat1[,c(4,7,8,10,11)])
firez = firez[complete.cases(firez),]
firez$ig_date = gsub("/","-",firez$ig_date)
firez$year = year(as.Date(firez$ig_date))
firez$month = month(as.Date(firez$ig_date))
firez$date = as.Date(paste(firez$year, firez$month, "01", sep="-"))

# Add fire persistance legacy effects
firez$effect_date = firez$date

# 0.5 year legacy (this is to allow a window for the fire x ppt interaction only)
firedates_0.5ylegacy = rbind(firez, 
                             cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(1)),
                             cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(2)),
                             cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(3)),
                             cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(4)),
                             cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(5)))

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
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(11)))

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
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(23)))

firedates_2ylegacy = data.frame(year = year(firedates_2ylegacy$effect_date),
                                month = month(firedates_2ylegacy$effect_date),
                                site = firedates_2ylegacy$site,
                                fire_pa_2ylegacy = 1,
                                ws_fire_area_m2_2ylegacy = firedates_2ylegacy$ws_fire_area_m2,
                                fire_perc_ws_2ylegacy = firedates_2ylegacy$fire_perc_ws)

firedates_2ylegacy = firedates_2ylegacy %>% 
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
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(35)))

firedates_3ylegacy = data.frame(year = year(firedates_3ylegacy$effect_date),
                                month = month(firedates_3ylegacy$effect_date),
                                site = firedates_3ylegacy$site,
                                fire_pa_3ylegacy = 1,
                                ws_fire_area_m2_3ylegacy = firedates_3ylegacy$ws_fire_area_m2,
                                fire_perc_ws_3ylegacy = firedates_3ylegacy$fire_perc_ws)

firedates_3ylegacy = firedates_3ylegacy %>% 
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
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(47)))

firedates_4ylegacy = data.frame(year = year(firedates_4ylegacy$effect_date),
                                month = month(firedates_4ylegacy$effect_date),
                                site = firedates_4ylegacy$site,
                                fire_pa_4ylegacy = 1,
                                ws_fire_area_m2_4ylegacy = firedates_4ylegacy$ws_fire_area_m2,
                                fire_perc_ws_4ylegacy = firedates_4ylegacy$fire_perc_ws)

firedates_4ylegacy = firedates_4ylegacy %>% 
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
                           cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(59)))

firedates_5ylegacy = data.frame(year = year(firedates_5ylegacy$effect_date),
                                month = month(firedates_5ylegacy$effect_date),
                                site = firedates_5ylegacy$site,
                                fire_pa_5ylegacy = 1,
                                ws_fire_area_m2_5ylegacy = firedates_5ylegacy$ws_fire_area_m2,
                                fire_perc_ws_5ylegacy = firedates_5ylegacy$fire_perc_ws)

firedates_5ylegacy = firedates_5ylegacy %>% 
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

# set all NAs to equal 0
dat14[,14:31][is.na(dat14[,14:31])] = 0

# plot legacy effects
qplot(index, fire_pa_0.5ylegacy, data=dat14, colour=site, geom="path", facets = "region")
qplot(index, ws_fire_area_m2_1ylegacy, data=dat14, colour=site, geom="path", facets = "region")
qplot(index, fire_perc_ws_1ylegacy, data=dat14, colour=site, geom="path", facets = "region")

# plot legacy effects
qplot(index, fire_pa_5ylegacy, data=dat14, colour=site, geom="path", facets = "region")
qplot(index, ws_fire_area_m2_5ylegacy, data=dat14, colour=site, geom="path", facets = "region")
qplot(index, fire_perc_ws_5ylegacy, data=dat14, colour=site, geom="path", facets = "region")

#### Add fire x precip legacy effects ####

# for fire perc burn ws
dat14$fire_perc_ws_ppt = dat14$fire_perc_ws_0.5ylegacy*dat14$cumulative_precip_mm
dat14$fire_perc_ws_ppt_1ylegacy = dat14$fire_perc_ws_1ylegacy*dat14$cumulative_precip_mm
dat14$fire_perc_ws_ppt_2ylegacy = dat14$fire_perc_ws_2ylegacy*dat14$cumulative_precip_mm
dat14$fire_perc_ws_ppt_3ylegacy = dat14$fire_perc_ws_3ylegacy*dat14$cumulative_precip_mm
dat14$fire_perc_ws_ppt_4ylegacy = dat14$fire_perc_ws_4ylegacy*dat14$cumulative_precip_mm
dat14$fire_perc_ws_ppt_5ylegacy = dat14$fire_perc_ws_5ylegacy*dat14$cumulative_precip_mm

# plot interaction effects
qplot(index, fire_perc_ws_ppt, data=dat14, 
      colour=site, geom="path", facets = "region")
qplot(index, fire_perc_ws_ppt_5ylegacy, data=dat14, 
      colour=site, geom="path", facets = "region")

#### Examine data closely ####

# precip (0-1000 mm)
ggplot(dat14, aes(x = index, y = cumulative_precip_mm, color = site)) +
  geom_point() +
  facet_wrap(.~site) +
  theme(legend.position = "none")

# Sp. Cond. (0-3000 uS)
ggplot(dat14, aes(x = index, y = mean_cond_uScm, color = site)) +
  geom_point() +
  facet_wrap(.~site) +
  theme(legend.position = "none")

# Wow, the VC sites have much lower Sp. Cond., so standardizing makes
# a lot of sense in the model script.

#### Export data ready for modeling ####

saveRDS(dat14, "data_working/marss_data_sb_vc_081423.rds")

# End of script.
