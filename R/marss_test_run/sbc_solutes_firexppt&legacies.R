# SB MARSS models for other solutes
# Script started August 26, 2022
# Heili Lowman, Alex Webster

#### READ ME ####

# The following script will prepare data and run MARSS analyses at SBC sites for the CRASS project.
# MARSS modeling sections each contain 1 model or 1 model comparison and restart the script and import data anew each time.
# Be sure to load libraries before starting modeling sections. 
# Much of this code has been copied over from the "sbc_vcnp_firexppt&legacies.R" script written by Alex, but any of the NM site processing/modeling has been removed.

#### Load packages - ***do this first for all sections!*** ####
library(tidyverse)
library(lubridate)
library(MARSS)
library(naniar) 
library(beepr)

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))

#### Compile data ####

# Load datasets
# Stream Chemistry - all SBC sites
# see script "sbc_chem_vwm.R" for the code used to tidy/generate this dataset
chem_nh4 <- read_csv("data_working/SB_NH4_VWM_011723.csv")
chem_no3 <- read_csv("data_working/SB_NO3_VWM_011723.csv")
chem_po4 <- read_csv("data_working/SB_PO4_VWM_011723.csv")

# Precipitation - all SBC sites
precip <- readRDS("data_working/SBprecip_edited_120121.rds")

# Fire Events - all SBC sites
fire <- readRDS("data_working/SBfire_edited_011323.rds")

# Site Location information - all SBC sites
location <- read_csv("data_raw/sbc_sites_stream_hydro.csv")

# I first need to identify the matching stream sites for the precip data
precip_ed <- precip %>%
  mutate(sitecode_match = factor(case_when(
    sitecode == "GV202" ~ "GV01",
    sitecode == "HO202" ~ "HO00", # switched back from "BARA" since Q data doesn't
    # start until 2002 either, using HO202 instead of 201 because record is longer
    # sitecode == "RG202" ~ "RG01",
    # sitecode == "SMPA" ~ "SP02",
    # sitecode == "GORY" ~ "AT07",
    sitecode == "CAWTP" ~ "AB00",
    sitecode == "STFS" ~ "MC06", # BOGA doesn't start until 2005
    sitecode == "ELDE" ~ "RS02"))) %>%
  dplyr::rename(cumulative_precip_mm = c_precip_mm,
                sitecode_precip = sitecode) %>%
  mutate(Day = 1) %>% # new column of "days"
  mutate(Date = make_date(Year, Month, Day))

# check edited fire data
sitez = c("AB00", "GV01", "HO00", "MC06", "RS02")

ggplot(fire %>% filter(site %in% sitez), 
       aes(x = date, y = fire_pa)) +
  geom_point() +
  labs(x = "Date") +
  facet_wrap(.~site, scales = "free",
             ncol = 1) +
  theme_bw() # looks good

ggplot(fire %>% filter(site %in% sitez), 
       aes(x = date, y = fire_perc_ws)) +
  geom_point() +
  labs(x = "Date") +
  facet_wrap(.~site, scales = "free",
             ncol = 1) +
  theme_bw() # also looks good

## Timeframe Selection

# Examine precip data coverage for MARSS timescale delineation
precip_ed %>%
  drop_na(sitecode_match) %>%
  ggplot(aes(x = Year, y = sitecode_match, color = sitecode_match)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none")

# So, covariate data, which cannot be missing, can run from a maximum of 9/2002
# to 10/2017, which will slightly increase the amount of data used.

# To avoid strange gaps in data, I'm going to start with the fire data,
# since I know it extends from 9/1/2002 to 10/1/2017. This will ensure all other 
# data are joined to these dates in full (which was causing problems earlier).
fire_precip <- left_join(fire, precip_ed, by = c("site" = "sitecode_match", "date" = "Date"))

fire_precip <- fire_precip %>%
  mutate(year = year(date),
         month = month(date)) # and the Year/Month didn't populate, 
# so adding in new columns

# Trim down to sites of interest.
dat <- fire_precip %>%
  filter(site %in% c("AB00", "GV01", "HO00", "MC06", "RS02"))

# Then, left join with chemistry so as not to lose any data.
dat_nh4 <- left_join(dat, chem_nh4, by = c("site" = "site_code", 
                                                   "year" = "Year", 
                                                   "month" = "Month")) %>%
  dplyr::rename("cum_Q_L_nh4" = "cum_Q_L",
         "c_x_Q_nh4" = "c_x_Q",
         "cum_Q_storms_nh4" = "cum_Q_storms",
         "vwm_nh4" = "vwm_monthly")

dat_no3 <- left_join(dat_nh4, chem_no3, by = c("site" = "site_code", 
                                           "year" = "Year", 
                                           "month" = "Month")) %>%
  dplyr::rename("cum_Q_L_no3" = "cum_Q_L",
                "c_x_Q_no3" = "c_x_Q",
                "cum_Q_storms_no3" = "cum_Q_storms",
                "vwm_no3" = "vwm_monthly")

dat_po4 <- left_join(dat_no3, chem_po4, by = c("site" = "site_code", 
                                               "year" = "Year", 
                                               "month" = "Month")) %>%
  dplyr::rename("cum_Q_L_po4" = "cum_Q_L",
                "c_x_Q_po4" = "c_x_Q",
                "cum_Q_storms_po4" = "cum_Q_storms",
                "vwm_po4" = "vwm_monthly")

dat_full <- dat_po4

# Adding in dummy covariates by season
n_months <- dat_full$date %>%
  unique() %>%
  length()

seas_1 <- sin(2 * pi * seq(n_months) / 12)
seas_2 <- cos(2 * pi * seq(n_months) / 12)

dat_full <- dat_full %>%
  mutate(Season1 = rep(seas_1, 5),
         Season2 = rep(seas_2, 5),
         index = rep(seq(1,182), 5))

# AJW: replace NaNs with NAs
dat_full[is.nan(dat_full)] = NA

# And, inspect dataset for missing covariate data.
sum(is.na(dat_full$fire)) # 0
sum(is.na(dat_full$cumulative_precip_mm)) # 3
# Need to trim out October 2017.

dat_full <- dat_full %>%
  filter(date != as_date(as.character("2017-10-01")))

# Trim down to columns of interest
dat_select <- dat_full %>%
  select(year, month, index, site, ws_area_m2,
         cumulative_precip_mm, 
         fire_ID, ig_date, fire_pa, ws_fire_area_m2, fire_perc_ws,
         vwm_nh4, vwm_no3, vwm_po4, 
         Season1, Season2) %>%
  mutate(region = "SB")

# And export to save progress
#saveRDS(dat_select, "data_working/marss_data_sb_N_P_vwm_011723.rds")

#### Import compiled data and add fire x ppt interactions and legacy effects ####

dat1 <- readRDS("data_working/marss_data_sb_N_P_vwm_011723.rds")

# reorganize
names(dat1)
dat2 = dat1[,c(1,2,3,17,4,5, #"year" "month" "index" "region" "site" "ws_area_m2"
               6, # cumulative_precip_mm
               7:11, # fire info
               15,16, # seasonal effects
               12:14)] # analytes

qplot(index, fire_pa, data=dat2, colour=site, geom="path", facets = "region")

# create df of sites x fire info
firez = unique(dat1[,c(4,7,8,10,11)])
firez = firez[complete.cases(firez),] # removes NAs
firez$ig_date = gsub("/","-",firez$ig_date) # reformat date
firez$year = year(as.Date(firez$ig_date)) # create year column
firez$month = month(as.Date(firez$ig_date)) # create month column
firez$date = as.Date(paste(firez$year, firez$month, "01", sep="-")) # reformat date again

#### Add fire persistence legacy effects ####

# This code is re-run as of 1/17/23.

# Create new effect date column.
firez$effect_date = firez$date

# 0.5 year/6mo legacy (this is to allow a window for the fire x ppt interaction only)
firedates_0.5ylegacy = rbind(firez, 
                             cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(1)),
                             cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(2)),
                             cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(3)),
                             cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(4)),
                             cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(5)))

firedates_0.5ylegacy = data.frame(
  year = year(firedates_0.5ylegacy$effect_date),
  month = month(firedates_0.5ylegacy$effect_date),
  site = firedates_0.5ylegacy$site,
  fire_pa_0.5ylegacy = 1,
  ws_fire_area_m2_0.5ylegacy = firedates_0.5ylegacy$ws_fire_area_m2,
  fire_perc_ws_0.5ylegacy = firedates_0.5ylegacy$fire_perc_ws)

firedates_0.5ylegacy = firedates_0.5ylegacy %>% 
  group_by(year, month, site, fire_pa_0.5ylegacy) %>% 
  dplyr::summarise(ws_fire_area_m2_0.5ylegacy = sum(ws_fire_area_m2_0.5ylegacy),
  fire_perc_ws_0.5ylegacy = sum(fire_perc_ws_0.5ylegacy)) %>%
  ungroup()

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

firedates_1ylegacy = firedates_1ylegacy %>% 
  group_by(year, month, site, fire_pa_1ylegacy) %>% 
  dplyr::summarise(ws_fire_area_m2_1ylegacy = sum(ws_fire_area_m2_1ylegacy),
            fire_perc_ws_1ylegacy = sum(fire_perc_ws_1ylegacy)) %>%
  ungroup()

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
  dplyr::summarise(ws_fire_area_m2_2ylegacy = sum(ws_fire_area_m2_2ylegacy),
            fire_perc_ws_2ylegacy = sum(fire_perc_ws_2ylegacy)) %>%
  ungroup()

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
  dplyr::summarise(ws_fire_area_m2_3ylegacy = sum(ws_fire_area_m2_3ylegacy),
            fire_perc_ws_3ylegacy = sum(fire_perc_ws_3ylegacy)) %>%
  ungroup()

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
  dplyr::summarise(ws_fire_area_m2_4ylegacy = sum(ws_fire_area_m2_4ylegacy),
            fire_perc_ws_4ylegacy = sum(fire_perc_ws_4ylegacy)) %>%
  ungroup()

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
  dplyr::summarise(ws_fire_area_m2_5ylegacy = sum(ws_fire_area_m2_5ylegacy),
            fire_perc_ws_5ylegacy = sum(fire_perc_ws_5ylegacy)) %>%
  ungroup()

# join legacy dates to data
dat3 = left_join(dat2, firedates_0.5ylegacy, by=c("year","month","site"))
dat4 = left_join(dat3, firedates_1ylegacy, by=c("year","month","site"))
dat5 = left_join(dat4, firedates_2ylegacy, by=c("year","month","site"))
dat6 = left_join(dat5, firedates_3ylegacy, by=c("year","month","site"))
dat7 = left_join(dat6, firedates_4ylegacy, by=c("year","month","site"))
dat8 = left_join(dat7, firedates_5ylegacy, by=c("year","month","site"))

# replace NAs in fire columns with 0s
dat8[,18:35][is.na(dat8[,18:35])] = 0

# plot 1 yr legacy effects
qplot(index, fire_pa_0.5ylegacy, data=dat8, colour=site, geom="path", facets = "region")
qplot(index, ws_fire_area_m2_1ylegacy, data=dat8, colour=site, geom="path", facets = "region")
qplot(index, fire_perc_ws_1ylegacy, data=dat8, colour=site, geom="path", facets = "region")

# plot 5 yr legacy effects
qplot(index, fire_pa_5ylegacy, data=dat8, colour=site, geom="path", facets = "region")
qplot(index, ws_fire_area_m2_5ylegacy, data=dat8, colour=site, geom="path", facets = "region")
qplot(index, fire_perc_ws_5ylegacy, data=dat8, colour=site, geom="path", facets = "region")
# looking good!

#### Add fire x ppt legacy effects ####

# for fire pa
dat8$fire_pa_ppt = dat8$fire_pa_0.5ylegacy*dat8$cumulative_precip_mm
dat8$fire_pa_ppt_1ylegacy = dat8$fire_pa_1ylegacy*dat8$cumulative_precip_mm
dat8$fire_pa_ppt_2ylegacy = dat8$fire_pa_2ylegacy*dat8$cumulative_precip_mm
dat8$fire_pa_ppt_3ylegacy = dat8$fire_pa_3ylegacy*dat8$cumulative_precip_mm
dat8$fire_pa_ppt_4ylegacy = dat8$fire_pa_4ylegacy*dat8$cumulative_precip_mm
dat8$fire_pa_ppt_5ylegacy = dat8$fire_pa_5ylegacy*dat8$cumulative_precip_mm

# for fire area
dat8$ws_fire_area_m2_ppt = dat8$ws_fire_area_m2_0.5ylegacy*dat8$cumulative_precip_mm
dat8$ws_fire_area_m2_ppt_1ylegacy = dat8$ws_fire_area_m2_1ylegacy*dat8$cumulative_precip_mm
dat8$ws_fire_area_m2_ppt_2ylegacy = dat8$ws_fire_area_m2_2ylegacy*dat8$cumulative_precip_mm
dat8$ws_fire_area_m2_ppt_3ylegacy = dat8$ws_fire_area_m2_3ylegacy*dat8$cumulative_precip_mm
dat8$ws_fire_area_m2_ppt_4ylegacy = dat8$ws_fire_area_m2_4ylegacy*dat8$cumulative_precip_mm
dat8$ws_fire_area_m2_ppt_5ylegacy = dat8$ws_fire_area_m2_5ylegacy*dat8$cumulative_precip_mm

# for fire perc burn ws
dat8$fire_perc_ws_ppt = dat8$fire_perc_ws_0.5ylegacy*dat8$cumulative_precip_mm
dat8$fire_perc_ws_ppt_1ylegacy = dat8$fire_perc_ws_1ylegacy*dat8$cumulative_precip_mm
dat8$fire_perc_ws_ppt_2ylegacy = dat8$fire_perc_ws_2ylegacy*dat8$cumulative_precip_mm
dat8$fire_perc_ws_ppt_3ylegacy = dat8$fire_perc_ws_3ylegacy*dat8$cumulative_precip_mm
dat8$fire_perc_ws_ppt_4ylegacy = dat8$fire_perc_ws_4ylegacy*dat8$cumulative_precip_mm
dat8$fire_perc_ws_ppt_5ylegacy = dat8$fire_perc_ws_5ylegacy*dat8$cumulative_precip_mm

# plot 6 mo interaction effects
qplot(index, fire_pa_ppt, data=dat8, colour=site, geom="path", facets = "region")
qplot(index, ws_fire_area_m2_ppt, data=dat8, colour=site, geom="path", facets = "region")
qplot(index, fire_perc_ws_ppt, data=dat8, colour=site, geom="path", facets = "region")

# plot 1 yr interaction effects
qplot(index, fire_pa_ppt_1ylegacy, data=dat8, colour=site, geom="path", facets = "region")
qplot(index, ws_fire_area_m2_ppt_1ylegacy, data=dat8, colour=site, geom="path", facets = "region")
qplot(index, fire_perc_ws_ppt_1ylegacy, data=dat8, colour=site, geom="path", facets = "region")
# yay, this looks fantastic :)

#### Examine data closely ####

# precip
ggplot(dat8, aes(x = index, y = cumulative_precip_mm, color = site)) +
  geom_point() +
  facet_wrap(.~site) +
  theme(legend.position = "none")

# NH4 (0-60 uM)
ggplot(dat8, aes(x = index, y = vwm_nh4, color = site)) +
  geom_point() +
  facet_wrap(.~site) +
  theme(legend.position = "none")

# NO3 (0-600 uM)
ggplot(dat8, aes(x = index, y = vwm_no3, color = site)) +
  geom_point() +
  facet_wrap(.~site) +
  theme(legend.position = "none")

# PO4 (0-15 uM)
ggplot(dat8, aes(x = index, y = vwm_po4, color = site)) +
  geom_point() +
  facet_wrap(.~site) +
  theme(legend.position = "none")

#### Remove outliers ####

# Examined chem_reg file from the sbc_chem_compilation.R script to compare weekly solute values to
# measured monthly means to make sure the values made sense and they appear to.

# Raw, unaveraged values for the watersheds in question ranged from:
# 0 - 600 uM for NH4 (GV01 has max value)
# 0 - 2000 uM for NO3 (RS02 has max value, with RG01 close behind)
# 0 - 120 uM for PO4 (AT07 has max value)

# In most cases the monthly averages were a dampened mirror of the trends displayed in the raw values.
# A few additional notes:
# The fire effect on NH4 in Gaviota appears VERY strong in the raw data.
# There is a period towards the end of the RG01 PO4 data where there are less values, so one "monthly"
# value is in fact just a raw data point, so it seems a bit higher/like an outlier, but I'm choosing to
# keep it in since it's still a real data point.

# Keeping all data points for now.

#### Export data ready for modeling ####

saveRDS(dat8, "data_working/marss_data_sb_011723.rds")

#### Plot data ready for modeling ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions
dat <- readRDS("data_working/marss_data_sb_011723.rds")

# add date
dat$date = as.Date(paste(dat$year, dat$month, "01", sep="-"))

# examine data once more to be sure sites are appropriate to model
# NH4 in all sites
(nh4fig <- ggplot(dat %>% filter(site %in% c("AB00", "GV01", 
                                             "HO00", "MC06", 
                                             "RS02")), 
                  aes(x = date, y = vwm_nh4)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red"))

# Export for CRASS slides
# ggsave(("SB_NH4_011723.png"),
#        path = "figures",
#        width = 30,
#        height = 15,
#        units = "cm"
#        )

# NO3 in all sites
(no3fig <- ggplot(dat %>% filter(site %in% c("AB00", "GV01", 
                                             "HO00", "MC06", 
                                             "RS02")), 
                  aes(x = date, y = vwm_no3)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red"))

# Export for CRASS slides
# ggsave(("SB_NO3_011723.png"),
#        path = "figures",
#        width = 30,
#        height = 15,
#        units = "cm"
#        )

# PO4 in all sites
(po4fig <- ggplot(dat %>% filter(site %in% c("AB00", "GV01", 
                                             "HO00", "MC06", 
                                             "RS02")), 
                 aes(x = date, y = vwm_po4)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red"))

# Export for CRASS slides
# ggsave(("SB_PO4_011723.png"),
#        path = "figures",
#        width = 30,
#        height = 15,
#        units = "cm"
#        )

# select sites for modeling
# these have the longest most complete ts for NH4/NO3/PO4 and have coverage before and after fires: 
# AB00, GV01, HO00, MC06, RS02
sitez = c("AB00", "GV01", "HO00", "MC06", "RS02")
dat = dat[dat$site %in% sitez,]
# AT07 was removed at John's request and due to lack of Q data
# RG01 has a fire on the last data point
# SP02 has only 6 pre-fire data points

# Export data for later.
saveRDS(dat, "data_working/marss_data_sb_5sites_011723.rds")

# Examine remaining covariates with the filtered dataset

# Fire area in all sites - check to be sure lines up with fire dates
ggplot(dat, aes(x = date, y = ws_fire_area_m2)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date),
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")
  
# Fire % of watershed burned in all sites - check to be sure lines up with fire dates
ggplot(dat, aes(x = date, y = fire_perc_ws)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")

# ppt in all sites - check for outliers
ggplot(dat, aes(x = date, y = cumulative_precip_mm)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")

# Fire p/a X ppt in all sites - check if ppt interacts w/ fire in 6 mo timeframe
ggplot(dat, aes(x = date, y = fire_pa_ppt)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")

# Fire area X ppt in all sites - again check for interaction with precip
ggplot(dat, aes(x = date, y = ws_fire_area_m2_ppt)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")

# Fire % watershed burned X ppt in all sites - again check for interaction
ggplot(dat, aes(x = date, y = fire_perc_ws_ppt)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")

#### MARSS NH4: ppt, fire pa (2m win), fire pa (6m win) x ppt, no legacy effects ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_082622.rds")

# select sites
# include these sites only (6 total - these have the longest most complete ts):
# AB00, AT07, GV01, HO00, MC06, RS02 = SB
sitez = c("AB00", "AT07", "GV01", "HO00", "MC06", "RS02")
dat = dat[dat$site %in% sitez,]
table(dat$site)
table(dat$site,dat$fire_pa)

# pivot wider for MARSS format
dat_nh4 <- dat %>%
  select(
    site, index, 
    mean_nh4_uM, 
    cumulative_precip_mm, 
    fire_pa, fire_pa_6m_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_nh4_uM, cumulative_precip_mm, 
                    fire_pa, fire_pa_6m_ppt))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:7)
cov_cols = c(8:25)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 264
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# yes

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,0,
# 0,"XXXX_AT07",0,0,0,0,
# 0,0,"XXXX_GV01",0,0,0,
# 0,0,0,"XXXX_HO00",0,0,
# 0,0,0,0,"XXXX_MC06",0,
# 0,0,0,0,0,"XXXX_RS02"

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,
  0,"cumulative_precip_mm_AT07",0,0,0,0,
  0,0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_pa_AB00",0,0,0,0,0,
  0,"fire_pa_AT07",0,0,0,0,
  0,0,"fire_pa_GV01",0,0,0,
  0,0,0,"fire_pa_HO00",0,0,
  0,0,0,0,"fire_pa_MC06",0,
  0,0,0,0,0,"fire_pa_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_pa_6m_ppt_AB00",0,0,0,0,0,
  0,"fire_pa_6m_ppt_AT07",0,0,0,0,
  0,0,"fire_pa_6m_ppt_GV01",0,0,0,
  0,0,0,"fire_pa_6m_ppt_HO00",0,0,
  0,0,0,0,"fire_pa_6m_ppt_MC06",0,
  0,0,0,0,0,"fire_pa_6m_ppt_RS02"), 6,18)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_08262022_6state_nh4_fire2mpa_fire6mpaxPPT_nolegacies_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_08262022_6state_nh4_fire2mpa_fire6mpaxPPT_nolegacies_mBFGS.rds")

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

#           dAIC df
# fit        0.0 36
# null.fit 68.7 18
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No and Yes

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. I mean some fall outside the CIs, but overall look ok.

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Meh, for the most part. HO00 looks worst.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00","AT07", "GV01","HO00","MC06","RS02")

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

CIs_fit_ed$Parameter = rep(c("Cum. Ppt", "Fire p/a","Cum. Ppt * Fire p/a (6m)"),6)

CIs_fit_ed$Region = c(rep(c("SB"),6*3)) # **CHECK ORDER OF SITES FOR THIS!!**

# plot results
(RESULTS_ALL_nh4 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Ammonium (NH4) MARSS modeling results - 08/26/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

#### MARSS NH4: ppt, % burn (2m win), % burn (6m win) x ppt, no legacy effects ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_5sites_011723.rds")

# pivot wider for MARSS format
# selecting interaction term with no lag
dat_nh4 <- dat %>%
  select(
    site, index, 
    vwm_nh4, 
    cumulative_precip_mm, 
    fire_perc_ws, fire_perc_ws_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_nh4, cumulative_precip_mm, 
                    fire_perc_ws, fire_perc_ws_ppt))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:6)
cov_cols = c(7:21)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 299
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Yay!!

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,
# 0,"XXXX_GV01",0,0,0,
# 0,0,"XXXX_HO00",0,0,
# 0,0,0,"XXXX_MC06",0,
# 0,0,0,0,"XXXX_RS02"

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_perc_ws_AB00",0,0,0,0,
  0,"fire_perc_ws_GV01",0,0,0,
  0,0,"fire_perc_ws_HO00",0,0,
  0,0,0,"fire_perc_ws_MC06",0,
  0,0,0,0,"fire_perc_ws_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_AB00",0,0,0,0,
  0,"fire_perc_ws_ppt_GV01",0,0,0,
  0,0,"fire_perc_ws_ppt_HO00",0,0,
  0,0,0,"fire_perc_ws_ppt_MC06",0,
  0,0,0,0,"fire_perc_ws_ppt_RS02"), 5,15)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_011723_5state_nh4_percburn2m_percburn6mxppt_nolegacies_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_011723_5state_nh4_percburn2m_percburn6mxppt_nolegacies_mBFGS.rds")

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

#           dAIC df
# fit        0.0 30
# null.fit  77.5 15
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yess!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks fine - keep outliers in.

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Meh, for the most part. Some have a strange curve around theoretical quantile = 1. AB00 looks best, HO00 looks worst.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns. HO00 again has a few that cross the threshold but leaving them for now.

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00", "GV01", "HO00", "MC06", "RS02")

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

CIs_fit_ed$Parameter = rep(c("Cum. Ppt", "% Ws Burned (2m)","Cum. Ppt * % Ws Burned (6m)"),5)

CIs_fit_ed$Region = c(rep(c("SB"),5*3)) # **CHECK ORDER OF SITES FOR THIS!!**

# plot results
(RESULTS_ALL_nh4 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Ammonium (NH4) MARSS modeling results - 01/17/2023\nImmediate Effects") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

# interesting that the less developed sites immediately have responses in the interaction term.

# and save out plot
# ggsave(("figures/MARSS_SB_5state_nh4_011723.png"),
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### MARSS NH4: ppt, % burn (1y legacy), % burn x ppt 1y legacy ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_5sites_011723.rds")

# pivot wider for MARSS format
dat_nh4 <- dat %>%
  select(
    site, index, 
    vwm_nh4, 
    cumulative_precip_mm, 
    fire_perc_ws_1ylegacy, fire_perc_ws_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(vwm_nh4, cumulative_precip_mm, 
                    fire_perc_ws_1ylegacy, fire_perc_ws_ppt_1ylegacy))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:6)
cov_cols = c(7:21)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 299
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,
# 0,"XXXX_GV01",0,0,0,
# 0,0,"XXXX_HO00",0,0,
# 0,0,0,"XXXX_MC06",0,
# 0,0,0,0,"XXXX_RS02"

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_1ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_1ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_1ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_1ylegacy_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02"),5,15)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_011723_5state_nh4_percburn1ylegacy_percburnxppt1ylegacy_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_011723_5state_nh4_percburn1ylegacy_percburnxppt1ylegacy_mBFGS.rds")

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

#           dAIC df
# fit        0.0 30
# null.fit 136.9 15
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok! Similar to past model results.

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look ok - better than immediate fire effect qq plots. HO00 still wonky though.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00", "GV01", "HO00", "MC06", "RS02")

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

CIs_fit_ed$Parameter = rep(c("Cum. Ppt", "% Ws Burned (1y)","Cum. Ppt * % Ws Burned (1y)"),5)

CIs_fit_ed$Region = c(rep(c("SB"),5*3))

# plot results
(RESULTS_ALL_nh4 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Ammonium (NH4) MARSS modeling results - 01/17/2023\n1 Year Legacy") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

# and save out plot
# ggsave(("figures/MARSS_SB_5state_nh4_1y_011723.png"),
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### MARSS NH4: ppt, % burn (2y legacy), % burn x ppt 2y legacy ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_6sites_090722.rds")

# NOTE: REMOVE AT07 SITE FROM 2,3,4,&5y LEGACY MODELS - not enough data.

# pivot wider for MARSS format
dat_nh4 <- dat %>%
  # remove AT07 site
  filter(site != "AT07") %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    mean_nh4_uM, 
    cumulative_precip_mm, 
    fire_perc_ws_2ylegacy, fire_perc_ws_ppt_2ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_nh4_uM, cumulative_precip_mm, 
                    fire_perc_ws_2ylegacy, fire_perc_ws_ppt_2ylegacy))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:6)
cov_cols = c(7:21)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 203
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,
# 0,"XXXX_GV01",0,0,0,
# 0,0,"XXXX_HO00",0,0,
# 0,0,0,"XXXX_MC06",0,
# 0,0,0,0,"XXXX_RS02""

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_perc_ws_2ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_2ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_2ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_2ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_2ylegacy_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_2ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_ppt_2ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_ppt_2ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_ppt_2ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_ppt_2ylegacy_RS02"),5,15)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_09072022_5state_nh4_percburn2ylegacy_percburnxppt2ylegacy_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_nh4_percburn2ylegacy_percburnxppt2ylegacy_mBFGS.rds")

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

#           dAIC df
# fit        0.0 30
# null.fit  40.6 15
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Check with Alex about the clustering happening at GV01/HO00/MC06?

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look ok. MC06 looks worst.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00", "GV01", "HO00", "MC06", "RS02")

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
                             "% Ws Burned (2y)",
                             "Cum. Ppt * % Ws Burned (2y)"),5)

CIs_fit_ed$Region = c(rep(c("SB"),5*3))

# plot results
(RESULTS_ALL_nh4 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Ammonium (NH4) MARSS modeling results - 09/07/2022\n2 Year Legacy") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

# and save out plot
# ggsave(("figures/MARSS_SB_5state_nh4_2y_090722.png"),
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### MARSS NH4: ppt, % burn (3y legacy), % burn x ppt 3y legacy ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_6sites_090722.rds")

# NOTE: REMOVE AT07 SITE FROM 2,3,4,&5y LEGACY MODELS - not enough data.

# pivot wider for MARSS format
dat_nh4 <- dat %>%
  # remove AT07 site
  filter(site != "AT07") %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    mean_nh4_uM, 
    cumulative_precip_mm, 
    fire_perc_ws_3ylegacy, fire_perc_ws_ppt_3ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_nh4_uM, cumulative_precip_mm, 
                    fire_perc_ws_3ylegacy, fire_perc_ws_ppt_3ylegacy))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:6)
cov_cols = c(7:21)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 203
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,
# 0,"XXXX_GV01",0,0,0,
# 0,0,"XXXX_HO00",0,0,
# 0,0,0,"XXXX_MC06",0,
# 0,0,0,0,"XXXX_RS02""

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_perc_ws_3ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_3ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_3ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_3ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_3ylegacy_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_3ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_ppt_3ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_ppt_3ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_ppt_3ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_ppt_3ylegacy_RS02"),5,15)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_09072022_5state_nh4_percburn3ylegacy_percburnxppt3ylegacy_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_nh4_percburn3ylegacy_percburnxppt3ylegacy_mBFGS.rds")

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

#           dAIC df
# fit        0.0 30
# null.fit  38.1 15
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Check with Alex about the clustering happening at GV01/HO00/MC06?

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look ok. What was previously a dent is now turning into a horizontal line sometimes.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00", "GV01", "HO00", "MC06", "RS02")

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
                             "% Ws Burned (3y)",
                             "Cum. Ppt * % Ws Burned (3y)"),5)

CIs_fit_ed$Region = c(rep(c("SB"),5*3))

# plot results
(RESULTS_ALL_nh4 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Ammonium (NH4) MARSS modeling results - 09/07/2022\n3 Year Legacy") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

# and save out plot
# ggsave(("figures/MARSS_SB_5state_nh4_3y_090722.png"),
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### MARSS NH4: ppt, % burn (4y legacy), % burn x ppt 4y legacy ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_6sites_090722.rds")

# NOTE: REMOVE AT07 SITE FROM 2,3,4,&5y LEGACY MODELS - not enough data.

# pivot wider for MARSS format
dat_nh4 <- dat %>%
  # remove AT07 site
  filter(site != "AT07") %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    mean_nh4_uM, 
    cumulative_precip_mm, 
    fire_perc_ws_4ylegacy, fire_perc_ws_ppt_4ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_nh4_uM, cumulative_precip_mm, 
                    fire_perc_ws_4ylegacy, fire_perc_ws_ppt_4ylegacy))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:6)
cov_cols = c(7:21)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 203
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,
# 0,"XXXX_GV01",0,0,0,
# 0,0,"XXXX_HO00",0,0,
# 0,0,0,"XXXX_MC06",0,
# 0,0,0,0,"XXXX_RS02""

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_perc_ws_4ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_4ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_4ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_4ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_4ylegacy_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_4ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_ppt_4ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_ppt_4ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_ppt_4ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_ppt_4ylegacy_RS02"),5,15)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_09072022_5state_nh4_percburn4ylegacy_percburnxppt4ylegacy_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_nh4_percburn4ylegacy_percburnxppt4ylegacy_mBFGS.rds")

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

#           dAIC df
# fit        0.0 30
# null.fit  39.5 15
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Check with Alex about the clustering happening at GV01/HO00/MC06?

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look ok. Now HO00 looks worse than MC06.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00", "GV01", "HO00", "MC06", "RS02")

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
                             "% Ws Burned (4y)",
                             "Cum. Ppt * % Ws Burned (4y)"),5)

CIs_fit_ed$Region = c(rep(c("SB"),5*3))

# plot results
(RESULTS_ALL_nh4 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Ammonium (NH4) MARSS modeling results - 09/07/2022\n4 Year Legacy") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

# and save out plot
# ggsave(("figures/MARSS_SB_5state_nh4_4y_090722.png"),
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### MARSS NH4: ppt, % burn (5y legacy), % burn x ppt 5y legacy ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_6sites_090722.rds")

# NOTE: REMOVE AT07 SITE FROM 2,3,4,&5y LEGACY MODELS - not enough data.

# pivot wider for MARSS format
dat_nh4 <- dat %>%
  # remove AT07 site
  filter(site != "AT07") %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    mean_nh4_uM, 
    cumulative_precip_mm, 
    fire_perc_ws_5ylegacy, fire_perc_ws_ppt_5ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_nh4_uM, cumulative_precip_mm, 
                    fire_perc_ws_5ylegacy, fire_perc_ws_ppt_5ylegacy))

# indicate column #s of response and predictor vars
names(dat_nh4)
resp_cols = c(2:6)
cov_cols = c(7:21)

# log and scale transform response var
dat_nh4_log = dat_nh4
dat_nh4_log[,resp_cols] = log10(dat_nh4_log[,resp_cols])
dat_nh4_log[,resp_cols] = scale(dat_nh4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_nh4_log[,resp_cols])) # 0
sum(is.na(dat_nh4_log[,resp_cols])) # 203
range(dat_nh4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_nh4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_nh4_log[,cov_cols])) # 0
sum(is.na(dat_nh4_log[,cov_cols])) # 0
sum(is.infinite(dat_nh4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_nh4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,
# 0,"XXXX_GV01",0,0,0,
# 0,0,"XXXX_HO00",0,0,
# 0,0,0,"XXXX_MC06",0,
# 0,0,0,0,"XXXX_RS02""

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_perc_ws_5ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_5ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_5ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_5ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_5ylegacy_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_5ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_ppt_5ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_ppt_5ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_ppt_5ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_ppt_5ylegacy_RS02"),5,15)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_09072022_5state_nh4_percburn5ylegacy_percburnxppt5ylegacy_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_nh4_percburn5ylegacy_percburnxppt5ylegacy_mBFGS.rds")

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

#           dAIC df
# fit        0.0 30
# null.fit  40.3 15
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Check with Alex about the clustering happening at GV01/HO00/MC06?

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look ok. HO00 still looks worst.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00", "GV01", "HO00", "MC06", "RS02")

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
                             "% Ws Burned (5y)",
                             "Cum. Ppt * % Ws Burned (5y)"),5)

CIs_fit_ed$Region = c(rep(c("SB"),5*3))

# plot results
(RESULTS_ALL_nh4 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Ammonium (NH4) MARSS modeling results - 09/07/2022\n5 Year Legacy") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

# and save out plot
# ggsave(("figures/MARSS_SB_5state_nh4_5y_090722.png"),
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### Summary fig of NH4 results ####

# import model fits
noleg_fit = readRDS(file = "data_working/marss_test_run/fit_09072022_6state_nh4_percburn2m_percburn6mxppt_nolegacies_mBFGS.rds")
leg1y_fit = readRDS(file = "data_working/marss_test_run/fit_09072022_6state_nh4_percburn1ylegacy_percburnxppt1ylegacy_mBFGS.rds")
leg2y_fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_nh4_percburn2ylegacy_percburnxppt2ylegacy_mBFGS.rds") # dropped AT07
leg3y_fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_nh4_percburn3ylegacy_percburnxppt3ylegacy_mBFGS.rds")
leg4y_fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_nh4_percburn4ylegacy_percburnxppt4ylegacy_mBFGS.rds")
leg5y_fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_nh4_percburn5ylegacy_percburnxppt5ylegacy_mBFGS.rds")

# extract necessary confidence interval info
noleg_est <- MARSSparamCIs(noleg_fit)
leg1y_est <- MARSSparamCIs(leg1y_fit)
leg2y_est <- MARSSparamCIs(leg2y_fit)
leg3y_est <- MARSSparamCIs(leg3y_fit)
leg4y_est <- MARSSparamCIs(leg4y_fit)
leg5y_est <- MARSSparamCIs(leg5y_fit)

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
CIs = rbind(
  noleg_CI, leg1y_CI, leg2y_CI, leg3y_CI, leg4y_CI, leg5y_CI
)

# add col for site names
CIs$Stream = gsub("_","",str_sub(CIs$Parameter, start= -4))

# simplify parm names
CIs$Parm_simple = c(rep("Ppt",6),
                    rep("Perc. burn",6),
                    rep("Ppt x Perc. burn",6),
                    
                    rep("Ppt",6),
                    rep("Perc. burn",6),
                    rep("Ppt x Perc. burn",6),
                    
                    rep("Ppt",5),
                    rep("Perc. burn",5),
                    rep("Ppt x Perc. burn",5),
                    
                    rep("Ppt",5),
                    rep("Perc. burn",5),
                    rep("Ppt x Perc. burn",5),
                    
                    rep("Ppt",5),
                    rep("Perc. burn",5),
                    rep("Ppt x Perc. burn",5),
                    
                    rep("Ppt",5),
                    rep("Perc. burn",5),
                    rep("Ppt x Perc. burn",5))

# plot results
(RESULTS_ALL_d <- ggplot(CIs, aes(x=factor(Parm_simple, 
                                           levels = c("Ppt x Perc. burn","Perc. burn","Ppt")), 
                                  Est., color=Stream)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0) +
    geom_point(position=position_dodge(width=0.3), size=3) + 
    theme_bw()+
    theme(plot.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20)) +
    geom_hline(aes(yintercept=0), linetype="dashed") +
    coord_flip(ylim = c(-1.25, 1.25)) + 
    labs(y = "", 
         x="",
         title = expression(paste(NH[4]^{"+"}," MARSS modeling results for Santa Barbara sites"))) +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(.~Model, scales = "free", nrow = 1))

# Export plot.
# ggsave(("MARSS_SB_NH4_091522.png"),
#        path = "figures",
#        width = 80,
#        height = 10,
#        units = "cm"
#        )

# Save these CIs as NH4 specific, so that we can create an enormous facetted plot at the very end.
CIs_NH4 <- CIs
CIs_NH4$Solute <- "NH4"

#### MARSS NO3: ppt, fire pa (2m win), fire pa (6m win) x ppt, no legacy effects ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_082622.rds")

# select sites
# include these sites only (6 total - these have the longest most complete ts):
# AB00, AT07, GV01, HO00, MC06, RS02 = SB
sitez = c("AB00", "AT07", "GV01", "HO00", "MC06", "RS02")
dat = dat[dat$site %in% sitez,]
table(dat$site)
table(dat$site,dat$fire_pa)

# pivot wider for MARSS format
dat_no3 <- dat %>%
  select(
    site, index, 
    mean_no3_uM, 
    cumulative_precip_mm, 
    fire_pa, fire_pa_6m_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_no3_uM, cumulative_precip_mm, 
                    fire_pa, fire_pa_6m_ppt))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:7)
cov_cols = c(8:25)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 256
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# yes

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,0,
# 0,"XXXX_AT07",0,0,0,0,
# 0,0,"XXXX_GV01",0,0,0,
# 0,0,0,"XXXX_HO00",0,0,
# 0,0,0,0,"XXXX_MC06",0,
# 0,0,0,0,0,"XXXX_RS02"

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,
  0,"cumulative_precip_mm_AT07",0,0,0,0,
  0,0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_pa_AB00",0,0,0,0,0,
  0,"fire_pa_AT07",0,0,0,0,
  0,0,"fire_pa_GV01",0,0,0,
  0,0,0,"fire_pa_HO00",0,0,
  0,0,0,0,"fire_pa_MC06",0,
  0,0,0,0,0,"fire_pa_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_pa_6m_ppt_AB00",0,0,0,0,0,
  0,"fire_pa_6m_ppt_AT07",0,0,0,0,
  0,0,"fire_pa_6m_ppt_GV01",0,0,0,
  0,0,0,"fire_pa_6m_ppt_HO00",0,0,
  0,0,0,0,"fire_pa_6m_ppt_MC06",0,
  0,0,0,0,0,"fire_pa_6m_ppt_RS02"), 6,18)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_08262022_6state_no3_fire2mpa_fire6mpaxPPT_nolegacies_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_08262022_6state_no3_fire2mpa_fire6mpaxPPT_nolegacies_mBFGS.rds")

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

#           dAIC df
# fit        0.0 36
# null.fit 109.3 18
# RESULT: covar model is better than null, phew

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall within the CIs? No temporal patterns. Yes, appear to fall within the CIs.

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yessiree!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Some fall outside the CIs, but overall look ok.

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# These look great!!

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No visible patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00","AT07", "GV01","HO00","MC06","RS02")

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

CIs_fit_ed$Parameter = rep(c("Cum. Ppt", "Fire p/a","Cum. Ppt * Fire p/a (6m)"),6)

CIs_fit_ed$Region = c(rep(c("SB"),6*3)) # **CHECK ORDER OF SITES FOR THIS!!**

# plot results
(RESULTS_ALL_no3 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Nitrate (NO3) MARSS modeling results - 08/26/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

#### MARSS NO3: ppt, % burn (2m win), % burn (6m win) x ppt, no legacy effects ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_6sites_090722.rds")

# pivot wider for MARSS format
dat_no3 <- dat %>%
  select(
    site, index, 
    mean_no3_uM, 
    cumulative_precip_mm, 
    fire_perc_ws, fire_perc_ws_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_no3_uM, cumulative_precip_mm, 
                    fire_perc_ws, fire_perc_ws_ppt))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:7)
cov_cols = c(8:25)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 256
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,0,
# 0,"XXXX_AT07",0,0,0,0,
# 0,0,"XXXX_GV01",0,0,0,
# 0,0,0,"XXXX_HO00",0,0,
# 0,0,0,0,"XXXX_MC06",0,
# 0,0,0,0,0,"XXXX_RS02"

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,
  0,"cumulative_precip_mm_AT07",0,0,0,0,
  0,0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_perc_ws_AB00",0,0,0,0,0,
  0,"fire_perc_ws_AT07",0,0,0,0,
  0,0,"fire_perc_ws_GV01",0,0,0,
  0,0,0,"fire_perc_ws_HO00",0,0,
  0,0,0,0,"fire_perc_ws_MC06",0,
  0,0,0,0,0,"fire_perc_ws_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_AB00",0,0,0,0,0,
  0,"fire_perc_ws_ppt_AT07",0,0,0,0,
  0,0,"fire_perc_ws_ppt_GV01",0,0,0,
  0,0,0,"fire_perc_ws_ppt_HO00",0,0,
  0,0,0,0,"fire_perc_ws_ppt_MC06",0,
  0,0,0,0,0,"fire_perc_ws_ppt_RS02"), 6,18)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_09072022_6state_no3_percburn2m_percburn6mxppt_nolegacies_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_09072022_6state_no3_percburn2m_percburn6mxppt_nolegacies_mBFGS.rds")

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

#           dAIC df
# fit        0.0 36
# null.fit 110.6 18
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yess!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# These look good!

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00", "AT07", "GV01", "HO00", "MC06", "RS02")

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

CIs_fit_ed$Parameter = rep(c("Cum. Ppt", "% Ws Burned (2m)","Cum. Ppt * % Ws Burned (6m)"),6)

CIs_fit_ed$Region = c(rep(c("SB"),6*3))

# plot results
(RESULTS_ALL_no3 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Nitrate (NO3) MARSS modeling results - 09/07/2022\nImmediate Effects") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

# and save out plot
# ggsave(("figures/MARSS_SB_6state_no3_090722.png"),
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### MARSS NO3: ppt, % burn (1y legacy), % burn x ppt 1y legacy ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_6sites_090722.rds")

# pivot wider for MARSS format
dat_no3 <- dat %>%
  select(
    site, index, 
    mean_no3_uM, 
    cumulative_precip_mm, 
    fire_perc_ws_1ylegacy, fire_perc_ws_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_no3_uM, cumulative_precip_mm, 
                    fire_perc_ws_1ylegacy, fire_perc_ws_ppt_1ylegacy))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:7)
cov_cols = c(8:25)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 256
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,0,
# 0,"XXXX_AT07",0,0,0,0,
# 0,0,"XXXX_GV01",0,0,0,
# 0,0,0,"XXXX_HO00",0,0,
# 0,0,0,0,"XXXX_MC06",0,
# 0,0,0,0,0,"XXXX_RS02"

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,
  0,"cumulative_precip_mm_AT07",0,0,0,0,
  0,0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_AT07",0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_GV01",0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_HO00",0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_MC06",0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_AT07",0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_MC06",0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02"),6,18)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_09072022_6state_no3_percburn1ylegacy_percburnxppt1ylegacy_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_09072022_6state_no3_percburn1ylegacy_percburnxppt1ylegacy_mBFGS.rds")

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

#           dAIC df
# fit        0.0 36
# null.fit 132.7 18
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Check with Alex about the clustering happening at GV01/HO00/MC06?

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look good!

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00", "AT07", "GV01", "HO00", "MC06", "RS02")

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

CIs_fit_ed$Parameter = rep(c("Cum. Ppt", "% Ws Burned (1y)","Cum. Ppt * % Ws Burned (1y)"),6)

CIs_fit_ed$Region = c(rep(c("SB"),6*3))

# plot results
(RESULTS_ALL_no3 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Nitrate (NO3) MARSS modeling results - 09/07/2022\n1 Year Legacy") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

# and save out plot
# ggsave(("figures/MARSS_SB_6state_no3_1y_090722.png"),
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### MARSS NO3: ppt, % burn (2y legacy), % burn x ppt 2y legacy ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_6sites_090722.rds")

# NOTE: REMOVE AT07 SITE FROM 2,3,4,&5y LEGACY MODELS - not enough data.

# pivot wider for MARSS format
dat_no3 <- dat %>%
  # remove AT07 site
  filter(site != "AT07") %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    mean_no3_uM, 
    cumulative_precip_mm, 
    fire_perc_ws_2ylegacy, fire_perc_ws_ppt_2ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_no3_uM, cumulative_precip_mm, 
                    fire_perc_ws_2ylegacy, fire_perc_ws_ppt_2ylegacy))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:6)
cov_cols = c(7:21)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 196
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,
# 0,"XXXX_GV01",0,0,0,
# 0,0,"XXXX_HO00",0,0,
# 0,0,0,"XXXX_MC06",0,
# 0,0,0,0,"XXXX_RS02""

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_perc_ws_2ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_2ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_2ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_2ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_2ylegacy_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_2ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_ppt_2ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_ppt_2ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_ppt_2ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_ppt_2ylegacy_RS02"),5,15)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_09072022_5state_no3_percburn2ylegacy_percburnxppt2ylegacy_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_no3_percburn2ylegacy_percburnxppt2ylegacy_mBFGS.rds")

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

#           dAIC df
# fit        0.0 30
# null.fit 134.6 15
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Check with Alex about the clustering happening at GV01/HO00/MC06?

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look good.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00", "GV01", "HO00", "MC06", "RS02")

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
                             "% Ws Burned (2y)",
                             "Cum. Ppt * % Ws Burned (2y)"),5)

CIs_fit_ed$Region = c(rep(c("SB"),5*3))

# plot results
(RESULTS_ALL_no3 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Nitrate (NO3) MARSS modeling results - 09/07/2022\n2 Year Legacy") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

# and save out plot
# ggsave(("figures/MARSS_SB_5state_no3_2y_090722.png"),
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### MARSS NO3: ppt, % burn (3y legacy), % burn x ppt 3y legacy ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_6sites_090722.rds")

# NOTE: REMOVE AT07 SITE FROM 2,3,4,&5y LEGACY MODELS - not enough data.

# pivot wider for MARSS format
dat_no3 <- dat %>%
  # remove AT07 site
  filter(site != "AT07") %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    mean_no3_uM, 
    cumulative_precip_mm, 
    fire_perc_ws_3ylegacy, fire_perc_ws_ppt_3ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_no3_uM, cumulative_precip_mm, 
                    fire_perc_ws_3ylegacy, fire_perc_ws_ppt_3ylegacy))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:6)
cov_cols = c(7:21)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 196
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,
# 0,"XXXX_GV01",0,0,0,
# 0,0,"XXXX_HO00",0,0,
# 0,0,0,"XXXX_MC06",0,
# 0,0,0,0,"XXXX_RS02""

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_perc_ws_3ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_3ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_3ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_3ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_3ylegacy_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_3ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_ppt_3ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_ppt_3ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_ppt_3ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_ppt_3ylegacy_RS02"),5,15)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_09072022_5state_no3_percburn3ylegacy_percburnxppt3ylegacy_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_no3_percburn3ylegacy_percburnxppt3ylegacy_mBFGS.rds")

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

#           dAIC df
# fit        0.0 30
# null.fit 131.5 15
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Check with Alex about the clustering happening at GV01/HO00/MC06?

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# All still look good.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00", "GV01", "HO00", "MC06", "RS02")

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
                             "% Ws Burned (3y)",
                             "Cum. Ppt * % Ws Burned (3y)"),5)

CIs_fit_ed$Region = c(rep(c("SB"),5*3))

# plot results
(RESULTS_ALL_no3 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Nitrate (NO3) MARSS modeling results - 09/07/2022\n3 Year Legacy") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

# and save out plot
# ggsave(("figures/MARSS_SB_5state_no3_3y_090722.png"),
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### MARSS NO3: ppt, % burn (4y legacy), % burn x ppt 4y legacy ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_6sites_090722.rds")

# NOTE: REMOVE AT07 SITE FROM 2,3,4,&5y LEGACY MODELS - not enough data.

# pivot wider for MARSS format
dat_no3 <- dat %>%
  # remove AT07 site
  filter(site != "AT07") %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    mean_no3_uM, 
    cumulative_precip_mm, 
    fire_perc_ws_4ylegacy, fire_perc_ws_ppt_4ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_no3_uM, cumulative_precip_mm, 
                    fire_perc_ws_4ylegacy, fire_perc_ws_ppt_4ylegacy))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:6)
cov_cols = c(7:21)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 196
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,
# 0,"XXXX_GV01",0,0,0,
# 0,0,"XXXX_HO00",0,0,
# 0,0,0,"XXXX_MC06",0,
# 0,0,0,0,"XXXX_RS02""

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_perc_ws_4ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_4ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_4ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_4ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_4ylegacy_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_4ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_ppt_4ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_ppt_4ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_ppt_4ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_ppt_4ylegacy_RS02"),5,15)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_09072022_5state_no3_percburn4ylegacy_percburnxppt4ylegacy_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_no3_percburn4ylegacy_percburnxppt4ylegacy_mBFGS.rds")

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

#           dAIC df
# fit        0.0 30
# null.fit 130.1 15
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Check with Alex about the clustering happening at GV01/HO00/MC06?

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# All still look good.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00", "GV01", "HO00", "MC06", "RS02")

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
                             "% Ws Burned (4y)",
                             "Cum. Ppt * % Ws Burned (4y)"),5)

CIs_fit_ed$Region = c(rep(c("SB"),5*3))

# plot results
(RESULTS_ALL_no3 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Nitrate (NO3) MARSS modeling results - 09/07/2022\n4 Year Legacy") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

# and save out plot
# ggsave(("figures/MARSS_SB_5state_no3_4y_090722.png"),
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### MARSS NO3: ppt, % burn (5y legacy), % burn x ppt 5y legacy ####

# Set up data for MARSS
# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_6sites_090722.rds")

# NOTE: REMOVE AT07 SITE FROM 2,3,4,&5y LEGACY MODELS - not enough data.

# pivot wider for MARSS format
dat_no3 <- dat %>%
  # remove AT07 site
  filter(site != "AT07") %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    mean_no3_uM, 
    cumulative_precip_mm, 
    fire_perc_ws_5ylegacy, fire_perc_ws_ppt_5ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_no3_uM, cumulative_precip_mm, 
                    fire_perc_ws_5ylegacy, fire_perc_ws_ppt_5ylegacy))

# indicate column #s of response and predictor vars
names(dat_no3)
resp_cols = c(2:6)
cov_cols = c(7:21)

# log and scale transform response var
dat_no3_log = dat_no3
dat_no3_log[,resp_cols] = log10(dat_no3_log[,resp_cols])
dat_no3_log[,resp_cols] = scale(dat_no3_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_no3_log[,resp_cols])) # 0
sum(is.na(dat_no3_log[,resp_cols])) # 196
range(dat_no3_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_no3_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_no3_log[,cov_cols])) # 0
sum(is.na(dat_no3_log[,cov_cols])) # 0
sum(is.infinite(dat_no3_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_no3_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,
# 0,"XXXX_GV01",0,0,0,
# 0,0,"XXXX_HO00",0,0,
# 0,0,0,"XXXX_MC06",0,
# 0,0,0,0,"XXXX_RS02""

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_perc_ws_5ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_5ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_5ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_5ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_5ylegacy_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_5ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_ppt_5ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_ppt_5ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_ppt_5ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_ppt_5ylegacy_RS02"),5,15)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_09072022_5state_no3_percburn5ylegacy_percburnxppt5ylegacy_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_no3_percburn5ylegacy_percburnxppt5ylegacy_mBFGS.rds")

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

#           dAIC df
# fit        0.0 30
# null.fit 126.9 15
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Check with Alex about the clustering happening at GV01/HO00/MC06?

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# All still look good.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00", "GV01", "HO00", "MC06", "RS02")

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
                             "% Ws Burned (5y)",
                             "Cum. Ppt * % Ws Burned (5y)"),5)

CIs_fit_ed$Region = c(rep(c("SB"),5*3))

# plot results
(RESULTS_ALL_no3 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Nitrate (NO3) MARSS modeling results - 09/07/2022\n5 Year Legacy") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

# and save out plot
# ggsave(("figures/MARSS_SB_5state_no3_5y_090722.png"),
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### Summary fig of NO3 results ####

# import model fits
noleg_fit = readRDS(file = "data_working/marss_test_run/fit_09072022_6state_no3_percburn2m_percburn6mxppt_nolegacies_mBFGS.rds")
leg1y_fit = readRDS(file = "data_working/marss_test_run/fit_09072022_6state_no3_percburn1ylegacy_percburnxppt1ylegacy_mBFGS.rds")
leg2y_fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_no3_percburn2ylegacy_percburnxppt2ylegacy_mBFGS.rds") # dropped AT07
leg3y_fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_no3_percburn3ylegacy_percburnxppt3ylegacy_mBFGS.rds")
leg4y_fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_no3_percburn4ylegacy_percburnxppt4ylegacy_mBFGS.rds")
leg5y_fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_no3_percburn5ylegacy_percburnxppt5ylegacy_mBFGS.rds")

# extract necessary confidence interval info
noleg_est <- MARSSparamCIs(noleg_fit)
leg1y_est <- MARSSparamCIs(leg1y_fit)
leg2y_est <- MARSSparamCIs(leg2y_fit)
leg3y_est <- MARSSparamCIs(leg3y_fit)
leg4y_est <- MARSSparamCIs(leg4y_fit)
leg5y_est <- MARSSparamCIs(leg5y_fit)

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
CIs = rbind(
  noleg_CI, leg1y_CI, leg2y_CI, leg3y_CI, leg4y_CI, leg5y_CI
)

# add col for site names
CIs$Stream = gsub("_","",str_sub(CIs$Parameter, start= -4))

# simplify parm names
CIs$Parm_simple = c(rep("Ppt",6),
                    rep("Perc. burn",6),
                    rep("Ppt x Perc. burn",6),
                    
                    rep("Ppt",6),
                    rep("Perc. burn",6),
                    rep("Ppt x Perc. burn",6),
                    
                    rep("Ppt",5),
                    rep("Perc. burn",5),
                    rep("Ppt x Perc. burn",5),
                    
                    rep("Ppt",5),
                    rep("Perc. burn",5),
                    rep("Ppt x Perc. burn",5),
                    
                    rep("Ppt",5),
                    rep("Perc. burn",5),
                    rep("Ppt x Perc. burn",5),
                    
                    rep("Ppt",5),
                    rep("Perc. burn",5),
                    rep("Ppt x Perc. burn",5))

# plot results
(RESULTS_ALL_e <- ggplot(CIs, aes(x=factor(Parm_simple, 
                                           levels = c("Ppt x Perc. burn","Perc. burn","Ppt")), 
                                  Est., color=Stream)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0) +
    geom_point(position=position_dodge(width=0.3), size=3) + 
    theme_bw()+
    theme(plot.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20)) +
    geom_hline(aes(yintercept=0), linetype="dashed") +
    coord_flip(ylim = c(-1.25, 1.25)) + 
    labs(y = "", 
         x="",
         title = expression(paste(NO[3]^{"-"}," MARSS modeling results for Santa Barbara sites"))) +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(.~Model, scales = "free", nrow = 1))

# Export plot.
# ggsave(("MARSS_SB_NO3_091522.png"),
#        path = "figures",
#        width = 80,
#        height = 10,
#        units = "cm"
#        )

# Save these CIs as NO3 specific, so that we can create an enormous facetted plot at the very end.
CIs_NO3 <- CIs
CIs_NO3$Solute <- "NO3"

#### MARSS PO4: ppt, fire pa (2m win), fire pa (6m win) x ppt, no legacy effects ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat = readRDS("data_working/marss_data_sb_082622.rds")

# select sites
# include these sites only (6 total - these have the longest most complete ts):
# AB00, AT07, GV01, HO00, MC06, RS02 = SB
sitez = c("AB00", "AT07", "GV01", "HO00", "MC06", "RS02")
dat = dat[dat$site %in% sitez,]
table(dat$site)
table(dat$site,dat$fire_pa)

# pivot wider for MARSS format
dat_po4 <- dat %>%
  select(
    site, index, 
    mean_po4_uM, 
    cumulative_precip_mm, 
    fire_pa, fire_pa_6m_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_po4_uM, cumulative_precip_mm, 
                    fire_pa, fire_pa_6m_ppt))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:7)
cov_cols = c(8:25)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 257
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# yes

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,0,
# 0,"XXXX_AT07",0,0,0,0,
# 0,0,"XXXX_GV01",0,0,0,
# 0,0,0,"XXXX_HO00",0,0,
# 0,0,0,0,"XXXX_MC06",0,
# 0,0,0,0,0,"XXXX_RS02"

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,
  0,"cumulative_precip_mm_AT07",0,0,0,0,
  0,0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_pa_AB00",0,0,0,0,0,
  0,"fire_pa_AT07",0,0,0,0,
  0,0,"fire_pa_GV01",0,0,0,
  0,0,0,"fire_pa_HO00",0,0,
  0,0,0,0,"fire_pa_MC06",0,
  0,0,0,0,0,"fire_pa_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_pa_6m_ppt_AB00",0,0,0,0,0,
  0,"fire_pa_6m_ppt_AT07",0,0,0,0,
  0,0,"fire_pa_6m_ppt_GV01",0,0,0,
  0,0,0,"fire_pa_6m_ppt_HO00",0,0,
  0,0,0,0,"fire_pa_6m_ppt_MC06",0,
  0,0,0,0,0,"fire_pa_6m_ppt_RS02"), 6,18)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_08262022_6state_po4_fire2mpa_fire6mpaxPPT_nolegacies_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_08262022_6state_no3_fire2mpa_fire6mpaxPPT_nolegacies_mBFGS.rds")

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

#           dAIC df
# fit        0.0 36
# null.fit  82.2 18
# RESULT: covar model is better than null, yay!

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall within the CIs? No temporal patterns. Yes, appear to fall within the CIs.

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yup!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. A few fall outside the CIs, but overall look ok.

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Just like NO3, these look great!!

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No patterns.

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00","AT07", "GV01","HO00","MC06","RS02")

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

CIs_fit_ed$Parameter = rep(c("Cum. Ppt", "Fire p/a","Cum. Ppt * Fire p/a (6m)"),6)

CIs_fit_ed$Region = c(rep(c("SB"),6*3)) # **CHECK ORDER OF SITES FOR THIS!!**

# plot results
(RESULTS_ALL_po4 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Phosphate (PO4) MARSS modeling results - 08/26/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

#### MARSS PO4: ppt, % burn (2m win), % burn (6m win) x ppt, no legacy effects ####

# Set up data for MARSS

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_6sites_090722.rds")

# pivot wider for MARSS format
dat_po4 <- dat %>%
  select(
    site, index, 
    mean_po4_uM, 
    cumulative_precip_mm, 
    fire_perc_ws, fire_perc_ws_ppt) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_po4_uM, cumulative_precip_mm, 
                    fire_perc_ws, fire_perc_ws_ppt))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:7)
cov_cols = c(8:25)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 257
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,0,
# 0,"XXXX_AT07",0,0,0,0,
# 0,0,"XXXX_GV01",0,0,0,
# 0,0,0,"XXXX_HO00",0,0,
# 0,0,0,0,"XXXX_MC06",0,
# 0,0,0,0,0,"XXXX_RS02"

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,
  0,"cumulative_precip_mm_AT07",0,0,0,0,
  0,0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_perc_ws_AB00",0,0,0,0,0,
  0,"fire_perc_ws_AT07",0,0,0,0,
  0,0,"fire_perc_ws_GV01",0,0,0,
  0,0,0,"fire_perc_ws_HO00",0,0,
  0,0,0,0,"fire_perc_ws_MC06",0,
  0,0,0,0,0,"fire_perc_ws_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_AB00",0,0,0,0,0,
  0,"fire_perc_ws_ppt_AT07",0,0,0,0,
  0,0,"fire_perc_ws_ppt_GV01",0,0,0,
  0,0,0,"fire_perc_ws_ppt_HO00",0,0,
  0,0,0,0,"fire_perc_ws_ppt_MC06",0,
  0,0,0,0,0,"fire_perc_ws_ppt_RS02"), 6,18)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_09072022_6state_po4_percburn2m_percburn6mxppt_nolegacies_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_09072022_6state_po4_percburn2m_percburn6mxppt_nolegacies_mBFGS.rds")

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

#           dAIC df
# fit        0.0 36
# null.fit  80.2 18
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yes!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# These look AHmazing!

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

# Check HO00 with Alex?

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00", "AT07", "GV01", "HO00", "MC06", "RS02")

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

CIs_fit_ed$Parameter = rep(c("Cum. Ppt", "% Ws Burned (2m)","Cum. Ppt * % Ws Burned (6m)"),6)

CIs_fit_ed$Region = c(rep(c("SB"),6*3))

# plot results
(RESULTS_ALL_po4 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Phosphate (PO4) MARSS modeling results - 09/07/2022\nImmediate Effects") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

# and save out plot
# ggsave(("figures/MARSS_SB_6state_po4_090722.png"),
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### MARSS PO4: ppt, % burn (1y legacy), % burn x ppt 1y legacy ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_6sites_090722.rds")

# pivot wider for MARSS format
dat_po4 <- dat %>%
  select(
    site, index, 
    mean_po4_uM, 
    cumulative_precip_mm, 
    fire_perc_ws_1ylegacy, fire_perc_ws_ppt_1ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_po4_uM, cumulative_precip_mm, 
                    fire_perc_ws_1ylegacy, fire_perc_ws_ppt_1ylegacy))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:7)
cov_cols = c(8:25)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 257
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,0,
# 0,"XXXX_AT07",0,0,0,0,
# 0,0,"XXXX_GV01",0,0,0,
# 0,0,0,"XXXX_HO00",0,0,
# 0,0,0,0,"XXXX_MC06",0,
# 0,0,0,0,0,"XXXX_RS02"

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,0,
  0,"cumulative_precip_mm_AT07",0,0,0,0,
  0,0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_perc_ws_1ylegacy_AB00",0,0,0,0,0,
  0,"fire_perc_ws_1ylegacy_AT07",0,0,0,0,
  0,0,"fire_perc_ws_1ylegacy_GV01",0,0,0,
  0,0,0,"fire_perc_ws_1ylegacy_HO00",0,0,
  0,0,0,0,"fire_perc_ws_1ylegacy_MC06",0,
  0,0,0,0,0,"fire_perc_ws_1ylegacy_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_1ylegacy_AB00",0,0,0,0,0,
  0,"fire_perc_ws_ppt_1ylegacy_AT07",0,0,0,0,
  0,0,"fire_perc_ws_ppt_1ylegacy_GV01",0,0,0,
  0,0,0,"fire_perc_ws_ppt_1ylegacy_HO00",0,0,
  0,0,0,0,"fire_perc_ws_ppt_1ylegacy_MC06",0,
  0,0,0,0,0,"fire_perc_ws_ppt_1ylegacy_RS02"),6,18)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_09072022_6state_po4_percburn1ylegacy_percburnxppt1ylegacy_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_09072022_6state_po4_percburn1ylegacy_percburnxppt1ylegacy_mBFGS.rds")

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

#           dAIC df
# fit        0.0 36
# null.fit 101.6 18
# RESULT: covar model is better than null, thank goodness

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Check with Alex about the clustering happening at GV01/HO00/MC06?

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look excellent!

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns.

# Check HO00/AT07?

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00", "AT07", "GV01", "HO00", "MC06", "RS02")

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

CIs_fit_ed$Parameter = rep(c("Cum. Ppt", "% Ws Burned (1y)","Cum. Ppt * % Ws Burned (1y)"),6)

CIs_fit_ed$Region = c(rep(c("SB"),6*3))

# plot results
(RESULTS_ALL_po4 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Phosphate (PO4) MARSS modeling results - 09/07/2022\n1 Year Legacy") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

# and save out plot
# ggsave(("figures/MARSS_SB_6state_po4_1y_090722.png"),
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### MARSS PO4: ppt, % burn (2y legacy), % burn x ppt 2y legacy ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_6sites_090722.rds")

# NOTE: REMOVE AT07 SITE FROM 2,3,4,&5y LEGACY MODELS - not enough data.

# pivot wider for MARSS format
dat_po4 <- dat %>%
  # remove AT07 site
  filter(site != "AT07") %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    mean_po4_uM, 
    cumulative_precip_mm, 
    fire_perc_ws_2ylegacy, fire_perc_ws_ppt_2ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_po4_uM, cumulative_precip_mm, 
                    fire_perc_ws_2ylegacy, fire_perc_ws_ppt_2ylegacy))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:6)
cov_cols = c(7:21)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 196
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,
# 0,"XXXX_GV01",0,0,0,
# 0,0,"XXXX_HO00",0,0,
# 0,0,0,"XXXX_MC06",0,
# 0,0,0,0,"XXXX_RS02""

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_perc_ws_2ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_2ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_2ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_2ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_2ylegacy_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_2ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_ppt_2ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_ppt_2ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_ppt_2ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_ppt_2ylegacy_RS02"),5,15)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_09072022_5state_po4_percburn2ylegacy_percburnxppt2ylegacy_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_po4_percburn2ylegacy_percburnxppt2ylegacy_mBFGS.rds")

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

#           dAIC df
# fit        0.0 30
# null.fit  96.2 15
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Check with Alex about the clustering happening at GV01/HO00/MC06?

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look fantastic, per usual.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns. HO00 still a lag at 11?

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00", "GV01", "HO00", "MC06", "RS02")

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
                             "% Ws Burned (2y)",
                             "Cum. Ppt * % Ws Burned (2y)"),5)

CIs_fit_ed$Region = c(rep(c("SB"),5*3))

# plot results
(RESULTS_ALL_po4 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Phosphate (PO4) MARSS modeling results - 09/07/2022\n2 Year Legacy") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

# and save out plot
# ggsave(("figures/MARSS_SB_5state_po4_2y_090722.png"),
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### MARSS PO4: ppt, % burn (3y legacy), % burn x ppt 3y legacy ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_6sites_090722.rds")

# NOTE: REMOVE AT07 SITE FROM 2,3,4,&5y LEGACY MODELS - not enough data.

# pivot wider for MARSS format
dat_po4 <- dat %>%
  # remove AT07 site
  filter(site != "AT07") %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    mean_po4_uM, 
    cumulative_precip_mm, 
    fire_perc_ws_3ylegacy, fire_perc_ws_ppt_3ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_po4_uM, cumulative_precip_mm, 
                    fire_perc_ws_3ylegacy, fire_perc_ws_ppt_3ylegacy))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:6)
cov_cols = c(7:21)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 197
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,
# 0,"XXXX_GV01",0,0,0,
# 0,0,"XXXX_HO00",0,0,
# 0,0,0,"XXXX_MC06",0,
# 0,0,0,0,"XXXX_RS02""

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_perc_ws_3ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_3ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_3ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_3ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_3ylegacy_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_3ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_ppt_3ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_ppt_3ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_ppt_3ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_ppt_3ylegacy_RS02"),5,15)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_09072022_5state_po4_percburn3ylegacy_percburnxppt3ylegacy_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_po4_percburn3ylegacy_percburnxppt3ylegacy_mBFGS.rds")

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

#           dAIC df
# fit        0.0 30
# null.fit  91.1 15
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Check with Alex about the clustering happening at GV01/HO00/MC06?

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look fantastic, per usual.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns. HO00 still a lag at 11?

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00", "GV01", "HO00", "MC06", "RS02")

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
                             "% Ws Burned (3y)",
                             "Cum. Ppt * % Ws Burned (3y)"),5)

CIs_fit_ed$Region = c(rep(c("SB"),5*3))

# plot results
(RESULTS_ALL_po4 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Phosphate (PO4) MARSS modeling results - 09/07/2022\n3 Year Legacy") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

# and save out plot
# ggsave(("figures/MARSS_SB_5state_po4_3y_090722.png"),
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### MARSS PO4: ppt, % burn (4y legacy), % burn x ppt 4y legacy ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_6sites_090722.rds")

# NOTE: REMOVE AT07 SITE FROM 2,3,4,&5y LEGACY MODELS - not enough data.

# pivot wider for MARSS format
dat_po4 <- dat %>%
  # remove AT07 site
  filter(site != "AT07") %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    mean_po4_uM, 
    cumulative_precip_mm, 
    fire_perc_ws_4ylegacy, fire_perc_ws_ppt_4ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_po4_uM, cumulative_precip_mm, 
                    fire_perc_ws_4ylegacy, fire_perc_ws_ppt_4ylegacy))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:6)
cov_cols = c(7:21)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 197
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,
# 0,"XXXX_GV01",0,0,0,
# 0,0,"XXXX_HO00",0,0,
# 0,0,0,"XXXX_MC06",0,
# 0,0,0,0,"XXXX_RS02""

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_perc_ws_4ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_4ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_4ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_4ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_4ylegacy_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_4ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_ppt_4ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_ppt_4ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_ppt_4ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_ppt_4ylegacy_RS02"),5,15)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_09072022_5state_po4_percburn4ylegacy_percburnxppt4ylegacy_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_po4_percburn4ylegacy_percburnxppt4ylegacy_mBFGS.rds")

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

#           dAIC df
# fit        0.0 30
# null.fit  83.6 15
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Check with Alex about the clustering happening at GV01/HO00/MC06?

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look fantastic, per usual.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns. HO00 still a lag at 11?

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00", "GV01", "HO00", "MC06", "RS02")

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
                             "% Ws Burned (4y)",
                             "Cum. Ppt * % Ws Burned (4y)"),5)

CIs_fit_ed$Region = c(rep(c("SB"),5*3))

# plot results
(RESULTS_ALL_po4 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Phosphate (PO4) MARSS modeling results - 09/07/2022\n4 Year Legacy") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

# and save out plot
# ggsave(("figures/MARSS_SB_5state_po4_4y_090722.png"),
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### MARSS PO4: ppt, % burn (5y legacy), % burn x ppt 5y legacy ####

# Set up data for MARSS

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects for selected sites
dat = readRDS("data_working/marss_data_sb_6sites_090722.rds")

# NOTE: REMOVE AT07 SITE FROM 2,3,4,&5y LEGACY MODELS - not enough data.

# pivot wider for MARSS format
dat_po4 <- dat %>%
  # remove AT07 site
  filter(site != "AT07") %>%
  # and then continue creating dataset for model
  select(
    site, index, 
    mean_po4_uM, 
    cumulative_precip_mm, 
    fire_perc_ws_5ylegacy, fire_perc_ws_ppt_5ylegacy) %>% 
  pivot_wider(
    names_from = site, 
    values_from = c(mean_po4_uM, cumulative_precip_mm, 
                    fire_perc_ws_5ylegacy, fire_perc_ws_ppt_5ylegacy))

# indicate column #s of response and predictor vars
names(dat_po4)
resp_cols = c(2:6)
cov_cols = c(7:21)

# log and scale transform response var
dat_po4_log = dat_po4
dat_po4_log[,resp_cols] = log10(dat_po4_log[,resp_cols])
dat_po4_log[,resp_cols] = scale(dat_po4_log[,resp_cols])

# check for NaNs (not allowed) and NAs (allowed in response but not predictors)
sum(is.nan(dat_po4_log[,resp_cols])) # 0
sum(is.na(dat_po4_log[,resp_cols])) # 197
range(dat_po4_log[,resp_cols], na.rm = T)

# Pull out only response var
dat_dep <- t(dat_po4_log[,c(resp_cols)])
row.names(dat_dep)

# check covars for nas, nans, or infs (none allowed)
sum(is.nan(dat_po4_log[,cov_cols])) # 0
sum(is.na(dat_po4_log[,cov_cols])) # 0
sum(is.infinite(dat_po4_log[,cov_cols])) # 0

# Make covariate inputs
dat_cov <- dat_po4_log[,c(cov_cols)]
# check for cols with all zeros
any(colSums(dat_cov)==0) # FALSE
# scale and transpose
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)
# check for nas, nans, or infs
sum(is.nan(dat_cov)) # 0
sum(is.na(dat_cov)) # 0
sum(is.infinite(dat_cov)) # 0
# are any rows identical? 
# this can cause model convergence issues
dat_cov[duplicated(dat_cov),]
# no

##### Make C Matrix 

# "XXXX_AB00",0,0,0,0,
# 0,"XXXX_GV01",0,0,0,
# 0,0,"XXXX_HO00",0,0,
# 0,0,0,"XXXX_MC06",0,
# 0,0,0,0,"XXXX_RS02""

CC <- matrix(list( 
  # precip by site: cumulative_precip_mm
  "cumulative_precip_mm_AB00",0,0,0,0,
  0,"cumulative_precip_mm_GV01",0,0,0,
  0,0,"cumulative_precip_mm_HO00",0,0,
  0,0,0,"cumulative_precip_mm_MC06",0,
  0,0,0,0,"cumulative_precip_mm_RS02",
  # fire_pa: fire effect in 2 m window
  "fire_perc_ws_5ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_5ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_5ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_5ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_5ylegacy_RS02",
  # fire_pa_6m_ppt: interaction of cum. ppt with fire p/a in 6 m window
  "fire_perc_ws_ppt_5ylegacy_AB00",0,0,0,0,
  0,"fire_perc_ws_ppt_5ylegacy_GV01",0,0,0,
  0,0,"fire_perc_ws_ppt_5ylegacy_HO00",0,0,
  0,0,0,"fire_perc_ws_ppt_5ylegacy_MC06",0,
  0,0,0,0,"fire_perc_ws_ppt_5ylegacy_RS02"),5,15)

##### Model setup for MARSS 

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

##### Fit MARSS model 

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1, safe=TRUE), fit=TRUE)

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_09072022_5state_po4_percburn5ylegacy_percburnxppt5ylegacy_mBFGS.rds")

##### Diagnoses 

# If you start here, make sure you run the parts of the script above to prepare data for MARSS. It is needed for diagnoses along with the model fit!

# import model fit
fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_po4_percburn5ylegacy_percburnxppt5ylegacy_mBFGS.rds")

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

#           dAIC df
# fit        0.0 30
# null.fit  87.8 15
# RESULT: covar model is better than null

### **** Autoplot diagnoses: VIEW AND RESPOND TO Qs BELOW **** ###
autoplot.marssMLE(fit)

# Plots 1 (xtT) & 2 (fitted.ytT): Do fitted values seem reasonable? Yes

# Plot 3 (model.resids.ytt1): Do resids have temporal patterns? Do 95% of resids fall withing the CIs? No temporal patterns, Yes most fall within CIs.

# Check with Alex about the clustering happening at GV01/HO00/MC06?

# Plot 4 (std.model.resids.ytT): These should all equal zero because we have nothing in the observation model (it is "turned off"). Yep!

# Plot 5 (std.state.resids.xtT): These residuals can be used to detect outliers. Looks ok!

# Plot 6 (qqplot.std.model.resids.ytt1: Are resids normal?
# These are qq plots that should look like a straight line. Datasets with many missing values will not be normal - this isn't a violation per se, but rather you must look at residuals with those associated with missing values removed. 
# Look fantastic, per usual.

# Plot 7 (acf.std.model.resids.ytt1): Do resids have temporal autocorrelation?
# What you don't want is a consistent lag, esp at 1, 6, or 12. Patterns are bad (esp. sinusoidal), random is good. Patterns suggest a seasonal effect is needed.
# No discernible patterns. HO00 still a lag at 11?

### Overall ###
# None of these diagnoses look prohibitively bad.

##### Plot Results 

### Plot coef and coef estimates ###
## estimates
# hessian method is faster but not ideal for final results - should bootstrap for final
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
my_list <- c("AB00", "GV01", "HO00", "MC06", "RS02")

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
                             "% Ws Burned (5y)",
                             "Cum. Ppt * % Ws Burned (5y)"),5)

CIs_fit_ed$Region = c(rep(c("SB"),5*3))

# plot results
(RESULTS_ALL_po4 <- ggplot(CIs_fit_ed, aes(Parameter, Est., color=Region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Phosphate (PO4) MARSS modeling results - 09/07/2022\n5 Year Legacy") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(Region~Site, scales = "free"))

# and save out plot
# ggsave(("figures/MARSS_SB_5state_po4_5y_090722.png"),
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### Summary fig of PO4 results ####

# import model fits
noleg_fit = readRDS(file = "data_working/marss_test_run/fit_09072022_6state_po4_percburn2m_percburn6mxppt_nolegacies_mBFGS.rds")
leg1y_fit = readRDS(file = "data_working/marss_test_run/fit_09072022_6state_po4_percburn1ylegacy_percburnxppt1ylegacy_mBFGS.rds")
leg2y_fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_po4_percburn2ylegacy_percburnxppt2ylegacy_mBFGS.rds") # dropped AT07
leg3y_fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_po4_percburn3ylegacy_percburnxppt3ylegacy_mBFGS.rds")
leg4y_fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_po4_percburn4ylegacy_percburnxppt4ylegacy_mBFGS.rds")
leg5y_fit = readRDS(file = "data_working/marss_test_run/fit_09072022_5state_po4_percburn5ylegacy_percburnxppt5ylegacy_mBFGS.rds")

# extract necessary confidence interval info
noleg_est <- MARSSparamCIs(noleg_fit)
leg1y_est <- MARSSparamCIs(leg1y_fit)
leg2y_est <- MARSSparamCIs(leg2y_fit)
leg3y_est <- MARSSparamCIs(leg3y_fit)
leg4y_est <- MARSSparamCIs(leg4y_fit)
leg5y_est <- MARSSparamCIs(leg5y_fit)

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
CIs = rbind(
  noleg_CI, leg1y_CI, leg2y_CI, leg3y_CI, leg4y_CI, leg5y_CI
)

# add col for site names
CIs$Stream = gsub("_","",str_sub(CIs$Parameter, start= -4))

# simplify parm names
CIs$Parm_simple = c(rep("Ppt",6),
                    rep("Perc. burn",6),
                    rep("Ppt x Perc. burn",6),
                    
                    rep("Ppt",6),
                    rep("Perc. burn",6),
                    rep("Ppt x Perc. burn",6),
                    
                    rep("Ppt",5),
                    rep("Perc. burn",5),
                    rep("Ppt x Perc. burn",5),
                    
                    rep("Ppt",5),
                    rep("Perc. burn",5),
                    rep("Ppt x Perc. burn",5),
                    
                    rep("Ppt",5),
                    rep("Perc. burn",5),
                    rep("Ppt x Perc. burn",5),
                    
                    rep("Ppt",5),
                    rep("Perc. burn",5),
                    rep("Ppt x Perc. burn",5))

# plot results
(RESULTS_ALL_f <- ggplot(CIs, aes(x=factor(Parm_simple, 
                                           levels = c("Ppt x Perc. burn","Perc. burn","Ppt")), 
                                  Est., color=Stream)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0) +
    geom_point(position=position_dodge(width=0.3), size=3) + 
    theme_bw()+
    theme(plot.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20)) +
    geom_hline(aes(yintercept=0), linetype="dashed") +
    coord_flip(ylim = c(-1.25, 1.25)) + 
    labs(y = "", 
         x="",
         title = expression(paste(PO[4]^{"-3"}," MARSS modeling results for Santa Barbara sites"))) +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(.~Model, scales = "free", nrow = 1))

# Export plot.
# ggsave(("MARSS_SB_PO4_091522.png"),
#        path = "figures",
#        width = 80,
#        height = 10,
#        units = "cm"
#        )

# Save these CIs as NO3 specific, so that we can create an enormous facetted plot at the very end.
CIs_PO4 <- CIs
CIs_PO4$Solute <- "PO4"

#### Summary fig of ALL results ####

# NOTE: You must have run the full sections of all 3 previous "Summary fig" code chunks before proceeding with the below.

# Bind all solutes together
CIs_ALL = rbind(
  CIs_NH4, CIs_NO3, CIs_PO4
)

# plot results
(RESULTS_ALL <- ggplot(CIs_ALL, aes(x=factor(Parm_simple, 
                                           levels = c("Ppt x Perc. burn","Perc. burn","Ppt")), 
                                  Est., color=Stream)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.3), width=0) +
    geom_point(position=position_dodge(width=0.3), size=3) + 
    theme_bw()+
    theme(plot.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20)) +
    geom_hline(aes(yintercept=0), linetype="dashed") +
    coord_flip(ylim = c(-0.5, 0.5)) + 
    labs(y = "", 
         x="",
         title = "All MARSS modeling results for Santa Barbara sites") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_grid(Solute~Model))

# Export plot.
# ggsave(("MARSS_SB_NH4_NO3_PO4_091522_zoomin.png"),
#        path = "figures",
#        width = 50,
#        height = 20,
#        units = "cm"
#        )

# End of script.