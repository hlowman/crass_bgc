# Nutrient MARSS models
# with fire x ppt interactions and legacy effects
# as well as 4 state structure for CA sites
# Script started August 8, 2023 by Heili Lowman

# This script will prep data for running in the nutrient MARSS models.

#### Setup ####

# Load packages.
library(tidyverse)
library(lubridate)
library(MARSS)
library(naniar) 

# Load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))

# Load data.
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

#### Compile Data ####

# I first need to identify the matching stream sites for the precip data
precip_ed <- precip %>%
  mutate(sitecode_match = factor(case_when(
    sitecode == "GV202" ~ "GV01",
    sitecode == "HO202" ~ "HO00", # switched back from "BARA" since Q data doesn't
    # start until 2002 either, using HO202 instead of 201 because record is longer
    sitecode == "CAWTP" ~ "AB00",
    sitecode == "ELDE" ~ "RS02"))) %>%
  dplyr::rename(cumulative_precip_mm = c_precip_mm,
                sitecode_precip = sitecode) %>%
  mutate(Day = 1) %>% # new column of "days"
  mutate(Date = make_date(Year, Month, Day))

# Visually double-check edited fire data.
sites <- c("AB00", "GV01", "HO00", "RS02")

ggplot(fire %>% filter(site %in% sites), 
       aes(x = date, y = fire_perc_ws)) +
  geom_point() +
  labs(x = "Date") +
  facet_wrap(.~site, scales = "free",
             ncol = 1) +
  theme_bw() # looks good based on fire dates in "sbc_fire_compilation.R" script

# Examine precip data coverage for MARSS timescale delineation.
precip_ed %>%
  drop_na(sitecode_match) %>%
  ggplot(aes(x = Year, y = sitecode_match, color = sitecode_match)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") # 2002-2017

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
  filter(site %in% c("AB00", "GV01", "HO00", "RS02"))

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

# No longer adding seasonal dummy covariates.
# But adding index values in for future use.
dat_full <- dat_full %>%
  mutate(index = rep(seq(1,182), 4))

# Replace NaNs with NAs
dat_full[is.nan(dat_full)] <- NA

# And, inspect dataset for missing covariate data.
sum(is.na(dat_full$fire_perc_ws)) # 0
sum(is.na(dat_full$cumulative_precip_mm)) # 2
# Need to trim out October 2017.

dat_full <- dat_full %>%
  filter(date != as_date(as.character("2017-10-01")))

# Trim down to columns of interest
dat_select <- dat_full %>%
  select(year, month, index, site, ws_area_m2,
         cumulative_precip_mm, 
         fire_ID, ig_date, fire_pa, ws_fire_area_m2, fire_perc_ws,
         vwm_nh4, vwm_no3, vwm_po4) %>%
  mutate(region = "SB")

#### Add fire legacy effects ####

# Reorganize
dat2 <- dat_select[,c(1,2,3,15,4,5, #"year" "month" "index" "region" "site" "ws_area_m2"
               6, # cumulative_precip_mm
               7:11, # fire info
               12:14)] # analytes

# Create df of sites x fire info
firez <- unique(dat_select[,c(4,7,8,10,11)])
firez <- firez[complete.cases(firez),] # removes NAs
firez$ig_date <- gsub("/","-",firez$ig_date) # reformat date
firez$year <- year(as.Date(firez$ig_date)) # create year column
firez$month <- month(as.Date(firez$ig_date)) # create month column
firez$date <- as.Date(paste(firez$year, firez$month, "01", sep="-")) # reformat date again

# Create new effect date column.
firez$effect_date <- firez$date

# 0.5 year/6mo legacy (this is to allow a window for the fire x ppt interaction only)
firedates_0.5ylegacy <- rbind(firez, 
                             cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(1)),
                             cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(2)),
                             cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(3)),
                             cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(4)),
                             cbind(firez[,1:8], "effect_date" = firez$effect_date %m+% months(5)))

firedates_0.5ylegacy <- data.frame(
  year = year(firedates_0.5ylegacy$effect_date),
  month = month(firedates_0.5ylegacy$effect_date),
  site = firedates_0.5ylegacy$site,
  fire_pa_0.5ylegacy = 1,
  ws_fire_area_m2_0.5ylegacy = firedates_0.5ylegacy$ws_fire_area_m2,
  fire_perc_ws_0.5ylegacy = firedates_0.5ylegacy$fire_perc_ws)

firedates_0.5ylegacy <- firedates_0.5ylegacy %>% 
  group_by(year, month, site, fire_pa_0.5ylegacy) %>% 
  dplyr::summarise(ws_fire_area_m2_0.5ylegacy = sum(ws_fire_area_m2_0.5ylegacy),
                   fire_perc_ws_0.5ylegacy = sum(fire_perc_ws_0.5ylegacy)) %>%
  ungroup()

# 1 year legacy
firedates_1ylegacy <- rbind(firez, 
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

firedates_1ylegacy <- data.frame(year = year(firedates_1ylegacy$effect_date),
                                month = month(firedates_1ylegacy$effect_date),
                                site = firedates_1ylegacy$site,
                                fire_pa_1ylegacy = 1,
                                ws_fire_area_m2_1ylegacy = firedates_1ylegacy$ws_fire_area_m2,
                                fire_perc_ws_1ylegacy = firedates_1ylegacy$fire_perc_ws)

firedates_1ylegacy <- firedates_1ylegacy %>% 
  group_by(year, month, site, fire_pa_1ylegacy) %>% 
  dplyr::summarise(ws_fire_area_m2_1ylegacy = sum(ws_fire_area_m2_1ylegacy),
                   fire_perc_ws_1ylegacy = sum(fire_perc_ws_1ylegacy)) %>%
  ungroup()

# 2 year legacy
firedates_2ylegacy <- rbind(firez, 
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

firedates_2ylegacy <- data.frame(year = year(firedates_2ylegacy$effect_date),
                                month = month(firedates_2ylegacy$effect_date),
                                site = firedates_2ylegacy$site,
                                fire_pa_2ylegacy = 1,
                                ws_fire_area_m2_2ylegacy = firedates_2ylegacy$ws_fire_area_m2,
                                fire_perc_ws_2ylegacy = firedates_2ylegacy$fire_perc_ws)

firedates_2ylegacy <- firedates_2ylegacy %>% 
  group_by(year, month, site, fire_pa_2ylegacy) %>% 
  dplyr::summarise(ws_fire_area_m2_2ylegacy = sum(ws_fire_area_m2_2ylegacy),
                   fire_perc_ws_2ylegacy = sum(fire_perc_ws_2ylegacy)) %>%
  ungroup()

# 3 year legacy
firedates_3ylegacy <- rbind(firez, 
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

firedates_3ylegacy <- data.frame(year = year(firedates_3ylegacy$effect_date),
                                month = month(firedates_3ylegacy$effect_date),
                                site = firedates_3ylegacy$site,
                                fire_pa_3ylegacy = 1,
                                ws_fire_area_m2_3ylegacy = firedates_3ylegacy$ws_fire_area_m2,
                                fire_perc_ws_3ylegacy = firedates_3ylegacy$fire_perc_ws)

firedates_3ylegacy <- firedates_3ylegacy %>% 
  group_by(year, month, site, fire_pa_3ylegacy) %>% 
  dplyr::summarise(ws_fire_area_m2_3ylegacy = sum(ws_fire_area_m2_3ylegacy),
                   fire_perc_ws_3ylegacy = sum(fire_perc_ws_3ylegacy)) %>%
  ungroup()

# 4 year legacy
firedates_4ylegacy <- rbind(firez, 
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

firedates_4ylegacy <- data.frame(year = year(firedates_4ylegacy$effect_date),
                                month = month(firedates_4ylegacy$effect_date),
                                site = firedates_4ylegacy$site,
                                fire_pa_4ylegacy = 1,
                                ws_fire_area_m2_4ylegacy = firedates_4ylegacy$ws_fire_area_m2,
                                fire_perc_ws_4ylegacy = firedates_4ylegacy$fire_perc_ws)

firedates_4ylegacy <- firedates_4ylegacy %>% 
  group_by(year, month, site, fire_pa_4ylegacy) %>% 
  dplyr::summarise(ws_fire_area_m2_4ylegacy = sum(ws_fire_area_m2_4ylegacy),
                   fire_perc_ws_4ylegacy = sum(fire_perc_ws_4ylegacy)) %>%
  ungroup()

# 5 year legacy
firedates_5ylegacy <- rbind(firez, 
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

firedates_5ylegacy <- data.frame(year = year(firedates_5ylegacy$effect_date),
                                month = month(firedates_5ylegacy$effect_date),
                                site = firedates_5ylegacy$site,
                                fire_pa_5ylegacy = 1,
                                ws_fire_area_m2_5ylegacy = firedates_5ylegacy$ws_fire_area_m2,
                                fire_perc_ws_5ylegacy = firedates_5ylegacy$fire_perc_ws)

firedates_5ylegacy <- firedates_5ylegacy %>% 
  group_by(year, month, site, fire_pa_5ylegacy) %>% 
  dplyr::summarise(ws_fire_area_m2_5ylegacy = sum(ws_fire_area_m2_5ylegacy),
                   fire_perc_ws_5ylegacy = sum(fire_perc_ws_5ylegacy)) %>%
  ungroup()

# join legacy dates to data
dat3 <- left_join(dat2, firedates_0.5ylegacy, by=c("year","month","site"))
dat4 <- left_join(dat3, firedates_1ylegacy, by=c("year","month","site"))
dat5 <- left_join(dat4, firedates_2ylegacy, by=c("year","month","site"))
dat6 <- left_join(dat5, firedates_3ylegacy, by=c("year","month","site"))
dat7 <- left_join(dat6, firedates_4ylegacy, by=c("year","month","site"))
dat8 <- left_join(dat7, firedates_5ylegacy, by=c("year","month","site"))

# replace NAs in fire columns with 0s
dat8[,16:33][is.na(dat8[,16:33])] = 0

#### Add fire x precip legacy effects ####

# for fire perc burn ws
dat8$fire_perc_ws_ppt <- dat8$fire_perc_ws_0.5ylegacy*dat8$cumulative_precip_mm
dat8$fire_perc_ws_ppt_1ylegacy <- dat8$fire_perc_ws_1ylegacy*dat8$cumulative_precip_mm
dat8$fire_perc_ws_ppt_2ylegacy <- dat8$fire_perc_ws_2ylegacy*dat8$cumulative_precip_mm
dat8$fire_perc_ws_ppt_3ylegacy <- dat8$fire_perc_ws_3ylegacy*dat8$cumulative_precip_mm
dat8$fire_perc_ws_ppt_4ylegacy <- dat8$fire_perc_ws_4ylegacy*dat8$cumulative_precip_mm
dat8$fire_perc_ws_ppt_5ylegacy <- dat8$fire_perc_ws_5ylegacy*dat8$cumulative_precip_mm

# Quick plot of 1yr legacy effects to be sure things are populating correctly.
qplot(index, fire_perc_ws_ppt_1ylegacy, data=dat8, colour=site, geom="path", facets = "region")

#### Examine data closely ####

# precip (0-1000 mm)
ggplot(dat8, aes(x = index, y = cumulative_precip_mm, color = site)) +
  geom_point() +
  facet_wrap(.~site) +
  theme(legend.position = "none")

# NH4 (0-150 uM)
ggplot(dat8, aes(x = index, y = vwm_nh4, color = site)) +
  geom_point() +
  facet_wrap(.~site) +
  theme(legend.position = "none")

# NO3 (0-600 uM)
ggplot(dat8, aes(x = index, y = vwm_no3, color = site)) +
  geom_point() +
  facet_wrap(.~site) +
  theme(legend.position = "none")

# PO4 (0-30 uM)
ggplot(dat8, aes(x = index, y = vwm_po4, color = site)) +
  geom_point() +
  facet_wrap(.~site) +
  theme(legend.position = "none")

#### Remove outliers ####

# Examined chem_reg file from the sbc_chem_compilation.R script to compare weekly solute values to
# measured monthly means to make sure the values made sense and they appear to.
chem_reg <- read_csv("data_raw/sbclter_stream_chemistry_allyears_registered_stations_20190628.csv")

# NH4 
ggplot(chem_reg %>%
         filter(nh4_uM > 0) %>%
         filter(site_code %in% c("AB00", "GV01", "HO00", "RS02")), 
       aes(x = timestamp_local, y = nh4_uM, color = site_code)) +
  geom_point() +
  facet_wrap(.~site_code) +
  theme(legend.position = "none")

# NO3
ggplot(chem_reg %>%
         filter(no3_uM > 0) %>%
         filter(site_code %in% c("AB00", "GV01", "HO00", "RS02")), 
       aes(x = timestamp_local, y = no3_uM, color = site_code)) +
  geom_point() +
  facet_wrap(.~site_code) +
  theme(legend.position = "none")

# PO4 
ggplot(chem_reg %>%
         filter(po4_uM > 0) %>%
         filter(site_code %in% c("AB00", "GV01", "HO00", "RS02")), 
       aes(x = timestamp_local, y = po4_uM, color = site_code)) +
  geom_point() +
  facet_wrap(.~site_code) +
  theme(legend.position = "none")

# Raw, unaveraged, and unweighted values for the watersheds in 
# question ranged from:
# 0 - 545 uM for NH4 (GV01 has max value)
# 0 - 1763 uM for NO3 (RS02 has max value)
# 0 - 59 uM for PO4 (GV01 has max value)

# Keeping all data points for now.

#### Export data ready for modeling ####

saveRDS(dat8, "data_working/marss_data_sb_080823.rds")

# End of script.
