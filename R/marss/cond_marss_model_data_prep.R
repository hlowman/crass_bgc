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
  filter(index > 26)

# Now, both datasets are the same size (724 observations).

# Join and export to save progress.
dat_select <- rbind(dat_select_sb, dat_select_vc_181)
saveRDS(dat_select, "data_working/marss_data_sb_vc_nolegacies_081123.rds")

# End of script.
