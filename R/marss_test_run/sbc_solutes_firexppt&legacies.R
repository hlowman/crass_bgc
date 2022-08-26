# SB MARSS models for other solutes
# Script started August 26, 2022
# Heili Lowman, Alex Webster

#### READ ME ####

# The following script will prepare data and run MARSS analyses at SBC sites for the CRASS project.
# MARSS modeling sections each contain 1 model or 1 model comparision and restart the script and import data anew each time.
# Be sure to load libraries before starting modeling sections. 
# Much of this code has been copied over from the "sbc_vcnp_firexppt&legacies.R" script written by Alex, but removing any of the NM site processing/modeling.

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
# Stream Chemistry - all sites
# see script "sbc_chem_compilation.R" for the code used to tidy/generate this dataset
chem <- readRDS("data_working/SBchem_edited_120321.rds")

# Precipitation - all sites
precip <- readRDS("data_working/SBprecip_edited_120121.rds")

# Fire Events - all sites
fire <- readRDS("data_working/SBfire_edited_072922.rds")

# Site Location information
location <- read_csv("data_raw/sbc_sites_stream_hydro.csv")

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
sitez = c("AB00", "GV01", "MC06", "RG01", "RS02", "HO00")

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

# Trim down to columns of interest
dat_select <- dat %>%
  select(year, month, index, site, ws_area_m2,
         cumulative_precip_mm, 
         fire_ID, ig_date, fire_pa, ws_fire_area_m2, fire_perc_ws,
         mean_nh4_uM, mean_no3_uM, mean_po4_uM, 
         Season1, Season2) %>%
  mutate(region = "SB")

# And export to save progress
saveRDS(dat_select, "data_working/marss_data_sb_N_P_082622.rds")
