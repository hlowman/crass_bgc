# MARSS model @ Santa Barbara Coastal LTER Sites
# October 24, 2021
# Heili Lowman

# The following script will run a MARSS analysis at SBC LTER sites for the CRASS project.
# M - Multi-variate
# AR - Auto-regressive
# S - State
# S - Space

# Note: As of November 3, I've been encountering issues with missing data so this first go-round, trying to make it work with just 2 sites at which I've verified the precip and fire covariate data is complete.

#### Setup ####

# Load packages
library(tidyverse)
library(lubridate)
library(here)
library(MARSS)

# Load datasets
# Stream Chemistry
chem <- readRDS("data_working/SBchem_edited_102421.rds")
# Precipitation
precip <- readRDS("data_working/SBprecip_edited_110321.rds")
# Fire Events
fire <- readRDS("data_working/SBfire_edited_110321.rds")
# Site Location information
location <- read_csv("data_raw/sbc_sites_stream_hydro.csv")

#### Filter & Joining Data ####

# And to join the datasets, I need to identify the matching stream sites for the precip data
precip_ed <- precip %>%
  mutate(sitecode_match = factor(case_when(
    sitecode == "GV202" ~ "GV01",
    sitecode == "HO201" ~ "HO00",
    sitecode == "RG202" ~ "RG01",
    sitecode == "TECA" ~ "TO02",
    sitecode == "GLAN" ~ "BC02",
    sitecode == "SMPA" ~ "SP02",
    sitecode == "GOFS" ~ "DV01",
    sitecode == "GORY" ~ "AT07",
    sitecode == "CAWTP" ~ "AB00",
    sitecode == "SBEB" ~ "MC00",
    sitecode == "ELDE" ~ "RS02"))) %>%
  rename(sitecode_precip = sitecode,
         site_precip = site)

# Joining together with the stream chemistry and wildfire datasets.
# left join with precip so as not to lose any data based on fewer chemistry measurements
precip_chem <- left_join(precip_ed, chem, by = c("Year", "Month", "sitecode_match" = "site_code"))
# and again left join with the larger dataset so as not to accidentally drop records
dat <- left_join(precip_chem, fire, by = c("sitecode_match" = "site_code", "Year", "Month")) %>%
  mutate(Day = 1) %>%
  mutate(Date = make_date(Year, Month, Day)) # And add dates
# Yay! :)

# And, for this first attempt to try and get the MARSS model working, I'm going to filter down to only HO00 and RG01 sites.
dat_2 <- dat %>%
  filter(site %in% c("HO00", "RG01"))

# Now, to check the timeframes of the data available so we can trim down to comparable timespans.
date_check <- dat_2 %>%
  group_by(site) %>%
  summarize(minDate = min(Date), maxDate = max(Date))

# Ok, so need to trim to start date of 9/2002 and end date of 07/2016.
dat_2_trim <- dat_2 %>%
  filter(Date >= "2002-09-01") %>%
  filter(Date <= "2016-07-01")

# And, inspect dataset for missing covariate data.
sum(is.na(dat_2_trim$fire)) # 0
sum(is.na(dat_2_trim$cumulative_precip_mm)) # 0

#### Model fit ####

# Data : Stream Chemistry analytes (NH4, NO3, TDN, TPN, PO4, TDP, TP, TPP, TPC, TSS, SpCond)
# Covariates : Year, Month, Precip, Fire

# Starting with NH4 for test run.

# Note: Not scaling for now, but this should also be added in later.
dat_nh4 <- dat_2_trim %>%
  select(site, Year, Month, mean_nh4_uM, cumulative_precip_mm, fire) %>%
  pivot_wider(names_from = site, values_from = c(mean_nh4_uM, cumulative_precip_mm, fire)) %>%
  #log() %>% # takes the log
  #scale(scale = FALSE) %>% # centers columns of a numeric matrix
  t() # transposes data

# STILL getting NAs?!?
# Ok, so I may need to start with a blank dataframe of the desired date range and then merge everything based on that...
#### STOPPED HERE NOVEMBER 3 ####

# Pull out only NH4 data
dat_nh4_ed <- dat_nh4[3:6,]

# Make covariate inputs
dat_cov <- dat_nh4[c(1:2,7:14),]

# Model setup
mod_list <- list(
  B = "identity",
  U = "zero",
  C = "unconstrained",
  c = dat_cov,
  Q = "diagonal and unequal",
  Z = "identity",
  A = "zero",
  R = "diagonal and equal"
)

# Fit model
fit_nh4 <- MARSS(y = dat_nh4_ed, model = mod_list,
                 control = list(maxit = 5000), method = "BFGS")
#MARSS: NaNs in data are being replaced with NAs.  There might be a problem if NaNs shouldn't be in the data.
#NA is the normal missing value designation.
#Errors were caught in checkModelList 
#model$c is a matrix. No NAs or Infs allowed in this case.
#Error: Stopped in checkModelList() due to specification problem(s).

# End of script.
