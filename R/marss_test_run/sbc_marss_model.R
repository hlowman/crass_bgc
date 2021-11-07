# MARSS model @ Santa Barbara Coastal LTER Sites
# October 24, 2021
# Heili Lowman

# The following script will run a MARSS analysis at SBC LTER sites for the CRASS project.
# M - Multi-variate
# AR - Auto-regressive
# S - State
# S - Space

# Note: As of November 3, I've been encountering issues with missing data so this first go-round, trying to make it work with just a few sites at which I've verified the precip and fire covariate data is complete.

#### Setup ####

# Load packages
library(tidyverse)
library(lubridate)
library(here)
library(MARSS)

# Load datasets
# Stream Chemistry - all sites
chem <- readRDS("data_working/SBchem_edited_110721.rds")
# Precipitation - all sites
precip <- readRDS("data_working/SBprecip_edited_110721.rds")
# Fire Events - HO00 and RG01 sites only so far
fire <- readRDS("data_working/SBfire_edited_110721.rds")
# Site Location information
location <- read_csv("data_raw/sbc_sites_stream_hydro.csv")

#### Filter & Joining Data ####

# First, to avoid strange gaps in data, I'm creating a dummy column with the dates of interest.
# And I need to do this for each site, since otherwise it'll allow gaps.
# Just being really careful since this was causing issues previously.
dates1 <- data.frame(seq(as.Date("2002/9/1"), by = "month", length.out = 166)) %>%
  rename(Date = 'seq.as.Date..2002.9.1....by....month...length.out...166.') %>%
  mutate(site = "HO00")

dates2 <- data.frame(seq(as.Date("2002/9/1"), by = "month", length.out = 166)) %>%
  rename(Date = 'seq.as.Date..2002.9.1....by....month...length.out...166.') %>%
  mutate(site = "RG01")
## by month from 9/1/2002 to 7/1/2016

dates <- rbind(dates1, dates2)

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
         site_precip = site) %>%
  mutate(Day = 1) %>% # new column of "days"
  mutate(Date = make_date(Year, Month, Day))

# First, join precip with the dates.
dates_precip <- left_join(dates, precip_ed, by = c("Date", "site" = "sitecode_match"))
# Then, left join precip with chemistry so as not to lose any data.
precip_chem <- left_join(dates_precip, chem, by = c("Year", "Month", "site"))
# And again left join with fire so as not to accidentally drop data.
dat <- left_join(precip_chem, fire, by = c("Date", "site"))

# And, for this first attempt to try and get the MARSS model working, I'm going to filter down to only HO00 and RG01 sites.
dat_2 <- dat %>%
  filter(site %in% c("HO00", "RG01"))

# The trimming below is no longer necessary since we created a primary column with this same date range at the beginning of this script.
# Now, to check the timeframes of the data available so we can trim down to comparable timespans.
# date_check <- dat_2 %>%
#   group_by(site) %>%
#   summarize(minDate = min(Date), maxDate = max(Date))

# Ok, so need to trim to start date of 9/2002 and end date of 07/2016.
# dat_2_trim <- dat_2 %>%
#   filter(Date >= "2002-09-01") %>%
#   filter(Date <= "2016-07-01")

# And, inspect dataset for missing covariate data.
sum(is.na(dat_2$fire)) # 0
sum(is.na(dat_2$cumulative_precip_mm)) # 0
# Great!

# Note for future me - be VERY careful with the joining above. Something weird was happening previously where precip data that IS present was simply dropping off.

#### Model fit ####

# Data : Stream Chemistry analytes (NH4, NO3, TDN, TPN, PO4, TDP, TP, TPP, TPC, TSS, SpCond)
# Covariates : Year, Month, Precip, Fire

# Starting with NH4 for test run.

# Note: Not scaling for now, but this should also be added in later.
dat_nh4 <- dat_2 %>%
  select(site, Year, Month, mean_nh4_uM, cumulative_precip_mm, fire) %>%
  pivot_wider(names_from = site, values_from = c(mean_nh4_uM, cumulative_precip_mm, fire)) %>%
  #log() %>% # takes the log
  #scale(scale = FALSE) %>% # centers columns of a numeric matrix
  t() # transposes data

# Pull out only NH4 data
dat_dep <- dat_nh4[3:4,]

# Make covariate inputs
dat_cov <- dat_nh4[c(1:2,5:8),]

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
fit <- MARSS(y = dat_dep, model = mod_list,
                 control = list(maxit = 5000), method = "BFGS")

# Still getting the following message:
#MARSS: NaNs in data are being replaced with NAs.  There might be a problem if NaNs shouldn't be in the data.
#NA is the normal missing value designation.

# It worked!!!

# End of script.
