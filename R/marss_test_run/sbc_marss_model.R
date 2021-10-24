# MARSS model @ Santa Barbara Coastal LTER Sites
# October 24, 2021
# Heili Lowman

# The following script will run a MARSS analysis at SBC LTER sites for the CRASS project.
# M - Multi-variate
# AR - Auto-regressive
# S - State
# S - Space

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
precip <- readRDS("data_working/SBprecip_edited_102421.rds")
# Site Location information
location <- read_csv("data_raw/sbc_sites_stream_hydro.csv")

#### Joining Data ####

# Fire dataset/parameters may change, so including it in the MARSS script for this iteration.

# Add watershed data to location dataframe using delineation at:
# https://databasin.org/maps/new/#datasets=6ad26ddb04ae4dbc9362303628270daf
location_streams <- location %>%
  filter(Data == "Stream Chemistry") %>%
  mutate(watershed = factor(case_when(
    sitecode == "ON02" ~ "Canada De Santa Anita",
    sitecode == "GV01" ~ "Canada De La Gaviota",
    sitecode == "HO00" ~ "Tajiguas Creek",
    sitecode == "RG01" ~ "Tajiguas Creek",
    sitecode == "TO02" ~ "Dos Pueblos Canyon",
    sitecode == "BC02" ~ "Dos Pueblos Canyon",
    sitecode == "DV01" ~ "Dos Pueblos Canyon",
    sitecode == "SP02" ~ "San Pedro Creek",
    sitecode == "AT07" ~ "Atascadero Creek",
    sitecode == "AB00" ~ "Mission Creek",
    sitecode == "MC00" ~ "Mission Creek",
    sitecode == "RS02" ~ "Mission Creek")))

# Joining watershed info to larger chemistry dataset
chem_location <- left_join(chem, location_streams, by = c("site_code" = "sitecode"))

# And now to add a fire column based on the following wildfire events:
# Fire Name (Start Date) - affected watersheds
# Gaviota (6/4/2004) - Canada De Santa Anita, Canada De La Gaviota, Tajiguas Creek
# Sherpa (6/15/2016) - Tajiguas Creek, Dos Pueblos Canyon
# Whittier (7/7/2017) - Tajiguas Creek, Dos Pueblos Canyon
# Gap (7/1/2008) - Dos Pueblos Canyon, San Pedro Creek
# Cave (11/25/2019) - Atascadero Creek
# Jesusita (5/5/2009) - Atascadero Creek, Mission Creek
# Tea (11/13/2008) - Mission Creek
# Thomas (12/4/2017) - Mission Creek

chem_fire <- chem_location %>%
  mutate(fire = case_when(Year == 2004 & Month == 6 & watershed == "Canada De Santa Anita" ~ 1,
                          Year == 2004 & Month == 6 & watershed == "Canada De La Gaviota" ~ 1,
                          Year == 2004 & Month == 6 & watershed == "Tajiguas Creek" ~ 1,
                          Year == 2016 & Month == 6 & watershed == "Tajiguas Creek" ~ 1,
                          Year == 2016 & Month == 6 & watershed == "Dos Pueblos Canyon" ~ 1,
                          Year == 2017 & Month == 7 & watershed == "Tajiguas Creek" ~ 1,
                          Year == 2017 & Month == 7 & watershed == "Dos Pueblos Canyon" ~ 1,
                          Year == 2008 & Month == 7 & watershed == "Dos Pueblos Canyon" ~ 1,
                          Year == 2008 & Month == 7 & watershed == "San Pedro Creek" ~ 1,
                          Year == 2019 & Month == 11 & watershed == "Atascadero Creek" ~ 1,
                          Year == 2009 & Month == 5 & watershed == "Atascadero Creek" ~ 1,
                          Year == 2009 & Month == 5 & watershed == "Mission Creek" ~ 1,
                          Year == 2008 & Month == 11 & watershed == "Mission Creek" ~ 1,
                          Year == 2017 & Month == 12 & watershed == "Mission Creek" ~ 1,
                          TRUE ~ 0)) %>%
  filter(site_code != "MC06") # And remove the USGS site

# And finally, to join the datasets, I need to identify the matching stream sites for the precip data
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
dat <- left_join(chem_fire, precip_ed, by = c("Year", "Month", "site_code" = "sitecode_match"))
# Yay! :)
# Note, Tecolotito Creek disappears because available chemistry and precipitation datasets
# do not match in terms of date.

#### Model fit ####

# Data : Stream Chemistry analytes (NH4, NO3, TDN, TPN, PO4, TDP, TP, TPP, TPC, TSS, SpCond)
# Covariates : Year, Month, Watershed, Precip, Fire

# Starting with NH4 for test run.

# Note: Not scaling for now, but this should also be added in at a later date.
dat_nh4 <- dat %>%
  select(Year:mean_nh4_uM, site, watershed, fire, cumulative_precip_mm) %>%
  pivot_wider(names_from = c(site, watershed), values_from = c(mean_nh4_uM, fire, cumulative_precip_mm)) %>%
  #log() %>% # takes the log
  #scale(scale = FALSE) %>% # centers columns of a numeric matrix
  t() # transposes data

# Pull out only NH4 data
dat_nh4_ed <- dat_nh4[3:14,]

# Make covariate inputs
dat_cov <- dat_nh4[c(1:2,15:38),]

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
