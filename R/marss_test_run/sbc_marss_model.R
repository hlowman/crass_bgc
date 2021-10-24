# MARSS model @ Santa Barbara Coastal LTER Sites
# October 24, 2021
# Heili Lowman

# The following script will run a MARSS analysis at SBC LTER sites for the CRASS project.

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

# End of script.

