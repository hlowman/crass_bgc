# Valles Caldera National Preserve (vcnp) Fire Data Assembly
# November 3, 2021
# Betsy Summers following Heili Lowman's workflow

# The following script will assemble the fire datasets available for the VCNP  stream sites into a single data file for use in the MARSS analysis for the CRASS project.

# the main fires of interest that overlap with the grab water quality data set are (date = month of ignition):
# 1) Las Conchas fire (Start date: 26 June, 2011)
# 2) Thompson Ridge fire (start date: 31 May, 2013)


# Load packages
library(plyr) # needs to be loaded prior to dplyr
library(tidyverse) # contains dplyr
library(lubridate)
library(here)


# May not need chem data but was thinking to merge with fire data
chem_reg_vcnp <- read_rds("data_working/VCNPchem_edited_110821.rds")

