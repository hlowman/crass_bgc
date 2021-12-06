# Valles Caldera National Preserve (vcnp) Fire Data Assembly
# November 3, 2021
# Betsy Summers following Heili Lowman's workflow

# The following script will assemble the fire datasets available for the VCNP  stream sites into a single data file for use in the MARSS analysis for the CRASS project.

# the main fires of interest that overlap with the grab water quality data set are (date = month of ignition):
# 1) Las Conchas fire (Start date: 26 June, 2011)
# 2) Thompson Ridge fire (start date: 31 May, 2013)


## Load packages
library(plyr) # needs to be loaded prior to dplyr
library(tidyverse) # contains dplyr
library(lubridate)
library(here)


## May not need chem data but was thinking to merge with fire data
#chem_reg_vcnp <- read_rds("data_working/VCNPchem_edited_110821.rds")
chem_reg_vcnp <- read_rds("data_working/VCNPchem_edited_120521.rds")

str(chem_reg_vcnp)
#list(chem_reg_vcnp$site_code)

chem_reg_vcnp <- chem_reg_vcnp %>% rename(site = "site_code") # rename to site to match sbc variable name

## Tidy site_code variable list:
# 1) Remove rows with Nas
chem_reg_vcnp <- chem_reg_vcnp %>% 
  filter(!is.na(site)) %>%
  filter(!is.na(Date)) %>%
  filter(site != "La Jara Well - HQ Drinking Water, Pre-treatment") # remove random site, not stream

# 2) Need to fix variable names because inconsistent naming (location of dashes)



## Fire names and start dates of burn - sites affected
# Las Conchas (2011-06-26) - EFJ, RSAW, RSA, IND, IND_BB
# Thompson Ridge (2013-05-31) - RED, EFJ, RSAW, SULF
# Prescribed burn (2016-05-11) - EFJ [Out of time frame]

## Add in columns for the fires based on start dates:
# 0 - denotes pre-fire months
# 1 - denotes post-fire months
# Each fire will get its own column, so that effects can potentially be additive
# Or so we can combine decay functions later

# currently use Date variable to determine when site is 1 or 0.
# As of 12/6/2021 - we are adding in dummy variable (1s) only for the year
# following the fire event
# Conchas = LasConchas
# Thompson = Thompson Ridge

namelist <- c("San Antonio Creek - Toledo", "San Antonio Creek- Toledo", "San Antonio Creek -Toledo")

# commented out those fires/sites that resulted in a column of all zeros
# which the MARSS model doesn't like
dates_fire_vcnp <- chem_reg_vcnp %>%
  mutate(RED_Thompson = ifelse(site == "Redondo Creek" & Date >= "2013-05-31" & Date < "2014-05-31", 1, 0),
         EFJ_Thompson = ifelse(site == "East Fork Jemez River" & Date >= "2013-05-31" & Date < "2014-05-31", 1, 0),
         #EFJ_Prescribe = ifelse(site == "East Fork Jemez River" & Date >= "2016-05-11" & Date < "2017-05-11", 1, 0),
         EFJ_Conchas = ifelse(site == "East Fork Jemez River" & Date >= "2011-06-26" & Date < "2012-06-26", 1, 0),
         RSAW_Thompson = ifelse(site == "San Antonio - West" & Date >= "2013-05-31" & Date < "2014-05-31", 1, 0),
         RSAW_Conchas = ifelse(site == "San Antonio - West" & Date >= "2011-06-26" & Date < "2012-06-26", 1, 0),
         RSA_Conchas = ifelse(site %in% namelist & Date >= "2011-06-26" & Date < "2012-06-26", 1, 0),
         #IND_Conchas = ifelse(site == "Indios Creek" & Date >= "2011-06-26" & Date < "2012-06-26", 1, 0),
         IN_BB_Conchas = ifelse(site == "Indios Creek - Post Fire (Below Burn)" & Date >= "2011-06-26" & Date < "2012-06-26", 1, 0),
         SULF_Thompson = ifelse(site == "Sulfur Creek" & Date >= "2013-05-31" & Date < "2014-05-31", 1, 0))

## Need information on what sites were within the burn area or downstream of burn area to better designate if a site gets a 1 or 0 after fire start date. 
# Sites not impacted by wildfire just receive a 0. 

# And export for MARSS script
saveRDS(dates_fire_vcnp, "data_working/VCNPfire_edited_120621.rds")

# End of script.
