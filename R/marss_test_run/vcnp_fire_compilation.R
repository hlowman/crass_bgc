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

# NOTE: THIS WAS ALL RE-DONE BY HEILI ON 2/24/22 HAVING FOUND A PREVIOUS ERROR
# WITH HOW THE FIRE DUMMY VARIABLES POPULATED.

# Create base dataframe of dates and sites (from June 2005 to March 2019). This timeframe
# was chosen based on available precipitation data.
# Create sequence of dates
d <- seq(as.Date("2005/6/1"), by = "month", length.out = 166)

# Repeat 7 times for each site
d7 <- rep(d, times = 7)
d7df <- data.frame(d7)

# Create repeated sequence of sites
s <- rep(c("EFJ", "RSAW", "RSA", "IND", "IND_BB", "RED", "SULF"), each=166)
sdf <- data.frame(s)

# Bind dates and sites together
dates <- cbind(d7df, sdf)
colnames(dates) <- c("date","site")

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
# As of 07/28/2022 - we are adding in dummy variables (1s) only for the month of ignition and no decay term
# Conchas = Las Conchas fire
# Thompson = Thompson Ridge fire

# commented out those fires/sites that resulted in a column of all zeros
# which the MARSS model doesn't like
dates_fire_vcnp <- dates %>%
  mutate(RED_Thompson = ifelse(site == "RED" & date == "2013-05-01", 1, 0),
         EFJ_Thompson = ifelse(site == "EFJ" & date == "2013-05-01", 1, 0),
         EFJ_Conchas = ifelse(site == "EFJ" & date == "2011-06-01", 1, 0),
         RSAW_Thompson = ifelse(site == "RSAW" & date == "2013-05-01", 1, 0),
         RSAW_Conchas = ifelse(site == "RSAW" & date == "2011-06-01", 1, 0),
         RSA_Conchas = ifelse(site == "RSA" & date == "2011-06-01", 1, 0),
         IND_Conchas = ifelse(site == "IND" & date == "2011-06-01", 1, 0),
         IND_BB_Conchas = ifelse(site == "IND_BB" & date == "2011-06-01", 1, 0),
         SULF_Thompson = ifelse(site == "SULF" & date == "2013-05-01", 1, 0))
  # # also adding in 4 year decay term for each of the fires
  # mutate(RED_Thompson_d = NA,
  #        EFJ_Thompson_d = NA,
  #        EFJ_Conchas_d = NA,
  #        RSA_Conchas_d = NA,
  #        RSAW_Thompson_d = NA,
  #        RSAW_Conchas_d = NA,
  #        IND_Conchas_d = NA,
  #        IND_BB_Conchas_d = NA,
  #        SULF_Thompson_d = NA)

# values <- rev(seq(1, 48, by = 1))
# decay <- exp(values)
# 
# dates_fire_vcnp$RED_Thompson_d[927:974] <- decay
# dates_fire_vcnp$EFJ_Thompson_d[97:144] <- decay
# dates_fire_vcnp$EFJ_Conchas_d[74:121] <- decay
# dates_fire_vcnp$RSA_Conchas_d[406:453] <- decay
# dates_fire_vcnp$RSAW_Thompson_d[263:310] <- decay
# dates_fire_vcnp$RSAW_Conchas_d[240:287] <- decay
# dates_fire_vcnp$IND_Conchas_d[572:619] <- decay
# dates_fire_vcnp$IND_BB_Conchas_d[738:785] <- decay
# dates_fire_vcnp$SULF_Thompson_d[1093:1140] <- decay

dates_fire_vcnp[is.na(dates_fire_vcnp)] = 0

# And export for MARSS script
#saveRDS(dates_fire_vcnp, "data_working/VCNPfire_edited_060622.rds")
saveRDS(dates_fire_vcnp, "data_working/VCNPfire_edited_072822.rds")

# End of script.
