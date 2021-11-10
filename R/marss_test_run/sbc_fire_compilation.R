# SB Fire Data Assembly
# November 3, 2021
# Heili Lowman

# The following script will assemble the fire datasets available for the SBC LTER stream sites into a single data file for use in the MARSS analysis for the CRASS project.

# NOTE: This only contains fire info for 2 sites currently, that are being included in the initial MARSS analysis.

# Load packages
library(plyr) # needs to be loaded prior to dplyr
library(tidyverse) # contains dplyr
library(lubridate)
library(here)

# Create base dataframe of dates and sites.
dates1 <- data.frame(seq(as.Date("2002/9/1"), by = "month", length.out = 166)) %>%
  rename(Date = 'seq.as.Date..2002.9.1....by....month...length.out...166.') %>%
  mutate(site = "HO00")

dates2 <- data.frame(seq(as.Date("2002/9/1"), by = "month", length.out = 166)) %>%
  rename(Date = 'seq.as.Date..2002.9.1....by....month...length.out...166.') %>%
  mutate(site = "RG01")
## by month from 9/1/2002 to 7/1/2016

dates <- rbind(dates1, dates2)

# Add watershed data to location dataframe using delineation at:
# https://databasin.org/maps/new/#datasets=6ad26ddb04ae4dbc9362303628270daf
dates_watersheds <- dates %>%
  mutate(watershed = factor(case_when(
    site == "ON02" ~ "Canada De Santa Anita",
    site == "GV01" ~ "Canada De La Gaviota",
    site == "HO00" ~ "Tajiguas Creek",
    site == "RG01" ~ "Tajiguas Creek",
    site == "TO02" ~ "Dos Pueblos Canyon",
    site == "BC02" ~ "Dos Pueblos Canyon",
    site == "DV01" ~ "Dos Pueblos Canyon",
    site == "SP02" ~ "San Pedro Creek",
    site == "AT07" ~ "Atascadero Creek",
    site == "AB00" ~ "Mission Creek",
    site == "MC06" ~ "Mission Creek",
    site == "RS02" ~ "Mission Creek")))

# And here is a list of the wildfire events in the SBC:
# Fire Name (Start Date) - affected watersheds
# Gaviota (6/4/2004) - Canada De Santa Anita, Canada De La Gaviota, Tajiguas Creek
# Sherpa (6/15/2016) - Tajiguas Creek, Dos Pueblos Canyon
# Whittier (7/7/2017) - Tajiguas Creek, Dos Pueblos Canyon
# Gap (7/1/2008) - Dos Pueblos Canyon, San Pedro Creek
# Cave (11/25/2019) - Atascadero Creek
# Jesusita (5/5/2009) - Atascadero Creek, Mission Creek
# Tea (11/13/2008) - Mission Creek
# Thomas (12/4/2017) - Mission Creek

# Now, add in columns for the fires
# 0 - denotes pre-fire months
# 1 - denotes post-fire months
# Each fire will get its own column, so that effects can potentially be additive
# Or so we can combine decay functions later

dates_fire <- dates_watersheds %>%
  mutate(HO00_Gaviota = ifelse(site == "HO00" & Date >= "2004-06-01", 1, 0),
         HO00_Sherpa = ifelse(site == "HO00" & Date >= "2016-06-01", 1, 0),
         HO00_Whittier = ifelse(site == "HO00" & Date >= "2017-07-01", 1, 0),
         RG01_Gaviota = ifelse(site == "RG01" & Date >= "2004-06-01", 1, 0),
         RG01_Sherpa = ifelse(site == "RG01" & Date >= "2016-06-01", 1, 0),
         RG01_Whittier = ifelse(site == "RG01" & Date >= "2017-07-01", 1, 0))

# And export for MARSS script
saveRDS(dates_fire, "data_working/SBfire_edited_111021.rds")

# End of script.
