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

# Load precipitation dataset from which dates will be pulled.
precip <- readRDS("data_working/SBprecip_edited_110321.rds")

# Filter for 2 sites used in initial MARSS analysis
precip_filter <- precip %>%
  filter(sitecode %in% c("HO201", "RG202")) %>%
  mutate(sitecode_match = factor(case_when(
    sitecode == "HO201" ~ "HO00",
    sitecode == "RG202" ~ "RG01")))

# And pull out only the site and date columns
dates_only <- precip_filter %>%
  select(sitecode_match, Year, Month) %>%
  mutate(Day = 1) %>% # adding a day to create a full date
  mutate(Date = make_date(Year, Month, Day))

# Add watershed data to location dataframe using delineation at:
# https://databasin.org/maps/new/#datasets=6ad26ddb04ae4dbc9362303628270daf
dates_watersheds <- dates_only %>%
  mutate(watershed = factor(case_when(
    sitecode_match == "ON02" ~ "Canada De Santa Anita",
    sitecode_match == "GV01" ~ "Canada De La Gaviota",
    sitecode_match == "HO00" ~ "Tajiguas Creek",
    sitecode_match == "RG01" ~ "Tajiguas Creek",
    sitecode_match == "TO02" ~ "Dos Pueblos Canyon",
    sitecode_match == "BC02" ~ "Dos Pueblos Canyon",
    sitecode_match == "DV01" ~ "Dos Pueblos Canyon",
    sitecode_match == "SP02" ~ "San Pedro Creek",
    sitecode_match == "AT07" ~ "Atascadero Creek",
    sitecode_match == "AB00" ~ "Mission Creek",
    sitecode_match == "MC06" ~ "Mission Creek",
    sitecode_match == "RS02" ~ "Mission Creek"))) %>%
  rename(site_code = sitecode_match)

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

# Now, add in a column for the fire variable
# 0 - denotes pre-fire months
# 1 - denotes post-fire months, beginning during the month of ignition
# 2+ - refers to additional post-fire event months

dates_fire <- dates_watersheds %>%
  mutate(fire = case_when(
    # Arroyo Hondo - Tajiguas Creek Watershed
    site_code == "HO00" & Date >= "2016-06-01" & Date < "2017-07-01" ~ 1, # Sherpa Fire
    site_code == "HO00" & Date >= "2017-07-01" ~ 2, # Whittier Fire
    # Refugio - also Tajiguas Creek Watershed
    site_code == "RG01" & Date >= "2016-06-01" & Date < "2017-07-01"~ 1, # Sherpa Fire
    site_code == "RG01" & Date >= "2017-07-01" ~ 2, # Whittier Fire
    TRUE ~ 0
  ))

# Trim it down for export
fire_trim <- dates_fire %>%
  select(site_code, Year, Month, watershed, fire)

# And export for MARSS script
saveRDS(fire_trim, "data_working/SBfire_edited_110321.rds")

# End of script.
