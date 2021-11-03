# SB Fire Data Assembly
# November 3, 2021
# Heili Lowman

# The following script will assemble the fire datasets available for the SBC LTER stream sites into a single data file for use in the MARSS analysis for the CRASS project.

# NOTE: This only contains fire info for 4 sites currently, that are being included in the initial MARSS analysis.

# Load packages
library(plyr) # needs to be loaded prior to dplyr
library(tidyverse) # contains dplyr
library(lubridate)
library(here)

# Load chemistry dataset from which dates will be pulled.
# Stream Chemistry
chem <- readRDS("data_working/SBchem_edited_102421.rds")

# Filter for 4 sites used in initial MARSS analysis
chem_filter <- chem %>%
  filter(site_code %in% c("AB00", "HO00", "MC06", "RG01"))

# And pull out only the site and date columns
dates_only <- chem_filter %>%
  select(site_code, Year, Month) %>%
  mutate(Day = 1) %>% # adding a day to create a full date
  mutate(Date = make_date(Year, Month, Day))

# Add watershed data to location dataframe using delineation at:
# https://databasin.org/maps/new/#datasets=6ad26ddb04ae4dbc9362303628270daf
dates_watersheds <- dates_only %>%
  mutate(watershed = factor(case_when(
    site_code == "ON02" ~ "Canada De Santa Anita",
    site_code == "GV01" ~ "Canada De La Gaviota",
    site_code == "HO00" ~ "Tajiguas Creek",
    site_code == "RG01" ~ "Tajiguas Creek",
    site_code == "TO02" ~ "Dos Pueblos Canyon",
    site_code == "BC02" ~ "Dos Pueblos Canyon",
    site_code == "DV01" ~ "Dos Pueblos Canyon",
    site_code == "SP02" ~ "San Pedro Creek",
    site_code == "AT07" ~ "Atascadero Creek",
    site_code == "AB00" ~ "Mission Creek",
    site_code == "MC06" ~ "Mission Creek",
    site_code == "RS02" ~ "Mission Creek")))

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
    # Arroyo Burro - Mission Creek Watershed
    site_code == "AB00" & Date >= "2008-11-01" & Date < "2009-05-01" ~ 1, # Tea Fire
    site_code == "AB00" & Date >= "2009-05-01" & Date < "2017-12-01" ~ 2, # Jesusita Fire
    site_code == "AB00" & Date >= "2017-12-01" ~ 3, # Thomas Fire
    # Arroyo Hondo - Tajiguas Creek Watershed
    site_code == "HO00" & Date >= "2016-06-01" & Date < "2017-07-01" ~ 1, # Sherpa Fire
    site_code == "HO00" & Date >= "2017-07-01" ~ 2, # Whittier Fire
    # Mission Creek - also Mission Creek Watershed
    site_code == "MC06" & Date >= "2008-11-01" & Date < "2009-05-01" ~ 1, # Tea Fire
    site_code == "MC06" & Date >= "2009-05-01" & Date < "2017-12-01" ~ 2, # Jesusita Fire
    site_code == "MC06" & Date >= "2017-12-01" ~ 3, # Thomas Fire
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
