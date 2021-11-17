# SB Fire Data Assembly
# November 3, 2021
# Heili Lowman

# The following script will assemble the fire datasets available for the SBC LTER stream sites into a single data file for use in the MARSS analysis for the CRASS project.

# Load packages
library(plyr) # needs to be loaded prior to dplyr
library(tidyverse) # contains dplyr
library(lubridate)
library(here)

# Create base dataframe of dates and sites (from September 2002 to July 2016). This timeframe
# was chosen based on available precipitation data.
# Create sequence of dates
d <- seq(as.Date("2002/9/1"), by = "month", length.out = 166)

# Repeat 8 times for each site
d8 <- rep(d, times = 8)
d8df <- data.frame(d8)

# Create repeated sequence of sites
s <- rep(c("AB00", "AT07", "GV01", "HO00", "MC06", "RG01", "RS02", "SP02"), each=166)
sdf <- data.frame(s)

# Bind dates and sites together
dates <- cbind(d8df, sdf)
colnames(dates)<- c("date","site")

# Add watershed data to location dataframe using delineation at:
# https://databasin.org/maps/new/#datasets=6ad26ddb04ae4dbc9362303628270daf
dates_watersheds <- dates %>%
  mutate(watershed = factor(case_when(
    site == "GV01" ~ "Canada De La Gaviota",
    site == "HO00" ~ "Tajiguas Creek",
    site == "RG01" ~ "Tajiguas Creek",
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

# Now, add in columns for the fires
# 0 - denotes pre-fire months
# 1 - denotes post-fire months
# Each fire will get its own column, so that effects can potentially be additive
# Or so we can combine decay functions later

dates_fire <- dates_watersheds %>%
  mutate(AB00_Tea = ifelse(site == "AB00" & date >= "2008-11-01", 1, 0),
    AB00_Jesusita = ifelse(site == "AB00" & date >= "2009-05-01", 1, 0),
    AT07_Jesusita = ifelse(site == "AT07" & date >= "2009-05-01", 1, 0),
    GV01_Gaviota = ifelse(site == "GV01" & date >= "2004-06-01", 1, 0),
    HO00_Gaviota = ifelse(site == "HO00" & date >= "2004-06-01", 1, 0),
    HO00_Sherpa = ifelse(site == "HO00" & date >= "2016-06-01", 1, 0),
    HO00_Whittier = ifelse(site == "HO00" & date >= "2017-07-01", 1, 0),
    MC06_Tea = ifelse(site == "MC06" & date >= "2008-11-01", 1, 0),
    MC06_Jesusita = ifelse(site == "MC06" & date >= "2009-05-01", 1, 0),
    RG01_Gaviota = ifelse(site == "RG01" & date >= "2004-06-01", 1, 0),
    RG01_Sherpa = ifelse(site == "RG01" & date >= "2016-06-01", 1, 0),
    RG01_Whittier = ifelse(site == "RG01" & date >= "2017-07-01", 1, 0),
    RS02_Tea = ifelse(site == "RS02" & date >= "2008-11-01", 1, 0),
    RS02_Jesusita = ifelse(site == "RS02" & date >= "2009-05-01", 1, 0),
    SP02_Gap = ifelse(site == "SP02" & date >= "2008-07-01", 1, 0))

# And export for MARSS script
saveRDS(dates_fire, "data_working/SBfire_edited_111721.rds")

# End of script.
