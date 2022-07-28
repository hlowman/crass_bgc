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

# Re-did watershed assignments on 7/28/22 based on
# watershed delineated back from sampling location.
dates_watersheds <- dates %>%
  mutate(watershed = factor(case_when(
    site == "GV01" ~ "Gaviota",
    site == "HO00" ~ "Arroyo Hondo",
    site == "RG01" ~ "Refugio",
    site == "SP02" ~ "San Pedro",
    site == "AT07" ~ "Atascadero",
    site == "AB00" ~ "Arroyo Burro",
    site == "MC06" ~ "Mission",
    site == "RS02" ~ "Rattlesnake")))

# And here is a list of the wildfire events in the SBC:
# Fire Name (Ignition Date) - affected watersheds
# Gaviota (6/5/2004) - Arroyo Hondo, Gaviota
# Sherpa (6/15/2016) - Refugio
# Gap (7/1/2008) - San Pedro
# Jesusita (5/5/2009) - Rattlesnake, Mission, Arroyo Burro, Atascadero
# Tea (11/13/2008) - Rattlesnake, Mission

# Cave, Whittier, and Thomas Fire are too late to be included.

# Now, add in columns for the fires
# 0 - denotes no-fire months
# 1 - denotes fire ignition months
# Each fire will get its own column, so that effects can potentially be additive
# Or so we can combine decay functions later

dates_fire <- dates_watersheds %>%
  mutate(AB00_Jesusita = ifelse(site == "AB00" & date == "2009-05-01", 1, 0),
    AT07_Jesusita = ifelse(site == "AT07" & date == "2009-05-01", 1, 0),
    GV01_Gaviota = ifelse(site == "GV01" & date == "2004-06-01", 1, 0),
    HO00_Gaviota = ifelse(site == "HO00" & date == "2004-06-01", 1, 0),
    MC06_Tea = ifelse(site == "MC06" & date == "2008-11-01", 1, 0),
    MC06_Jesusita = ifelse(site == "MC06" & date == "2009-05-01", 1, 0),
    RG01_Sherpa = ifelse(site == "RG01" & date == "2016-06-01", 1, 0),
    RS02_Tea = ifelse(site == "RS02" & date == "2008-11-01", 1, 0),
    RS02_Jesusita = ifelse(site == "RS02" & date == "2009-05-01", 1, 0),
    SP02_Gap = ifelse(site == "SP02" & date == "2008-07-01" , 1, 0)) # %>%
  # also adding in 4 year decay term for each of the fires
  #mutate(AB00_Tea_d = NA,
         # AB00_Jesusita_d = NA,
         # AT07_Jesusita_d = NA,
         # GV01_Gaviota_d = NA,
         # HO00_Gaviota_d = NA,
         # HO00_Sherpa_d = NA,
         # MC06_Tea_d = NA,
         # MC06_Jesusita_d = NA,
         # RG01_Gaviota_d = NA,
         # RG01_Sherpa_d = NA,
         # RS02_Tea_d = NA,
         # RS02_Jesusita_d = NA,
         # SP02_Gap_d = NA)

# values <- rev(seq(1, 48, by = 1))
# decay <- exp(values)
# 
# dates_fire$AB00_Tea_d[75:122] <- decay
# dates_fire$AB00_Jesusita_d[81:128] <- decay
# dates_fire$AT07_Jesusita_d[247:294] <- decay
# dates_fire$GV01_Gaviota_d[354:401] <- decay
# dates_fire$HO00_Gaviota_d[520:567] <- decay
# dates_fire$HO00_Sherpa_d[664] <- decay[1]
# dates_fire$MC06_Tea_d[739:786] <- decay
# dates_fire$MC06_Jesusita_d[745:792] <- decay
# dates_fire$RG01_Gaviota_d[852:899] <- decay
# dates_fire$RG01_Sherpa_d[996] <- decay[1]
# dates_fire$RS02_Tea_d[1071:1118] <- decay
# dates_fire$RS02_Jesusita_d[1077:1124] <- decay
# dates_fire$SP02_Gap_d[1233:1280] <- decay

dates_fire[is.na(dates_fire)] = 0

# And export for MARSS script
saveRDS(dates_fire, "data_working/SBfire_edited_072822.rds")

# End of script.
