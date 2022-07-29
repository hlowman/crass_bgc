# SB Fire Data Assembly
# November 3, 2021
# Heili Lowman

# The following script will assemble the fire datasets available for the SBC LTER stream sites into a single data file for use in the MARSS analysis for the CRASS project.

# And edited on 7/29/2022 by Alex to add fire areas and burn percent of watersheds, set fire effect to one column with an effect lasting 2 months, and remove decay effects. 

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

# load fire area data
fire_area = read.csv("data_raw/sbc_vc_ws_fire_areas.csv")
unique(fire_area$fire_name) #use these names below

# And here is a list of the wildfire events in the SBC:
# Fire Name (Ignition Date) - affected watersheds
# Jesusita (5/5/2009) - Rattlesnake, Mission, Arroyo Burro, Atascadero
# Gap (7/1/2008) - San Pedro
# Gaviota (6/5/2004) - Arroyo Hondo, Gaviota
# Tea (11/13/2008) - Rattlesnake, Mission
# Sherpa (6/15/2016) - Refugio
# Cave, Whittier, and Thomas Fire are too late to be included.

# Now, add in columns for the fires
# 1 - denotes month of ignition + one additional month

sites_JESUSITA = c("RS02", "MC06", "AB00", "AT07")
dates_JESUSITA = as.Date(c("2009-05-01","2009-06-01"))
sites_GAP = c("SP02")
dates_GAP = as.Date(c("2008-07-01","2008-08-01"))
sites_GAVIOTA = c("HO00","GV01")
dates_GAVIOTA = as.Date(c("2004-06-01","2004-07-01"))
sites_TEA = c("RS02","MC06")
dates_TEA = as.Date(c("2008-11-01","2008-12-01"))
sites_SHERPA = c("RG01")
dates_SHERPA = as.Date(c("2016-06-01","2016-07-01"))

dates_fire <- dates %>%
  mutate(fire_pa = ifelse(site %in% sites_JESUSITA & date %in% dates_JESUSITA, 1,
                          ifelse(site %in% sites_GAP & date %in% sites_GAP,1,
                                 ifelse(site %in% sites_GAVIOTA & date %in% dates_GAVIOTA,1,
                                        ifelse(site %in% sites_TEA & date %in% dates_TEA,1,
                                               ifelse(site %in% sites_SHERPA & date %in% dates_SHERPA,1,0))))))%>%
  mutate(fire_ID = ifelse(site %in% sites_JESUSITA & date %in% dates_JESUSITA, "JESUSITA",
                          ifelse(site %in% sites_GAP & date %in% sites_GAP,"GAP",
                                 ifelse(site %in% sites_GAVIOTA & date %in% dates_GAVIOTA,"GAVIOTA",
                                        ifelse(site %in% sites_TEA & date %in% dates_TEA,"TEA",
                                               ifelse(site %in% sites_SHERPA & date %in% dates_SHERPA,"SHERPA",NA))))))

# match names to dates_fire_vcnp
names(fire_area) = c("site","ws_area_m2","fire_ID","ig_date","ws_fire_area_m2", "fire_perc_ws")
# join watershed area to dates_fire_vcnp
fire_1 = left_join(dates_fire, unique(fire_area[,1:2]), by=c("site"))
# join fire area info to dates_fire_vcnp
fire_2 = left_join(fire_1, fire_area[,c(1,3:6)], by=c("site","fire_ID"))
# replace NAs in fire info cols with zero
fire_2[,7:8][is.na(fire_2[,7:8])] = 0

# And export for MARSS script
#saveRDS(dates_fire, "data_working/SBfire_edited_072822.rds")
saveRDS(dates_fire, "data_working/SBfire_edited_072922.rds")

# End of script.
