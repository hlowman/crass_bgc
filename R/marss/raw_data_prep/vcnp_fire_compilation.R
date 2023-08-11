# Valles Caldera National Preserve (vcnp) Fire Data Assembly
# November 3, 2021
# Betsy Summers following Heili Lowman's workflow

# The following script will assemble the fire datasets available for the VCNP 
# stream sites into a single data file for use in the MARSS analysis for 
# the CRASS project.

# the main fires of interest that overlap with the grab water quality data set
# are (date = month of ignition):
# 1) Las Conchas fire (Start date: 26 June, 2011)
# 2) Thompson Ridge fire (start date: 31 May, 2013)

## Load packages
library(plyr) # needs to be loaded prior to dplyr
library(tidyverse) # contains dplyr
library(lubridate)
library(here)

# Create base dataframe of dates and sites (from April 2005 to October 2020). This timeframe
# was chosen based on available precipitation and chem data.
# Create sequence of dates
d <- seq(as.Date("2005/4/1"), by = "month", length.out = 186)

# Repeat 7 times for each site
d7 <- rep(d, times = 7)
d7df <- data.frame(d7)

# Create repeated sequence of sites
s <- rep(c("EFJ", "RSAW", "RSA", "IND", "IND_BB", "RED", "SULF"), each = 186)
sdf <- data.frame(s)

# Bind dates and sites together
dates <- cbind(d7df, sdf)
colnames(dates) <- c("date","site")

# load fire area data
fire_area = read.csv("data_raw/sbc_vc_ws_fire_areas.csv")

## Fire names and start dates of burn - sites affected
# Las Conchas (2011-06-26) - EFJ, RSAW, RSA, IND, IND_BB
# Thompson Ridge (2013-05-31) - RED, EFJ, RSAW, SULF
# Prescribed burn (2016-05-11) - EFJ [deemed too small to include]

## Add in columns for the fires based on start dates:
# 1 - denotes month of ignition + one additional month
sites_THOMPSON = c("RED", "EFJ", "RSAW", "SULF")
dates_THOMPSON = as.Date(c("2013-05-01","2013-06-01"))
sites_CONCHAS = c("EFJ", "RSAW", "RSA", "IND", "IND_BB")
dates_CONCHAS = as.Date(c("2011-06-01","2011-07-01"))

dates_fire_vcnp <- dates %>%
  mutate(fire_pa = ifelse(site %in% sites_THOMPSON & date %in% dates_THOMPSON, 1,
                          ifelse(site %in% sites_CONCHAS & date %in% dates_CONCHAS,1,0)))%>%
  mutate(fire_ID = ifelse(site %in% sites_THOMPSON & date %in% dates_THOMPSON,"THOMPSON RIDGE",
                          ifelse(site %in% sites_CONCHAS & date %in% dates_CONCHAS,"LAS CONCHAS",NA)))

# match names to dates_fire_vcnp
names(fire_area) = c("site","ws_area_m2","fire_ID","ig_date","ws_fire_area_m2", "fire_perc_ws")
# join watershed area to dates_fire_vcnp
fire_vcnp_1 = left_join(dates_fire_vcnp, unique(fire_area[,1:2]), by=c("site"))
# join fire area info to dates_fire_vcnp
fire_vcnp_2 = left_join(fire_vcnp_1, fire_area[,c(1,3:6)], by=c("site","fire_ID"))
# replace NAs in fire info cols with zero
fire_vcnp_2[,7:8][is.na(fire_vcnp_2[,7:8])] = 0

# And export for MARSS script
saveRDS(fire_vcnp_2, "data_working/VCNPfire_edited_081123.rds")

# End of script.
