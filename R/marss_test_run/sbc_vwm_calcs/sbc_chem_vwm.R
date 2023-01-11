# SB Stream Chemistry Calculations
# January 3, 2023
# Heili Lowman

# The following script will assemble the stream discharge datasets,
# as well as the volume weighted means by storm created by R. Aguilera,
# and calculate the remaining volume (discharge) weighted means to be
# used in the CRASS MARSS models for the Santa Barbara sites.

#### Setup ####

# Load packages
library(plyr) # needs to be loaded prior to dplyr
library(tidyverse) # contains dplyr
library(lubridate)
library(here)
library(naniar)
library(readxl)
library(dataRetrieval)

# Load discharge datasets.
# Note - not processing AT07 data since we do not have discharge data
# for the full record during which we have chemistry data.
# Also, John suggested removing this site from our analysis.
Qab <- read_csv("data_raw/Discharge/SBC_LTER/AB00_allyears_2018-09-30.csv")
Qgv <- read_csv("data_raw/Discharge/SBC_LTER/GV01_allyears_2017-09-30.csv")
Qho <- read_csv("data_raw/Discharge/SBC_LTER/HO00_allyears_2018-09-30.csv")
Qrs <- read_csv("data_raw/Discharge/SBC_LTER/RS02_allyears_2018-09-30.csv")

# The MC06 record was short on the LTER website so downloading directly
# from USGS using their dataRetrieval package.

siteNo <- "11119745" # https://waterdata.usgs.gov/ca/nwis/dv?referred_module=sw&site_no=11119745
pCode <- "00060" # discharge
start.date <- "2000-11-01"
end.date <- "2018-06-30"

rockynook <- readNWISdv(siteNumbers = siteNo,
                       parameterCd = pCode,
                       startDate = start.date,
                       endDate = end.date) # only daily data available

# exported units - cubic feet per second - need to be converted to lps
rockynook <- rockynook %>%
  mutate(discharge_lps = X_00060_00003*28.32)
# https://pubs.usgs.gov/wdr/wdr-id-03-3/ID03v3-InsideBackCover.pdf

Qmc <- rockynook %>%
  select(Date, discharge_lps)

# Load full chemistry dataset.
chem_reg <- read_csv("data_raw/sbclter_stream_chemistry_allyears_registered_stations_20190628.csv")

# Load Rosana's chemistry datasets.
nh4_storm <- read_xlsx("data_raw/SB_VWM_Storm_Conc/NH4stormfluxes.xlsx")
no3_storm <- read_xlsx("data_raw/SB_VWM_Storm_Conc/NO3stormfluxes.xlsx")
po4_storm <- read_xlsx("data_raw/SB_VWM_Storm_Conc/PO4stormfluxes.xlsx")

#### Filter Sites ####

all_chem_filtered <- chem_reg %>%
  # filter for 5 sites of interest
  filter(site_code %in% 
           c("AB00", "GV01", "HO00", "MC06", "RS02")) %>%
  # and replace all -999 values with NA
  replace_with_na_all(condition = ~.x == -999.0)

nh4_storm_filtered <- nh4_storm %>%
  filter(Site %in% 
           c("AB00", "GV01", "HO00", "MC06", "RS02"))

no3_storm_filtered <- no3_storm %>%
  filter(Site %in% 
           c("AB00", "GV01", "HO00", "MC06", "RS02"))

po4_storm_filtered <- po4_storm %>%
  filter(Site %in% 
           c("AB00", "GV01", "HO00", "MC06", "RS02"))

#### Delineate Storms ####

# I will be using similar aggregation methods between the wet/dry seasons,
# but I first need to delineate which measurements are duplicates.

# So, I need to create a list of dates that are represented in
# Rosana's datasets and remove them from the SBC LTER chem dataset.

# NH4 

# Take original dataset and trim times off dates.
nh4_storm_filtered <- nh4_storm_filtered %>%
  mutate(startpt = with_tz(startutc, "US/Pacific"),
         endpt = with_tz(endutc, "US/Pacific")) %>%
  mutate(start_date = date(startpt),
         end_date = date(endpt)) %>%
  mutate(interval = interval(start = start_date, 
                             end = end_date, 
                             tzone = tz("US/Pacific")))

# Then create smaller dataset with only sites and intervals.
nh4_site_storms <- nh4_storm_filtered %>%
  select(Site, interval)

# Rosana's data is slightly different for each analyte, so I also need
# to split up the chem data by analyte.

nh4_chem_filtered <- all_chem_filtered %>%
  select(site_code, timestamp_local, nh4_uM) %>%
  mutate(date = as_datetime(as.character(date(timestamp_local))))

# Iterate over sites and then over intervals

abtest1 <- nh4_chem_filtered %>% filter(site_code == "AB00")
abtest2 <- nh4_site_storms %>% filter(Site == "AB00")

## within ANY of the intervals of a list
dates <- abtest1$date
lst <- as.list(abtest2$interval)
testvector <- dates %within% lst
#ok so this works.

site_code <- c(rep("AB00",2767), rep("GV01",2011), 
               rep("HO00",1175), rep("MC06",1744), 
               rep("RS02",2130))

overlap_out <- c()

overlap <- data.frame() # New, blank df to receive inputs.

for(sitevar in c("AB00", "GV01", "HO00", "MC06", "RS02")) {
  
  # filter data by site desired
  df <- nh4_chem_filtered %>%
    filter(site_code == sitevar) # SB data
  
  rdf <- nh4_site_storms %>%
    filter(Site == sitevar) # Rosana's data
  
  # assign vectors of dates and storm intervals and compare
  dates <- df$date
  lst <- as.list(rdf$interval)
  overlap_out <- dates %within% lst
  
  # create df of output
  overlap_df <- data.frame(overlap_out)
  
  # assign site column to output
  overlap_df$Site <- c(rep(paste0(sitevar), length(overlap_df$overlap_out)))
  
  # bind with previous outputs
  overlap <- rbind(overlap, overlap_df)
  
} # OMG WUT

no3_chem_filtered <- all_chem_filtered %>%
  select(site_code, timestamp_local, no3_uM)

po4_chem_filtered <- all_chem_filtered %>%
  select(site_code, timestamp_local, po4_uM)

# Create dataset of all storms. Using NO3 dataset since it's the largest.
storms <- no3_storm_filtered %>%
  select(Site, Storm_name, startutc, endutc) %>%
  mutate(month = month(startutc))

# Appears to be storms October through June, so treating July - September
# as the non-rainy season.

# Add column to primary database to help with seasonal delineation.
all_chem_filtered$Season <- case_when(month(all_chem_filtered$timestamp_local) 
                                      %in% c("7", "8", "9") ~ "Dry",
                                      TRUE ~ "Rainy")

#### Non-rainy Monthly VWMs ####

# Filter out only dry season and add a new date column.
dry_chem <- all_chem_filtered %>%
  filter(Season == "Dry") %>%
  mutate(Date = date(timestamp_local))

# Create stand alone dataset of dates
dry_dates <- dry_chem %>%
  select(site_code, Date)

# Delineate groupings of discharge data
Q <- Qab %>%
  group_by(year, month) %>%
  mutate()


# calculate the difference from one day to the next
for(i in 2:nrow(d)){
  d$diff_time[i] = difftime(time1 = d$Date[i], time2 = d$Date[(i-1)],
                            units = "days")
}

# delineate sequenced time frames based on day to day differences
for(i in 2:nrow(d)){
  if(d$diff_time[i] < 14){
    d$seq[i] = d$seq[(i-1)]
  } else {
    d$seq[i] = d$seq[(i-1)] + 1
  }
}

# End of script.
