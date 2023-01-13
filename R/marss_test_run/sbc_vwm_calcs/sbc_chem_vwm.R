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

#### Assemble Discharge Data ####

# Trim SB discharge datasets to match USGS format.
# a.k.a. take daily averages.

Qab_daily <- Qab %>%
  mutate(timestamp_date = date(timestamp_local)) %>%
  filter(discharge_lps != -999.0) %>%
  group_by(timestamp_date) %>%
  summarize(discharge_lps = round(mean(discharge_lps, na.rm = TRUE),
                                  digits = 4)) %>%
  ungroup() %>%
  mutate(site_code = "AB00")

Qgv_daily <- Qgv %>%
  mutate(timestamp_date = date(timestamp_local)) %>%
  filter(discharge_lps != -999.0) %>%
  group_by(timestamp_date) %>%
  summarize(discharge_lps = round(mean(discharge_lps, na.rm = TRUE),
                                  digits = 4)) %>%
  ungroup() %>%
  mutate(site_code = "GV01")

Qho_daily <- Qho %>%
  mutate(timestamp_date = date(timestamp_local)) %>%
  filter(discharge_lps != -999.0) %>%
  group_by(timestamp_date) %>%
  summarize(discharge_lps = round(mean(discharge_lps, na.rm = TRUE),
                                  digits = 4)) %>%
  ungroup() %>%
  mutate(site_code = "HO00")

Qmc_daily <- Qmc %>%
  mutate(site_code = "MC06") %>%
  rename("timestamp_date" = "Date")

Qrs_daily <- Qrs %>%
  mutate(timestamp_date = date(timestamp_local)) %>%
  filter(discharge_lps != -999.0) %>%
  group_by(timestamp_date) %>%
  summarize(discharge_lps = round(mean(discharge_lps, na.rm = TRUE),
                                  digits = 4)) %>%
  ungroup() %>%
  mutate(site_code = "RS02")

Q_all_daily <- rbind(Qab_daily, Qgv_daily, Qho_daily, Qmc_daily, Qrs_daily)

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
# Rosana's datasets and remove them from the SBC LTER chem and discharge data.

#### NH4 ####

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
#ok so this works as a standalone test. now to iterate.

overlap_nh4 <- data.frame() # New, blank df to receive inputs.

for(sitevar in c("AB00", "GV01", "HO00", "MC06", "RS02")) {
  
  # filter data by site desired
  df <- nh4_chem_filtered %>%
    filter(site_code == sitevar) # SB data
  
  rdf <- nh4_site_storms %>%
    filter(Site == sitevar) # Rosana's data
  
  # assign vectors of dates and storm intervals and compare
  dates <- df$date
  lst <- as.list(rdf$interval)
  overlap_storm <- dates %within% lst
  
  # create df of output
  overlap_df <- data.frame(overlap_storm)
  
  # carry over date information for easier binding later
  overlap_df$date <- dates
  
  # assign site column to output
  overlap_df$Site <- c(rep(paste0(sitevar), length(overlap_df$overlap_storm)))
  
  # bind with previous outputs
  overlap_nh4 <- rbind(overlap_nh4, overlap_df)
  
} # OMG WUT

# Arrange nutrient dataset by site.
nh4_chem_filtered <- arrange(nh4_chem_filtered, site_code)

# Add results of for loop.
nh4_chem_combined <- cbind(nh4_chem_filtered, overlap_nh4)

# Now to detect storm overlap with available discharge data.
# Note, I'll be doing this separately for each analyte since the data/storms
# in Rosana's dataset differ by analyte.

overlap_nh4_Q <- data.frame() # New, blank df to receive inputs.

for(sitevar in c("AB00", "GV01", "HO00", "MC06", "RS02")) {
  
  # filter data by site desired
  df <- Q_all_daily %>%
    filter(site_code == sitevar) # SB data
  
  rdf <- nh4_site_storms %>%
    filter(Site == sitevar) # Rosana's data
  
  # assign vectors of dates and storm intervals and compare
  dates <- df$timestamp_date
  lst <- as.list(rdf$interval)
  overlap_storm <- dates %within% lst
  
  # create df of output
  overlap_df <- data.frame(overlap_storm)
  
  # carry over date information for easier binding later
  overlap_df$date <- dates
  
  # assign site column to output
  overlap_df$Site <- c(rep(paste0(sitevar), length(overlap_df$overlap_storm)))
  
  # bind with previous outputs
  overlap_nh4_Q <- rbind(overlap_nh4_Q, overlap_df)
  
} # Yay!

# Add results of for loop.
Q_all_combined <- cbind(Q_all_daily, overlap_nh4_Q)

# Next, filter out overlapping days from nh4 and Q datasets.

nh4_not_storm <- nh4_chem_combined[,c(1,2,4,3,5)] %>%
  filter(overlap_storm == FALSE)

Q_not_storm <- Q_all_combined[,c(3,1,2,4)] %>%
  filter(overlap_storm == FALSE) %>%
  mutate(discharge_L = discharge_lps*86400) # number of seconds/day

# Calculate mean monthly averages of NH4 for non-storm SB days.
nh4_not_storm_monthly <- nh4_not_storm %>%
  mutate(Year = year(date),
         Month = month(date)) %>%
  group_by(site_code, Year, Month) %>%
  summarize(mean_nh4_uM = mean(nh4_uM, na.rm = TRUE)) %>%
  ungroup()

# Calculate cumulative monthly discharges for non-storm SB days.
Q_not_storm_monthly <- Q_not_storm %>%
  mutate(Year = year(timestamp_date),
         Month = month(timestamp_date)) %>%
  group_by(site_code, Year, Month) %>%
  summarize(cum_Q_L = sum(discharge_L)) %>%
  ungroup()

# Combine nh4 and Q datasets.
SB_nh4_Q <- left_join(nh4_not_storm_monthly, Q_not_storm_monthly) %>%
  mutate(c_x_Q = mean_nh4_uM*cum_Q_L)

# And now to summarize Rosana's data.
# First, need to pull out columns of interest, and add the vwm*Q calculation.
nh4_storm_trim <- nh4_storm_filtered %>%
  select(Site, start_date, Storm_name, vwm_micromol, totalQ_L) %>%
  mutate(vwm_x_Q = vwm_micromol*totalQ_L)
  
# Next, need to calculate monthly totals for both numerator (c*V) and 
# denominator (V) of Williams equation.
nh4_storm_summ <- nh4_storm_trim %>%
  mutate(Year = year(start_date),
         Month = month(start_date)) %>%
  group_by(Site, Year, Month) %>%
  summarize(cum_nh4_storms = sum(vwm_x_Q),
            cum_Q_storms = sum(totalQ_L)) %>%
  ungroup() %>%
  rename("site_code" = "Site")

# Combine LTER non-storm data with Rosana's storm data.
nh4_all <- full_join(SB_nh4_Q, nh4_storm_summ,
                     by = c("site_code", "Year", "Month"))

# But need to replace all the storm NA values with 0 because it's throwing
# off the vwm calculation below.
nh4_all$cum_nh4_storms[is.na(nh4_all$cum_nh4_storms)] <- 0
nh4_all$cum_Q_storms[is.na(nh4_all$cum_Q_storms)] <- 0 

# Omg, and is this the last calculation????
# Calculate VOLUME WEIGHTED MEANS.
nh4_all <- nh4_all %>%
  mutate(vwm_monthly = (cum_nh4_storms + c_x_Q)/(cum_Q_storms + cum_Q_L))

# Quick plot to see how this looks:
(plot1 <- ggplot(nh4_all %>%
                  mutate(Day = 1) %>%
                  mutate(Date = make_date(Year, Month, Day)), 
                aes(x = Date, y = vwm_monthly)) +
  geom_point() +
  geom_line() +
  labs(y = "V.W.M. NH4 uM") +
  theme_bw() +
  facet_wrap(.~site_code, nrow = 2, scales = "free"))

# ggsave(plot1,
#        filename = "figures/SB_NH4_VWM_011223.jpg",
#        width = 30,
#        height = 15,
#        units = "cm")

#### NO3 ####

# Take original dataset from Rosana and trim times off dates.
no3_storm_filtered <- no3_storm_filtered %>%
  mutate(startpt = with_tz(startutc, "US/Pacific"),
         endpt = with_tz(endutc, "US/Pacific")) %>%
  mutate(start_date = date(startpt),
         end_date = date(endpt)) %>%
  mutate(interval = interval(start = start_date, 
                             end = end_date, 
                             tzone = tz("US/Pacific")))

# Then create smaller dataset with only sites and intervals.
no3_site_storms <- no3_storm_filtered %>%
  select(Site, interval)

# Rosana's data is slightly different for each analyte, so I also need
# to split up the chem data by analyte.
no3_chem_filtered <- all_chem_filtered %>%
  select(site_code, timestamp_local, no3_uM) %>%  
  mutate(date = as_datetime(as.character(date(timestamp_local))))

# Iterate over sites and then over intervals
overlap_no3 <- data.frame() # New, blank df to receive inputs.

for(sitevar in c("AB00", "GV01", "HO00", "MC06", "RS02")) {
  
  # filter data by site desired
  df <- no3_chem_filtered %>%
    filter(site_code == sitevar) # SB data
  
  rdf <- no3_site_storms %>%
    filter(Site == sitevar) # Rosana's data
  
  # assign vectors of dates and storm intervals and compare
  dates <- df$date
  lst <- as.list(rdf$interval)
  overlap_storm <- dates %within% lst
  
  # create df of output
  overlap_df <- data.frame(overlap_storm)
  
  # carry over date information for easier binding later
  overlap_df$date <- dates
  
  # assign site column to output
  overlap_df$Site <- c(rep(paste0(sitevar), 
                           length(overlap_df$overlap_storm)))
  
  # bind with previous outputs
  overlap_no3 <- rbind(overlap_no3, overlap_df)
  
} # Yes!!

# Arrange nutrient dataset by site.
no3_chem_filtered <- arrange(no3_chem_filtered, site_code)

# Add results of for loop.
no3_chem_combined <- cbind(no3_chem_filtered, overlap_no3)

# Now to detect storm overlap with available discharge data.
# Note, I'll be doing this separately for each analyte since the data/storms
# in Rosana's dataset differ by analyte.

overlap_no3_Q <- data.frame() # New, blank df to receive inputs.

for(sitevar in c("AB00", "GV01", "HO00", "MC06", "RS02")) {
  
  # filter data by site desired
  df <- Q_all_daily %>%
    filter(site_code == sitevar) # SB data
  
  rdf <- no3_site_storms %>%
    filter(Site == sitevar) # Rosana's data
  
  # assign vectors of dates and storm intervals and compare
  dates <- df$timestamp_date
  lst <- as.list(rdf$interval)
  overlap_storm <- dates %within% lst
  
  # create df of output
  overlap_df <- data.frame(overlap_storm)
  
  # carry over date information for easier binding later
  overlap_df$date <- dates
  
  # assign site column to output
  overlap_df$Site <- c(rep(paste0(sitevar), 
                           length(overlap_df$overlap_storm)))
  
  # bind with previous outputs
  overlap_no3_Q <- rbind(overlap_no3_Q, overlap_df)
  
} # Yay!

# Add results of for loop.
Q_all_combined_no3 <- cbind(Q_all_daily, overlap_no3_Q)

# Next, filter out overlapping days from nh4 and Q datasets.

no3_not_storm <- no3_chem_combined[,c(1,2,4,3,5)] %>%
  filter(overlap_storm == FALSE)

Q_not_storm_no3 <- Q_all_combined_no3[,c(3,1,2,4)] %>%
  filter(overlap_storm == FALSE) %>%
  mutate(discharge_L = discharge_lps*86400) # number of seconds/day

# Calculate mean monthly averages of NO3 for non-storm SB days.
no3_not_storm_monthly <- no3_not_storm %>%
  mutate(Year = year(date),
         Month = month(date)) %>%
  group_by(site_code, Year, Month) %>%
  summarize(mean_no3_uM = mean(no3_uM, na.rm = TRUE)) %>%
  ungroup()

# Calculate cumulative monthly discharges for non-storm SB days.
Q_not_storm_monthly_no3 <- Q_not_storm_no3 %>%
  mutate(Year = year(timestamp_date),
         Month = month(timestamp_date)) %>%
  group_by(site_code, Year, Month) %>%
  summarize(cum_Q_L = sum(discharge_L)) %>%
  ungroup()

# Combine nh4 and Q datasets.
SB_no3_Q <- left_join(no3_not_storm_monthly, Q_not_storm_monthly_no3) %>%
  mutate(c_x_Q = mean_no3_uM*cum_Q_L)

# And now to summarize Rosana's data.
# First, need to pull out columns of interest, and add the vwm*Q calculation.
no3_storm_trim <- no3_storm_filtered %>%
  select(Site, start_date, Storm_name, vwm_micromol, totalQ_L) %>%
  mutate(vwm_x_Q = vwm_micromol*totalQ_L)

# Next, need to calculate monthly totals for both numerator (c*V) and 
# denominator (V) of Williams equation.
no3_storm_summ <- no3_storm_trim %>%
  mutate(Year = year(start_date),
         Month = month(start_date)) %>%
  group_by(Site, Year, Month) %>%
  summarize(cum_no3_storms = sum(vwm_x_Q),
            cum_Q_storms = sum(totalQ_L)) %>%
  ungroup() %>%
  rename("site_code" = "Site")

# Combine LTER non-storm data with Rosana's storm data.
no3_all <- full_join(SB_no3_Q, no3_storm_summ,
                     by = c("site_code", "Year", "Month"))

# But need to replace all the storm NA values with 0 because it's throwing
# off the vwm calculation below.
no3_all$cum_no3_storms[is.na(no3_all$cum_no3_storms)] <- 0
no3_all$cum_Q_storms[is.na(no3_all$cum_Q_storms)] <- 0 

# Calculate VOLUME WEIGHTED MEANS.
no3_all <- no3_all %>%
  mutate(vwm_monthly = (cum_no3_storms + c_x_Q)/(cum_Q_storms + cum_Q_L))

# Quick plot to see how this looks:
(plot2 <- ggplot(no3_all %>%
                   mutate(Day = 1) %>%
                   mutate(Date = make_date(Year, Month, Day)), 
                 aes(x = Date, y = vwm_monthly)) +
    geom_point() +
    geom_line() +
    labs(y = "V.W.M. NO3 uM") +
    theme_bw() +
    facet_wrap(.~site_code, nrow = 2, scales = "free"))

# ggsave(plot2,
#        filename = "figures/SB_NO3_VWM_011323.jpg",
#        width = 30,
#        height = 15,
#        units = "cm")

#### PO4 ####

po4_chem_filtered <- all_chem_filtered %>%
  select(site_code, timestamp_local, po4_uM)

# End of script.
