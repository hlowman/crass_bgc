# SB Stream Chemistry Assembly
# October 24, 2021
# Heili Lowman

# The following script will assemble the stream chemistry datasets available for the SBC LTER stream sites into a single data file for use in the MARSS analysis for the CRASS project.

# Load packages
library(plyr) # needs to be loaded prior to dplyr
library(tidyverse) # contains dplyr
library(lubridate)
library(here)
library(naniar)

# Load dataset
chem_reg <- read_csv("data_raw/sbclter_stream_chemistry_allyears_registered_stations_20190628.csv")
# After reviewing the metadata on the SBC LTER website, chose not to include non-registered sites,
# because there is no geolocation data available/published for them, like there is for the registered
# sites.
# Doing quick check of data availability though.
# chem_nreg <- read_csv("data_raw/sbclter_stream_chemistry_allyears_non_registered_stations_20190628.csv")
 
# chem_nreg <- chem_nreg %>%
#   mutate(Date = ymd_hms(timestamp_local)) %>%
#   mutate(Year = year(Date))

# timespans <- chem_nreg %>%
#   group_by(site_code) %>%
#   summarize(minYear = min(Year), maxYear = max(Year)) %>%
#   ungroup() %>%
#   mutate(Span = maxYear - minYear)

# multiYear_sites <- timespans %>%
#   filter(Span > 1)

# Should loop back with John to see if this Ventura River (VR)
# data would be suitable for the large river analysis?

# Note: -999 is the "NA" record used by the LTER
# Make edits for data assembly purposes
chem_full_ed <- chem_reg %>% # can't remove records because in wide format
  replace_with_na_all(condition = ~.x == -999.0) %>% # so instead replace with NA
  mutate(DateTime = ymd_hms(timestamp_local)) %>% # format dates
  mutate(Year = year(DateTime), Month = month(DateTime))

# Calculate mean monthly concentrations at all sites.
# Note: For this first iteration, not weighted by discharge.
chem_full_monthly <- chem_full_ed %>%
  group_by(site_code, Year, Month) %>%
  summarize(mean_nh4_uM = mean(nh4_uM, na.rm = TRUE),
            mean_no3_uM = mean(no3_uM, na.rm = TRUE),
            mean_po4_uM = mean(po4_uM, na.rm = TRUE),
            mean_tdn_uM = mean(tdn_uM, na.rm = TRUE),
            mean_tdp_uM = mean(tdp_uM, na.rm = TRUE),
            mean_tpc_uM = mean(tpc_uM, na.rm = TRUE),
            mean_tpn_uM = mean(tpn_uM, na.rm = TRUE),
            mean_tpp_uM = mean(tpp_uM, na.rm = TRUE),
            mean_tss_mgL = mean(tss_mgperLiter, na.rm = TRUE),
            mean_cond_uScm = mean(spec_cond_uSpercm, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(site = factor(site_code))

# And export for MARSS script
saveRDS(chem_full_monthly, "data_working/SBchem_edited_110721.rds")

# End of script.
