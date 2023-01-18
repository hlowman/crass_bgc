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
# Using limits of detection / 2 based on the methods detailed here:
# https://sbclter.msi.ucsb.edu/external/Land/Protocols/Stream_Chemistry/Melack_Schimel_20090529_SBCLTER_Laboratory_Analyses.pdf

# Make edits for data assembly purposes
chem_full_ed <- chem_reg %>% # can't remove records because in wide format
  replace_with_na_all(condition = ~.x == -999.0) %>% # so instead replace with NA
  mutate(nh4_uM_lod = ifelse(nh4_uM < 0.25, 0.25, nh4_uM), # NH4 LOD = 0.5uM
         no3_uM_lod = ifelse(no3_uM < 0.25, 0.25, no3_uM), # NO3 LOD = 0.5uM
         po4_uM_lod = ifelse(po4_uM < 0.15, 0.15, po4_uM), # PO4 LOD = 0.3uM
         tdn_uM_lod = ifelse(tdn_uM < 0.5, 0.5, tdn_uM), # TDN LOD = 1uM
         tdp_uM_lod = ifelse(tdp_uM < 0.5, 0.5, tdp_uM), # TDP LOD = 1uM
         # tpc_uM_lod = ifelse(tpc_uM <= 0, 0.25, tpc_uM), # TPC LOD = 2ug - need to convert later
         # tpn_uM_lod = ifelse(tpn_uM <= 0, 0.25, tpn_uM), # TPN LOD = 2ug - need to convert later
         tpp_uM_lod = ifelse(tpp_uM < 0.5, 0.5, tpp_uM) # TPP LOD = 1uM
         ) %>%
  mutate(DateTime = ymd_hms(timestamp_local)) %>% # format dates
  mutate(Year = year(DateTime), Month = month(DateTime)) %>%
  mutate(site = factor(site_code))

# Calculate mean monthly concentrations at all sites.
# Note: For this first iteration, not weighted by discharge.
chem_full_monthly <- chem_full_ed %>%
  group_by(site, Year, Month) %>%
  summarize(mean_nh4_uM = mean(nh4_uM_lod, na.rm = TRUE),
            mean_no3_uM = mean(no3_uM_lod, na.rm = TRUE),
            mean_po4_uM = mean(po4_uM_lod, na.rm = TRUE),
            mean_tdn_uM = mean(tdn_uM_lod, na.rm = TRUE),
            mean_tdp_uM = mean(tdp_uM_lod, na.rm = TRUE),
            mean_tpc_uM = mean(tpc_uM, na.rm = TRUE),
            mean_tpn_uM = mean(tpn_uM, na.rm = TRUE),
            mean_tpp_uM = mean(tpp_uM_lod, na.rm = TRUE),
            mean_tss_mgL = mean(tss_mgperLiter, na.rm = TRUE),
            mean_cond_uScm = mean(spec_cond_uSpercm, na.rm = TRUE)) %>%
  ungroup()

# And export for MARSS script
write_csv(chem_full_monthly, "data_working/SBchem_edited_120321.csv")
saveRDS(chem_full_monthly, "data_working/SBchem_edited_120321.rds")

# Additional plots made July 28,2022

chem_oursites_cond <- chem_full_ed %>%
  filter(site_code %in% 
           c("RS02", "MC06", "AB00", "AT07", "SP02", "RG01", "HO00", "GV01")) %>%
  select(site_code, DateTime, spec_cond_uSpercm)

(RS02_plot <- ggplot(chem_oursites_cond %>% filter(site_code == "RS02")) +
  geom_point(aes(x = DateTime, y = spec_cond_uSpercm)) +
  geom_vline(xintercept = as.POSIXct(as.Date("2008-11-13")), color = "red") +
  geom_vline(xintercept = as.POSIXct(as.Date("2009-05-05")), color = "red") +
    labs(x = "Date", y = "Conductivity (uS/cm)", title = "RS02") +
    theme_bw())

(MC06_plot <- ggplot(chem_oursites_cond %>% filter(site_code == "MC06")) +
    geom_point(aes(x = DateTime, y = spec_cond_uSpercm)) +
    geom_vline(xintercept = as.POSIXct(as.Date("2008-11-13")), color = "red") +
    geom_vline(xintercept = as.POSIXct(as.Date("2009-05-05")), color = "red") +
    labs(x = "Date", y = "Conductivity (uS/cm)", title = "MC06") +
    theme_bw())

(AB00_plot <- ggplot(chem_oursites_cond %>% filter(site_code == "AB00")) +
    geom_point(aes(x = DateTime, y = spec_cond_uSpercm)) +
    geom_vline(xintercept = as.POSIXct(as.Date("2009-05-05")), color = "red") +
    labs(x = "Date", y = "Conductivity (uS/cm)", title = "AB00") +
    theme_bw())

(AT07_plot <- ggplot(chem_oursites_cond %>% filter(site_code == "AT07")) +
    geom_point(aes(x = DateTime, y = spec_cond_uSpercm)) +
    geom_vline(xintercept = as.POSIXct(as.Date("2009-05-05")), color = "red") +
    labs(x = "Date", y = "Conductivity (uS/cm)", title = "AT07") +
    theme_bw())

(SP02_plot <- ggplot(chem_oursites_cond %>% filter(site_code == "SP02")) +
    geom_point(aes(x = DateTime, y = spec_cond_uSpercm)) +
    geom_vline(xintercept = as.POSIXct(as.Date("2008-07-01")), color = "red") +
    labs(x = "Date", y = "Conductivity (uS/cm)", title = "SP02") +
    theme_bw())

(RG01_plot <- ggplot(chem_oursites_cond %>% filter(site_code == "RG01")) +
    geom_point(aes(x = DateTime, y = spec_cond_uSpercm)) +
    geom_vline(xintercept = as.POSIXct(as.Date("2016-06-15")), color = "red") +
    labs(x = "Date", y = "Conductivity (uS/cm)", title = "RG01") +
    theme_bw())

(HO00_plot <- ggplot(chem_oursites_cond %>% filter(site_code == "HO00")) +
    geom_point(aes(x = DateTime, y = spec_cond_uSpercm)) +
    geom_vline(xintercept = as.POSIXct(as.Date("2004-06-05")), color = "red") +
    labs(x = "Date", y = "Conductivity (uS/cm)", title = "HO00") +
    theme_bw())

(GV01_plot <- ggplot(chem_oursites_cond %>% filter(site_code == "GV01")) +
    geom_point(aes(x = DateTime, y = spec_cond_uSpercm)) +
    geom_vline(xintercept = as.POSIXct(as.Date("2004-06-05")), color = "red") +
    labs(x = "Date", y = "Conductivity (uS/cm)", title = "GV01") +
    theme_bw())

library(patchwork)

(all_sites <- RS02_plot / MC06_plot / AB00_plot / AT07_plot / SP02_plot / RG01_plot / HO00_plot / GV01_plot)

ggsave(all_sites,
       filename = "figures/crass_sbc_sites_cond_fire.png",
       width = 15,
       height = 40,
       units = "cm"
       )

(combined_sites <- ggplot(chem_oursites_cond) +
  geom_line(aes(x = DateTime, y = spec_cond_uSpercm, color = site_code)) +
  labs(x = "Date", y = "Conductivity (uS/cm)", title = "GV01") +
  theme_bw())

# End of script.
