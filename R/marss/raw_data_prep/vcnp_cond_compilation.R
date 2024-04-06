# Valles Caldera National Preserve (vcnp)
# Specific Conductivity Data Compilation
# August 11, 2023

#### Read Me ####
# The purpose of this script is to combine sonde and grab data from the VCNP. 

#### Setup ####

# Load libraries
library(here)
library(lubridate)
library(tidyverse)
library(calecopal)

# Load data
# Sonde data.
sonde_dat <- readRDS("data_raw/VCNP_sonde.rds")

# Grab sample data.
grab_dat <- readRDS("data_working/VCNPchem_edited_120521.rds")
# This data was previously edited, but not trimmed or aggregated, in the script named
# "vcnp_chem_compilation.R" by Betsy.

#### QAQC sonde data ####

sonde_trim <- map_df(sonde_dat,
         ~as.data.frame(.x), .id="site_name") %>%
  # Filter for sites of interest.
  filter(site_name %in% c("EFJ", "RED", "RSA", "RSAW")) %>%
  # Select only conductivity data.
  select(site_name, datetime_NM, SPC_uS_Value) %>%
  # And add column to denote data source
  mutate(source = "Sonde") %>%
  # And rename to match the grab sample dataset
  rename(DateTime = datetime_NM,
         cond_uScm = SPC_uS_Value)

# Also, need to remove QAQC to remove outliers from sonde data.
# Determine mean and sd conductivity values for each site.
summary <- sonde_trim %>%
  group_by(site_name) %>%
  summarize(mean = mean(cond_uScm, na.rm = TRUE),
            sd = sd(cond_uScm, na.rm = TRUE)) %>%
  ungroup()

meanEFJ <- 96.51694
sdEFJ <- 37.80697
meanRED <- 61.84900
sdRED <- 37.91663
meanRSA <- 88.43305
sdRSA <- 20.54047
meanRSAW <- 117.55648
sdRSAW <- 24.41510

# And QAQC original dataset accordingly.
sonde_qaqc <- sonde_trim %>%
  # First, removing all values equal to or less than 25 (not reasonable)
  filter(cond_uScm > 25) %>%
  # And, create new column with edited data.
  mutate(cond_uScm_ed = case_when(site_name == "EFJ" & cond_uScm >= meanEFJ+(4*sdEFJ) ~ NA,
                                  site_name == "EFJ" & cond_uScm <= meanEFJ-(4*sdEFJ) ~ NA,
                                  site_name == "RED" & cond_uScm >= meanRED+(4*sdRED) ~ NA,
                                  site_name == "RED" & cond_uScm <= meanRED-(4*sdRED) ~ NA,
                                  site_name == "RSA" & cond_uScm >= meanRSA+(4*sdRSA) ~ NA,
                                  site_name == "RSA" & cond_uScm <= meanRSA-(4*sdRSA) ~ NA,
                                  site_name == "RSAW" & cond_uScm >= meanRSAW+(4*sdRSAW) ~ NA,
                                  site_name == "RSAW" & cond_uScm <= meanRSAW-(4*sdRSAW) ~ NA,
                                  TRUE ~ cond_uScm)) %>%
  select(site_name, DateTime, cond_uScm_ed, source) %>%
  rename(cond_uScm = cond_uScm_ed)

# Quick plot to examine outliers that may have been removed.
(fig_out <- ggplot(sonde_trim %>%
                     filter(site_name %in% c("EFJ", "RED", "RSAW")) %>%
                     # removing one outrageous measure to plot better
                     filter(cond_uScm < 1000), #%>%
                     #filter(DateTime > as_datetime("01-01-2011 00:00")) %>%
                     #filter(DateTime < as_datetime("01-01-2015 00:00")),
                   aes(x = DateTime, 
                       y = cond_uScm,
                       color = site_name)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = (meanEFJ+(4*sdEFJ)), color = "pink") +
  geom_hline(yintercept = (meanRED+(4*sdRED)), color = "lightgreen") +
  geom_hline(yintercept = (meanRSAW+(4*sdRSAW)), color = "cornflowerblue") +
  labs(x = "Date",
       y = "Sonde SC",
       color = "Stream") +
  theme_bw())

# Export.
# ggsave(fig_out,
#        filename = "figures/all_raw_sonde_w_outlierlines_031224.jpg",
#        width = 30,
#        height = 10,
#        units = "cm")

# Counting how many records each removes from the sites of interest
# to address reviewer comments.
sonde_no_out <- sonde_trim %>% 
  # Select for sites of interest.
  filter(site_name %in% c("EFJ", "RED", "RSAW")) %>% # 491,441 records
  # First, removing all values equal to or less than 25 (not reasonable)
  filter(cond_uScm > 25) %>% # 447,970 records (8.85% removed)
  # And, actually filter out outliers.
  mutate(cond_uScm_ed = case_when(site_name == "EFJ" & cond_uScm >= meanEFJ+(4*sdEFJ) ~ NA,
                                  site_name == "EFJ" & cond_uScm <= meanEFJ-(4*sdEFJ) ~ NA,
                                  site_name == "RED" & cond_uScm >= meanRED+(4*sdRED) ~ NA,
                                  site_name == "RED" & cond_uScm <= meanRED-(4*sdRED) ~ NA,
                                  site_name == "RSAW" & cond_uScm >= meanRSAW+(4*sdRSAW) ~ NA,
                                  site_name == "RSAW" & cond_uScm <= meanRSAW-(4*sdRSAW) ~ NA,
                                  TRUE ~ cond_uScm)) %>%
  drop_na(cond_uScm_ed) # 447,111 records (0.19% removed)

#### Trim sonde data ####

# First, we want to remove all dates for which we will not be using sonde data.
# So, we need to keep:
# (1) All EFJ sonde data.
# (2) RED data after fire (before fire sonde data looks bad, 
# i.e. drops to zero lots, likely due to sediment infill).
# (3) All RSA sonde data.
# (4) All RSAW sonde data (only available post-fire, 
# but picks up where grab samples left off).

# Remove pre-fire RED data (Thompson Ridge fire started 2013-05-31).
sonde_trimmed <- sonde_qaqc %>%
  mutate(cond_uScm_ed = case_when(site_name == "RED" & 
                                    DateTime < ymd_hms("2013-05-31 00:00:00") ~ NA,
                                  TRUE ~ cond_uScm)) %>%
  select(site_name, DateTime, cond_uScm_ed, source)

# Second, we want to select one weekly value for all available sonde data to summarize
# on a monthly basis alongside grab sample data.
sonde_weekly <- sonde_trimmed %>%
  # Add date delineations to help with grouping.
  mutate(year = year(DateTime),
    week = strftime(DateTime, format = "%V"))

sonde_select <- sonde_weekly %>%
  group_by(site_name, year, week) %>%
  # Select the 150th observation from each group - roughly midday Monday
  summarize(DateTime_select = DateTime[150],
    cond_uScm_select = cond_uScm_ed[150]) %>%
  ungroup() %>%
  select(site_name, DateTime_select,cond_uScm_select) %>%
  rename(DateTime = DateTime_select,
         cond_uScm = cond_uScm_select) %>%
  mutate(source = "sonde")

#### Trim grab sample data ####

# Filter dataset down to columns of interest.
grab_dat_qaqc <- grab_dat %>%
  mutate(site_name = case_when(site_code == "East Fork Jemez River" ~ "EFJ",
                               site_code == "Redondo Creek" ~ "RED",
                               site_code %in% c("San Antonio Creek - Toledo",
                                                "San Antonio Creek- Toledo",
                                                "San Antonio Creek -Toledo") ~ "RSA",
                               site_code == "San Antonio - West" ~ "RSAW")) %>%
  filter(site_name %in% c("EFJ", "RED", "RSA", "RSAW")) %>%
  select(site_name, DateTime, mean_cond_uScm) %>%
  rename(cond_uScm = mean_cond_uScm) %>%
  mutate(source = "grab")

# First, we want to remove all dates for which we will not be using grab data.
# So, we need to keep:
# (1) No EFJ grab data.
# (2) All RED grab data (pre- and immediately post-fire).
# (3) No RSA grab data.
# (4) All RSAW grab data (only available pre-fire).

grab_trimmed <- grab_dat_qaqc %>%
  filter(site_name %in% c("RED", "RSAW")) %>%
  # AND remove grab samples below 25 uS (not reasonable)
  filter(cond_uScm > 25)

#### Join sonde and grab sample data ####

both_sonde_grab <- full_join(sonde_select, grab_trimmed)

# Export "raw" (a.k.a. not monthly aggregated) data.
saveRDS(both_sonde_grab, "data_working/VCNP_conductivity_sonde_grab_091123.rds")

# Plot datasets to compare to one another.
(cond_fig <- ggplot(both_sonde_grab, 
                   aes(x = DateTime, y = cond_uScm, color = source)) +
    geom_point(alpha = 0.7) +
    geom_line(alpha = 0.6) +
    scale_color_manual(values = c("#0FB2D3", "#026779")) +
    labs(x = "Date-Time", y = "Spec. Conductivity (uS)") +
    theme_bw() +
    theme(legend.position = "none") +
    facet_grid(rows = vars(site_name)))

# Export.
# ggsave(cond_fig,
#        filename = "figures/all_grab_vs_sonde_091123.jpg",
#        width = 20,
#        height = 20,
#        units = "cm")

#### Monthly summary & export ####

monthly_sonde_grab <- both_sonde_grab %>%
  mutate(Year = year(DateTime),
         Month = month(DateTime)) %>%
  group_by(site_name, Year, Month) %>%
  summarize(mean_cond_uScm = mean(cond_uScm, na.rm = TRUE)) %>%
  ungroup()

# Export aggregated data.
saveRDS(monthly_sonde_grab, "data_working/VCNP_monthly_conductivity_sonde_grab_091123.rds")

# End of script.
