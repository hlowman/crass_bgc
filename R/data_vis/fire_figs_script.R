# Time series figures
# 6/16/2021
# Heili Lowman

# This script will create some of the initial figures, including fire start
# dates to examine different constituents by month and watershed through
# time.

#### Setup ####
# Load packages
library(tidyverse)
library(lubridate)
library(naniar)
library(viridis)

# Load data
dat <- read_csv("data/sbclter_stream_chemistry_allyears_registered_stations_20190628.csv") %>% # Load in LTER dataset
  replace_with_na_all(condition = ~.x == -999) # Replace all -999s with NAs.

# Pivot data
dat_long <- dat %>%
  pivot_longer(
    cols = nh4_uM:spec_cond_uSpercm,
    names_to = "analyte",
    values_to = "value")

# First, I need add in each watershed name so that I can assign the correct
# fire occurrences to each watershed.

dat_shed <- dat_long %>%
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
    site_code == "MC00" ~ "Mission Creek",
    site_code == "MC06" ~ "Mission Creek",
    site_code == "RS02" ~ "Mission Creek")))
# Used delineation at: https://databasin.org/maps/new/#datasets=6ad26ddb04ae4dbc9362303628270daf

# Now, also here are the start dates and affected watersheds
# for each of the fires.

# Gaviota (6/4/2004) - Canada De Santa Anita, Canada De La Gaviota, Tajiguas Creek
# Sherpa (6/15/2016) - Tajiguas Creek, Dos Pueblos Canyon
# Whittier (7/7/2017) - Tajiguas Creek, Dos Pueblos Canyon
# Gap (7/1/2008) - Dos Pueblos Canyon, San Pedro Creek
# Cave (11/25/2019) - Atascadero Creek
# Jesusita (5/5/2009) - Atascadero Creek, Mission Creek
# Tea (11/13/2008) - Mission Creek
# Thomas (12/4/2017) - Mission Creek

# Used delineation at: https://sbc-gis.maps.arcgis.com/apps/dashboards/316aef44f80743c88ee1af95ab2f64ed

# Making a series of 7 paneled boxplots, one for each watershed,
# with wildfire start dates denoted:

#### Canada De Santa Anita ####
fig1 <- dat_shed %>%
  # filter for first watershed
  filter(watershed == "Canada De Santa Anita") %>%
  # round all dates to the appropriate month & create new column for year
  mutate(month_date = floor_date(timestamp_local, unit = "month"),
         year = year(timestamp_local)) %>%
  # create paneled boxplot - change to yearly
  ggplot(aes(x = year, y = value, group = year)) +
  geom_boxplot(aes(color = analyte)) +
  geom_vline(xintercept = 2004, 
             color = "red") + # Gaviota fire date
  xlim(2003, 2010) + # set date limits so as to include the fire date
  scale_color_viridis(discrete = TRUE) + # color by analyte
  facet_grid(rows = vars(analyte), scales = "free") + # facet by analyte
  theme_bw() + # remove background
  theme(axis.title.y = element_blank(), # remove y axis title
        legend.position = "none") + # remove legend
  labs(x = "Date", 
       title = "Canada De Santa Anita Watershed",
       subtitle = "San Onofre Creek")

# Export figure
# ggsave("figures/crass_watershed1_yr.png",
#        fig1,
#        height = 8, width = 10,
#        units = "in")

#### Canada De La Gaviota ####
fig2 <- dat_shed %>%
  # filter for next watershed
  filter(watershed == "Canada De La Gaviota") %>%
  # round all dates to the appropriate month 
  mutate(month_date = floor_date(timestamp_local, unit = "month"),
         year = year(timestamp_local)) %>%
  # create paneled boxplot
  ggplot(aes(x = year, y = value, group = year)) +
  geom_boxplot(aes(color = analyte)) +
  geom_vline(xintercept = 2004, 
             color = "red") + # Gaviota fire date
  scale_color_viridis(discrete = TRUE) + # color by analyte
  facet_grid(rows = vars(analyte), scales = "free") + # facet by analyte
  theme_bw() + # remove background
  theme(axis.title.y = element_blank(), # remove y axis title
        legend.position = "none") + # remove legend
  labs(x = "Date", 
       title = "Canada De La Gaviota Watershed",
       subtitle = "Gaviota Creek")

# Export figure
# ggsave("figures/crass_watershed2_yr.png",
#        fig2,
#        height = 8, width = 10,
#        units = "in")

#### Tajiguas Creek ####
fig3 <- dat_shed %>%
  # filter for next watershed
  filter(watershed == "Tajiguas Creek") %>%
  mutate(month_date = floor_date(timestamp_local, unit = "month"),
         year = year(timestamp_local)) %>%
  # create paneled boxplot
  ggplot(aes(x = year, y = value, group = year)) +
  geom_boxplot(aes(color = analyte)) +
  geom_vline(xintercept = 2004, 
             color = "red") + # Gaviota fire date
  geom_vline(xintercept = 2016, 
             color = "red") + # Sherpa fire date
  geom_vline(xintercept = 2017, 
             color = "red") + # Whittier fire date
  scale_color_viridis(discrete = TRUE) + # color by analyte
  facet_grid(rows = vars(analyte), scales = "free") + # facet by analyte
  theme_bw() + # remove background
  theme(axis.title.y = element_blank(), # remove y axis title
        legend.position = "none") + # remove legend
  labs(x = "Date", 
       title = "Tajiguas Creek Watershed",
       subtitle = "Arroyo Hondo and Refugio Creeks")

# Export figure
# ggsave("figures/crass_watershed3_yr.png",
#        fig3,
#        height = 8, width = 10,
#        units = "in")

#### Dos Pueblos Canyon ####
fig4 <- dat_shed %>%
  # filter for next watershed
  filter(watershed == "Dos Pueblos Canyon") %>%
  mutate(month_date = floor_date(timestamp_local, unit = "month"),
         year = year(timestamp_local)) %>%
  # create paneled boxplot
  ggplot(aes(x = year, y = value, group = year)) +
  geom_boxplot(aes(color = analyte)) +
  geom_vline(xintercept = 2008, 
             color = "red") + # Gap fire date
  geom_vline(xintercept = 2016, 
             color = "red") + # Sherpa fire date
  geom_vline(xintercept = 2017, 
             color = "red") + # Whittier fire date
  scale_color_viridis(discrete = TRUE) + # color by analyte
  facet_grid(rows = vars(analyte), scales = "free") + # facet by analyte
  theme_bw() + # remove background
  theme(axis.title.y = element_blank(), # remove y axis title
        legend.position = "none") + # remove legend
  labs(x = "Date", 
       title = "Dos Pueblos Canyon Watershed",
       subtitle = "Tecolote, Bell Canyon, and Devereaux Creeks")

# Export figure
# ggsave("figures/crass_watershed4_yr.png",
#        fig4,
#        height = 8, width = 10,
#        units = "in")

#### San Pedro Creek ####
fig5 <- dat_shed %>%
  # filter for next watershed
  filter(watershed == "San Pedro Creek") %>%
  mutate(month_date = floor_date(timestamp_local, unit = "month"),
         year = year(timestamp_local)) %>%
  # create paneled boxplot
  ggplot(aes(x = year, y = value, group = year)) +
  geom_boxplot(aes(color = analyte)) +
  geom_vline(xintercept = 2008, 
             color = "red") + # Gap fire date
  scale_color_viridis(discrete = TRUE) + # color by analyte
  facet_grid(rows = vars(analyte), scales = "free") + # facet by analyte
  theme_bw() + # remove background
  theme(axis.title.y = element_blank(), # remove y axis title
        legend.position = "none") + # remove legend
  labs(x = "Date", 
       title = "San Pedro Creek Watershed",
       subtitle = "San Pedro Creek")

# Export figure
# ggsave("figures/crass_watershed5_yr.png",
#        fig5,
#        height = 8, width = 10,
#        units = "in")

#### Atascadero Creek ####
fig6 <- dat_shed %>%
  # filter for next watershed
  filter(watershed == "Atascadero Creek") %>%
  mutate(month_date = floor_date(timestamp_local, unit = "month"),
         year = year(timestamp_local)) %>%
  # create paneled boxplot
  ggplot(aes(x = year, y = value, group = year)) +
  geom_boxplot(aes(color = analyte)) +
  geom_vline(xintercept = 2009, 
             color = "red") + # Jesusita fire date
  geom_vline(xintercept = 2019, 
             color = "red") + # Cave fire date
  scale_color_viridis(discrete = TRUE) + # color by analyte
  facet_grid(rows = vars(analyte), scales = "free") + # facet by analyte
  theme_bw() + # remove background
  theme(axis.title.y = element_blank(), # remove y axis title
        legend.position = "none") + # remove legend
  labs(x = "Date", 
       title = "Atascadero Creek Watershed",
       subtitle = "Atascadero Creek")

# Export figure
# ggsave("figures/crass_watershed6_yr.png",
#        fig6,
#        height = 8, width = 10,
#        units = "in")

#### Mission Creek ####
fig7 <- dat_shed %>%
  # filter for next watershed
  filter(watershed == "Mission Creek") %>%
  mutate(month_date = floor_date(timestamp_local, unit = "month"),
         year = year(timestamp_local)) %>%
  # create paneled boxplot
  ggplot(aes(x = year, y = value, group = year)) +
  geom_boxplot(aes(color = analyte)) +
  geom_vline(xintercept = 2009, 
             color = "red") + # Jesusita fire date
  geom_vline(xintercept = 2008, 
             color = "red") + # Tea fire date
  geom_vline(xintercept = 2017, 
             color = "red") + # Thomas fire date
  scale_color_viridis(discrete = TRUE) + # color by analyte
  facet_grid(rows = vars(analyte), scales = "free") + # facet by analyte
  theme_bw() + # remove background
  theme(axis.title.y = element_blank(), # remove y axis title
        legend.position = "none") + # remove legend
  labs(x = "Date", 
       title = "Mission Creek Watershed",
       subtitle = "Arroyo Burro, Rattlesnake, and Mission Creeks")

# Export figure
# ggsave("figures/crass_watershed7_yr.png",
#        fig7,
#        height = 8, width = 10,
#        units = "in")

# End of script.
