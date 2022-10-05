### Merge Q with chems ###
## Currently working with a subset of chems (USGS_select_chems.csv)

## Inputs:
      # USGS_select_chems.csv
      # catchments_daily_discharge: output of query generated in script 02

library(here)
library(tidyverse)

### Read in data ###
select.chems <- read.csv(here("USGS_data", "USGS_select_chems.csv"))

### Merge Q w/ chems ###
## Prep data for merge
# Q is mean daily. chems are by date time. Will need to modify later if we want to examine sensor data.
select.chems$Date_UTC <- as.Date(select.chems$datetimeUTC, format = "%Y-%m-%d", tz = "UTC")

# Daily means
select.chems <- select.chems %>% group_by(usgs_site, CharacteristicName, Date_UTC) %>% 
                         summarize(mn_value_std = mean(value_std))
  
# Assuming Q Date is in UTC. Check w/ Stevan.
catchments_daily_discharge$Date_UTC <- as.Date(catchments_daily_discharge$Date, format = "%Y-%m-%d", tz = "UTC")

# Merge Q & chems
chemQ <- left_join(select.chems, catchments_daily_discharge, by = c("usgs_site", "Date_UTC"))

### Plots ###
## Discharge time series
chemQ %>% filter(grepl("092", usgs_site)) %>%
  ggplot(aes(x = Date_UTC, y = Flow)) +
  geom_point() +
  ylab("discharge (cfs)") +
  facet_wrap(~usgs_site, scales = "free_y")

## C-Q
chemQ %>% filter(CharacteristicName == "Nitrate") %>%
  filter(grepl("092", usgs_site)) %>%
  ggplot(aes(x = Flow, y = mn_value_std)) +
  geom_point() +
  ylab("Nitrate") +
  xlab("Q_cfs") +
  facet_wrap(~usgs_site, scales = "free")

chemQ %>% filter(CharacteristicName == "NO3_NO2") %>%
  filter(grepl("092", usgs_site)) %>%
  ggplot(aes(x = Flow, y = mn_value_std)) +
  geom_point() +
  ylab("NO3_NO2") +
  xlab("Q_cfs") +
  facet_wrap(~usgs_site, scales = "free")

chemQ %>% filter(CharacteristicName == "SPC") %>%
  filter(grepl("092", usgs_site)) %>%
  ggplot(aes(x = Flow, y = mn_value_std)) +
  geom_point() +
  ylab("SPC") +
  xlab("Q_cfs") +
  facet_wrap(~usgs_site, scales = "free")

chemQ %>% filter(CharacteristicName == "Turbidity") %>%
  filter(grepl("092", usgs_site)) %>%
  ggplot(aes(x = Flow, y = mn_value_std)) +
  geom_point() +
  ylab("Turbidity") +
  xlab("Q_cfs") +
  facet_wrap(~usgs_site, scales = "free")

write.csv(chemQ, here("USGS_data", "USGS_chem_Q.csv"), row.names = FALSE)


