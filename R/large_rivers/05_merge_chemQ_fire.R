### Merge chem, Q, & fire ###
## Calculate % burn for each catchment*fire
## Identify largest fire
## Identify most recent fire
## Categorize dates to pre/post fire, pre/post largest fire, pre/post most recent fire

## Inputs: 
      # Merged chems & Q: USGS_chem_Q.csv
      # Output of catchments_fires from script 02: catchment area and area burned for each fire

library(here)
library(tidyverse)

### Read in chemQ
chemQ <- read.csv(here("USGS_data", "USGS_chem_Q.csv"))
  
### Fire summary attributes ###
# Is Ig_Date local time or UTC? Assuming UTC here.
catchments_fires$Ig_Date <- as.Date(catchments_fires$Ig_Date, format = "%Y-%m-%d", tz = "UTC")

# Fraction catchment burned in each fire
catchments_fires <- catchments_fires %>% mutate(burn_ext = ws_burn_area_km2/catchment_area) 

# First handle singletons. The mutate summary commands will drop sites with a single observation, so first singletons (= only burned in one year) are found and rct & lg columns are created for them
singletons <- catchments_fires %>% group_by(usgs_site) %>%
  filter(n() == 1) %>%
  mutate(Date_rct = Ig_Date) %>%
  mutate(pctburn_rct = burn_ext) %>%
  mutate(areaburn_rct = ws_burn_area_km2) %>%
  mutate(Date_lg = Ig_Date) %>%
  mutate(pctburn_lg = burn_ext) %>%
  mutate(areaburn_lg = ws_burn_area_km2) %>%
  mutate(Event_ID_lg = Event_ID) %>%
  mutate(Event_ID_rct = Event_ID) %>%
  mutate(Incid_Name_lg = Incid_Name) %>%
  mutate(Incid_Name_rct = Incid_Name) %>%
  select(-c(Ig_Date, ws_burn_area_km2, burn_ext, total_burn_area_km2, Event_ID, Incid_Name))
# 527 catchments one fires

# Subset most recent fire each catchment
recent <- catchments_fires %>% group_by(usgs_site) %>%
  filter(n() > 1) %>%
  filter(Ig_Date == max(Ig_Date)) %>%
  select(-c(total_burn_area_km2))
# 499 catchments burned 2 or more times

# Subset largest fire each catchment
largest <- catchments_fires %>% group_by(usgs_site) %>%
  filter(n() > 1) %>%
  filter(burn_ext == max(burn_ext)) %>%
  select(-c(total_burn_area_km2))
# 4 missing catchments here... must be due to NAs in burn_ext

#rename to reflect max and recent burn years & burn areas
names(recent) <- c("usgs_site", "catchment_area", "Event_ID_rct", "Incid_Name_rct", "Date_rct", "areaburn_rct", "pctburn_rct")

names(largest) <- c("usgs_site", "catchment_area", "Event_ID_lg", "Incid_Name_lg", "Date_lg", "areaburn_lg", "pctburn_lg")

# Join most recent and largest
fires <- merge(recent, largest, by = c("usgs_site", "catchment_area"))
fires <- bind_rows(singletons, fires)

# Join to catchments_summary for cumulative burn extent
fires <- left_join(fires, catchments_summary)

# re-format dates

# Set fire year to 1980 for unburned catchments. This is the temporal extent of the MTBS fire data source
# note: using dplyr if_else because base ifelse strips off formatting
fires <- fires %>% ungroup() %>%
                   mutate_at(c("Date_rct", "Date_lg"), ~if_else(is.na(Date_rct), as.Date("1979-12-31"), .)) 

fires %>% select("usgs_site") %>% n_distinct()
# 1022 catchments

chemQ %>% select("usgs_site") %>% n_distinct()
# 1016. Not all catchments have data for the particular chems selected

fires %>% ggplot(aes(x = pctburn_lg)) +
          geom_histogram(color = "red") +
          geom_histogram(aes(x = pctburn_rct), color = "blue")
  
fires %>% ggplot(aes(x = pctburn_lg)) +
  geom_histogram(color = "red")

fires %>% ggplot(aes(x = pctburn_rct)) +
  geom_histogram(color = "blue")

fires %>% group_by(usgs_site) %>% filter(!is.na(pctburn_lg)) %>% nrow()
# 651 catchments with fires... seems like a lot?

fires %>% filter(is.na(pctburn_lg)) %>% nrow()
# 375... seems like too few? ratio of burned:unburned catchments seems skewed

write.csv(fires, here("USGS_data", "USGS_spatial_fires.csv"), row.names = FALSE)

### Merge fire data with chemQ ###
chemQsp <- left_join(chemQ, fires, by = "usgs_site")

# Quick plot to check
chemQsp %>% filter(!is.na(Incid_Name_lg) & CharacteristicName == "Nitrate" & Date_UTC > "2021-05-01") %>%
            ggplot(aes(x = pctburn_lg, y = mn_value_std)) +
            geom_point()

### Summary metrics ###
## Time since fire in fractional years
# recent fire

# largest fire

# pre/post burn column

### Plot C-Q pre/post burn ###
# cumulative fire vs. recent mean chem

######
# C-Q pre-post burn for burned catchments (how long post-burn window?): same year, 1 y, 2, y, 3 y, etc.
# time series for burned catchments

### Next scripts
## Then apply filters: catchment size, Now add LULC filter & replot C-Q