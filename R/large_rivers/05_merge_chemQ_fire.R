### Merge chem, Q, & fire ###
## Calculate % burn for each catchment*fire
## Identify largest fire
## Identify most recent fire
## Categorize dates to pre/post fire, pre/post largest fire, pre/post most recent fire
## Calculate years since most recent and years since largest fire

## Inputs: 
      # Merged chems & Q: USGS_chem_Q.csv
      # Output of catchments_fires from script 02: catchment area and area burned for each fire

library(here)
library(tidyverse)
library(lubridate)

### Read in chemQ
chemQ <- read.csv(here("USGS_data", "USGS_chem_Q.csv"))
 
############################### 
### Fire summary attributes ###
###############################
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

nburn <- fires %>% group_by(usgs_site) %>% filter(!is.na(pctburn_lg)) %>% nrow()
# 651 catchments with fires... seems like a lot?

nunburn <- fires %>% filter(is.na(pctburn_lg)) %>% nrow()
# 375 unburned catchments... seems like too few? ratio of burned:unburned catchments seems skewed

pctburn_hist <- fires %>% ggplot(aes(x = pctburn_lg)) +
                          geom_histogram(color = "red") +
                          geom_histogram(aes(x = pctburn_rct), color = "blue") +
                          annotate("text", label = paste("N burned catchments =", nburn), x = 0.75, y = 400) +
                          annotate("text", label = paste("N unburned catchments =", nunburn), x = 0.75, y = 375) +
                          ylab("count") +
                          xlab("percent catchment burned") +
                          theme_bw() 

ggsave(pctburn_hist, path = here("USGS_data", "plots"), file = "pctburn_hist.pdf", width = 10, height = 10, units = "in")
  
fires %>% ggplot(aes(x = pctburn_lg)) +
  geom_histogram(color = "red")

fires %>% ggplot(aes(x = pctburn_rct)) +
  geom_histogram(color = "blue")

write.csv(fires, here("USGS_data", "USGS_spatial_fires.csv"), row.names = FALSE)

### Merge fire data with chemQ ###
chemQsp <- left_join(chemQ, fires, by = "usgs_site")

# Quick plot to check
chemQsp %>% filter(!is.na(Incid_Name_lg) & CharacteristicName == "Nitrate" & Date_UTC > "2021-05-01") %>%
            ggplot(aes(x = pctburn_lg, y = mn_value_std)) +
            geom_point()

## Summary metrics
## Time since most recent and largest fires in fractional years
# recent fire
chemQsp <- chemQsp %>% mutate(yrs_recent = time_length(difftime(as.POSIXct(Date_UTC), as.POSIXct(Date_rct)), "years")) %>%
                       mutate(yrs_lg = time_length(difftime(as.POSIXct(Date_UTC), as.POSIXct(Date_lg)), "years")) 

# pre/post burn column
chemQsp <- chemQsp %>% mutate(post_rct = ifelse(Date_UTC < Date_rct, "pre", "post")) %>%
                       mutate(post_lg = ifelse(Date_UTC < Date_lg, "pre", "post")) 

## Export ###
write.csv(chemQsp, here("USGS_data", "chemQfire.csv"), row.names = FALSE)

### Data summaries ###
# >/= 18 obs 3y pre-fire & >/= 18 obs 3y post-fire chem+Q
# N catchments
# list of site IDs
# recent fire
Nrct <- chemQsp %>% filter(!is.na(Flow)) %>%
                    count(usgs_site, CharacteristicName, post_rct) %>%
                    pivot_wider(names_from = post_rct, values_from = n) %>%
                    filter(post >= 18 & pre >= 18) %>%
                    pivot_wider(names_from = CharacteristicName, values_from = c(post, pre))

# largest fire
Nlg <- chemQsp %>% filter(!is.na(Flow)) %>%
                   count(usgs_site, CharacteristicName, post_lg) %>%
                   pivot_wider(names_from = post_lg, values_from = n) %>%
                   filter(post >= 18 & pre >= 18) %>%
                   pivot_wider(names_from = CharacteristicName, values_from = c(post, pre))

# Export
write.csv(Nrct, here("USGS_data", "data_summaries", "Qchem_18obs_prepost_rctfire.csv"), row.names = FALSE)
write.csv(Nlg, here("USGS_data", "data_summaries", "Qchem_18obs_prepost_lgfire.csv"), row.names = FALSE)

## Example time series figs
# pre/post largest fire & has Q
SPC.pl <- chemQsp %>% filter(usgs_site %in% Nlg$usgs_site) %>%
                             filter(CharacteristicName == "SPC") %>%
                             ggplot(aes(x = Date_UTC, y = mn_value_std, group = post_lg)) +
                              geom_point(aes(color = post_lg)) +
                              facet_wrap(~usgs_site, scales = "free") +
                              ylab("SPC")

ggsave(SPC.pl, path = here("USGS_data", "plots"), file = "SPC_lg.pdf", units = "in")

NlgNO3NO2 <- Nlg %>% filter(!is.na(pre_NO3_NO2) & !is.na(post_NO3_NO2))
NO3NO2.pl <- chemQsp %>% filter(usgs_site %in% NlgNO3NO2$usgs_site) %>%
                         filter(CharacteristicName == "NO3_NO2") %>%
                         ggplot(aes(x = Date_UTC, y = mn_value_std, group = post_lg)) +
                          geom_point(aes(color = post_lg)) +
                          facet_wrap(~usgs_site, scales = "free") +
                          ylab("NO3_NO2")

ggsave(NO3NO2.pl, path = here("USGS_data", "plots"), file = "NO3NO2_lg.pdf", units = "in")

NlgNitrate <- Nlg %>% filter(!is.na(pre_Nitrate) & !is.na(post_Nitrate))
Nitrate.pl <- chemQsp %>% filter(usgs_site %in% NlgNitrate$usgs_site) %>%
                          filter(CharacteristicName == "Nitrate") %>%
                          ggplot(aes(x = Date_UTC, y = mn_value_std, group = post_lg)) +
                            geom_point(aes(color = post_lg)) +
                            facet_wrap(~usgs_site, scales = "free") +
                            ylab("Nitrate")

ggsave(Nitrate.pl, path = here("USGS_data", "plots"), file = "Nitrate_lg.pdf", units = "in")

NlgK <- Nlg %>% filter(!is.na(pre_Potassium) & !is.na(post_Potassium))
K.pl <- chemQsp %>% filter(usgs_site %in% NlgK$usgs_site) %>%
  filter(CharacteristicName == "Potassium") %>%
  ggplot(aes(x = Date_UTC, y = mn_value_std, group = post_lg)) +
  geom_point(aes(color = post_lg)) +
  facet_wrap(~usgs_site, scales = "free") +
  ylab("Potassium")

ggsave(K.pl, path = here("USGS_data", "plots"), file = "K_lg.pdf", units = "in")

Nlgturb <- Nlg %>% filter(!is.na(pre_Turbidity) & !is.na(post_Turbidity))
Turb.pl <- chemQsp %>% filter(usgs_site %in% Nlgturb$usgs_site) %>%
  filter(CharacteristicName == "Turbidity") %>%
  ggplot(aes(x = Date_UTC, y = mn_value_std, group = post_lg)) +
  geom_point(aes(color = post_lg)) +
  facet_wrap(~usgs_site, scales = "free") +
  ylab("Turbidity")

ggsave(Turb.pl, path = here("USGS_data", "plots"), file = "Turb_lg.pdf", units = "in")

NlgSO4 <- Nlg %>% filter(!is.na(pre_Sulfate) & !is.na(post_Sulfate))
SO4.pl <- chemQsp %>% filter(usgs_site %in% NlgSO4$usgs_site) %>%
  filter(CharacteristicName == "Sulfate") %>%
  ggplot(aes(x = Date_UTC, y = mn_value_std, group = post_lg)) +
  geom_point(aes(color = post_lg)) +
  facet_wrap(~usgs_site, scales = "free") +
  ylab("Sulfate")

ggsave(SO4.pl, path = here("USGS_data", "plots"), file = "SO4_lg.pdf", units = "in")

NlgCa <- Nlg %>% filter(!is.na(pre_Calcium) & !is.na(post_Calcium))
Ca.pl <- chemQsp %>% filter(usgs_site %in% NlgCa$usgs_site) %>%
  filter(CharacteristicName == "Calcium") %>%
  ggplot(aes(x = Date_UTC, y = mn_value_std, group = post_lg)) +
  geom_point(aes(color = post_lg)) +
  facet_wrap(~usgs_site, scales = "free") +
  ylab("Calcium")

ggsave(Ca.pl, path = here("USGS_data", "plots"), file = "Ca_lg.pdf", units = "in")

# C-Q
Ca.CQ.pl <- chemQsp %>% filter(usgs_site %in% NlgCa$usgs_site) %>%
                        filter(CharacteristicName == "Calcium") %>%
                        ggplot(aes(x = log(Flow), y = log(mn_value_std), group = post_lg)) +
                        geom_point(aes(color = post_lg)) +
                        facet_wrap(~usgs_site, scales = "free") +
                        ylab("Calcium")

NO3_NO2.CQ.pl <- chemQsp %>% filter(usgs_site %in% NlgNO3NO2$usgs_site) %>%
  filter(CharacteristicName == "NO3_NO2") %>%
  ggplot(aes(x = log(Flow), y = log(mn_value_std), group = post_lg)) +
  geom_point(aes(color = post_lg)) +
  facet_wrap(~usgs_site, scales = "free") +
  ylab("NO3_NO2")

ggsave(NO3_NO2.CQ.pl, path = here("USGS_data", "plots"), file = "NO3NO2.CQ.pdf", width = 14, height = 12, units = "in")

# Scatter
nNO3_NO2 <- chemQsp %>% filter(!grepl("Cultivated Crops", dominant_lulc)) %>%
                        filter(CharacteristicName == "NO3_NO2")
pairs(~log(mn_value_std) + Flow + catchment_area + pctburn_rct + pctburn_lg + number_fires + total_area_burned + number_wwtp + yrs_recent + yrs_lg, data = nNO3_NO2, cols = as.factor(nNO3_NO2$post_lg))
