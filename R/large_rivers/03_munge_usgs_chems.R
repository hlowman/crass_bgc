### Munging USGS chemistry data ###
## Harmonize units

## Inputs:
    # catchments_water_chemistry output from database query in script 02

## Outputs:
    # Subset of constituents, with common units: USGS_select_chems.csv

library(tidyverse)
library(here)

### Visualize and pare down chems ###
chems <- catchments_water_chemistry

### TKH notes to check w/ Stevan:
# DOC = 00681
# TOC = 00680	
# TDN = 00602
# DON = 00607
# NH4 = 00608
# NO3 = 00618
# NO3+NO2 = 00631
# PO4 (as PO4) = 00653
# ortho-P (as PO4) = 00660
# phosphorus (as P) = 00666
# ortho-P (as P) = 00671
# ortho-P (as P) = 51288
# ortho-P + condensed P (as PO4) = 52315
# NH4 (as NH4) = 71846
# Nitrate (as NO3) = 71851
# Phosphorus (as PO4) = 71888
# Nitrate (as NO3) = 91003
# Nitrate (as N) = 99121
# Orthophosphate (as PO4) = 99122
# Nitrate (as N), this might be sensor = 99133
# Nitrate, this might be sensor = 99136
# Nitrate (as N), this might be sensor = 99137
# NO3 + NO2 (as N) = 99889
# Silica dissolved (as SiO2) 00955
# Silica total (as SiO2) 00956
# Silica recoverable unfiltered (as SiO2) 00954

# Rename chemical constituents
chems <- chems %>% mutate(CharacteristicName = ifelse(CharacteristicName == "Specific conductance", "SPC", 
                                                  ifelse(CharacteristicName == "Ammonia and ammonium", "NH4",
                                                    ifelse(CharacteristicName == "Organic Nitrogen", "DON",
                                                      ifelse(CharacteristicName == "Total suspended solids", "TSS",
                                                        ifelse(CharacteristicName == "Inorganic carbon", "DIC",
                                                          ifelse(CharacteristicName == "Acidity, (H+)", "acidity",
                                                            ifelse(CharacteristicName == "Inorganic nitrogen (nitrate and nitrite and ammoni", "DIN",
                                                              ifelse(CharacteristicName == "Total dissolved solids", "TSS",
                                                                ifelse(CharacteristicName == "Nitrogen, mixed forms (NH3), (NH4), organic, (NO2) and (NO3)", "Nmisc",
                                                                  ifelse(CharacteristicName == "Inorganic nitrogen (nitrate and nitrite)", "NO3_NO2",
                                                                    ifelse(CharacteristicName == "Kjeldahl nitrogen", "KTDN",
                                                                      ifelse(CharacteristicName == "Total solids", "TS",
                                                                        ifelse(CharacteristicName == "UV 254", "abs254",
                                                                          ifelse(CharacteristicName == "Organic phosphorus", "orgP",
                                                                                 CharacteristicName)))))))))))))))
                          

# Working with ActivityStartDateTime as sampledatetime
chems <- chems %>% mutate(datetimeUTC = as.POSIXct(ActivityStartDateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

### Harmonize units ###
## Note: working with a subset of constituents to generate preliminary figures
print(n = 100, chems %>% group_by(CharacteristicName, ResultMeasure.MeasureUnitCode) %>% summarize(count = n()))

# column for converted units
chems$units_std <- NA

# column for values with common units
chems$value_std <- NA

## SPC
# select only uS/cm @25C. No way to find reference temperature for remaining records.
chems.SPC <- chems %>% filter(CharacteristicName == "SPC" & ResultMeasure.MeasureUnitCode == "uS/cm @25C") %>%
                       mutate(units_std = "us_cm_25C") %>%
                       mutate(value_std = ResultMeasureValue)

## Nitrate
# mg/l as N
# mg/l asNO3
#ug/l...assume NO3? only 33 instances
chems.NO3 <- chems %>% filter(CharacteristicName == "Nitrate") %>%
                       mutate(units_std = "uM") %>%
                       mutate(value_std = ifelse(ResultMeasure.MeasureUnitCode == "mg/l as N", ResultMeasureValue*(1000/14),
                                            ifelse(ResultMeasure.MeasureUnitCode == "mg/l asNO3", ResultMeasureValue*(1000/63), ResultMeasureValue/63)))

## NO3_NO2
# mg/l as N
chems.NO3_NO2 <- chems %>% filter(CharacteristicName == "NO3_NO2") %>%
                           mutate(units_std = "uM") %>%
                           mutate(value_std = ifelse(ResultMeasure.MeasureUnitCode == "mg/l as N", ResultMeasureValue*(1000/14),NA))

## TSS
# mg/l
# tons/ac ft. 1 US ton = 9.072*10^8 mg, 1 ac ft = 1.233*10^6L. Conversion factor = 735.7664
# tons/day Can't do anything with this before merging with Q. A lot of the data in these units.
chems.TSS <- chems %>% filter(CharacteristicName == "TSS") %>%
  mutate(units_std = "mgL") %>%
  mutate(value_std = ifelse(ResultMeasure.MeasureUnitCode == "mg/l", ResultMeasureValue,
                        ifelse(ResultMeasure.MeasureUnitCode == "tons/ac ft", ResultMeasureValue*735.7664, NA)))

## Turbidity
# FNU
# JTU: don't use. Historical visual method.
# NTU
# no conversion btw FNU & NTU. Should check that a site's measurements remain consistent throughout the period of record analyzed.
# Check for systematic difference btw FNU & NTU
chems %>% filter(CharacteristicName == "Turbidity") %>%
          group_by(CharacteristicName, ResultMeasure.MeasureUnitCode) %>% summarize(mn = mean(ResultMeasureValue, na.rm = TRUE), max = max(ResultMeasureValue, na.rm = TRUE), min = min(ResultMeasureValue, na.rm = TRUE), count = n())
# NTU at least 2* greater... leave these separate for now

chems.Turb <- chems %>% filter(CharacteristicName == "Turbidity") %>%
  mutate(units_std = ifelse(ResultMeasure.MeasureUnitCode == "NTU", "NTU",
                        ifelse(ResultMeasure.MeasureUnitCode == "FNU", "FNU", NA))) %>%
  mutate(value_std = ifelse(ResultMeasure.MeasureUnitCode == "JTU", NA, ResultMeasureValue))

## Phosphate
# Orthophosphate:
# mg/l as P
# mg/l asPO4
# ug/L as P

# Phosphate-phosphorus:
# mg/l

# Phosphorus:
# mg/l as P
# mg/l PO4

## Carbon
# Is this DOC?
# mg/l

## Potassium
chems.K <- chems %>% filter(CharacteristicName == "Potassium") %>%
  mutate(units_std = "uM") %>%
  mutate(value_std = ifelse(ResultMeasure.MeasureUnitCode == "mg/l", ResultMeasureValue*(1000/39.0983),
                            ifelse(ResultMeasure.MeasureUnitCode == "ug/l", ResultMeasureValue/39.0983, NA)))

## Sulfate
chems.SO4 <- chems %>% filter(CharacteristicName == "Sulfate") %>%
  mutate(units_std = "uM") %>%
  mutate(value_std = ifelse(ResultMeasure.MeasureUnitCode == "mg/l", ResultMeasureValue*(1000/96.06),
                            ifelse(ResultMeasure.MeasureUnitCode == "ug/l", ResultMeasureValue/96.06, NA)))

## Calcium
chems.Ca <- chems %>% filter(CharacteristicName == "Calcium") %>%
  mutate(units_std = "uM") %>%
  mutate(value_std = ifelse(ResultMeasure.MeasureUnitCode == "mg/l", ResultMeasureValue*(1000/40.078), NA))

select.chems <- rbind(chems.SPC, chems.NO3, chems.NO3_NO2, chems.TSS, chems.Turb, chems.K, chems.SO4, chems.Ca)

## How much data are there? 
# This is just a snapshot of sites with IDs beginning 092
select.chems %>% filter(CharacteristicName == "SPC")%>%
              filter(grepl("092", usgs_site)) %>%
          ggplot(aes(x = datetimeUTC, y = value_std)) +
            geom_point() +
            ylab("SPC") +
            facet_wrap(~usgs_site, scales = "free_y")
                 
select.chems %>% filter(CharacteristicName == "Nitrate") %>%
  filter(grepl("092", usgs_site)) %>%
  ggplot(aes(x = datetimeUTC, y = value_std)) +
  geom_point() +
  ylab("Nitrate") +
  facet_wrap(~usgs_site, scales = "free_y")

select.chems %>% filter(CharacteristicName == "NO3_NO2") %>%
  filter(grepl("092", usgs_site)) %>%
  ggplot(aes(x = datetimeUTC, y = value_std)) +
  geom_point() +
  ylab("NO3_NO2") +
  facet_wrap(~usgs_site, scales = "free_y")
# More SPC than nitrate

select.chems %>% filter(CharacteristicName == "TSS") %>%
  filter(grepl("092", usgs_site)) %>%
  ggplot(aes(x = datetimeUTC, y = value_std)) +
  geom_point() +
  ylab("TSS") +
  facet_wrap(~usgs_site, scales = "free_y")

select.chems %>% filter(CharacteristicName == "Turbidity") %>%
  filter(grepl("092", usgs_site)) %>%
  ggplot(aes(x = datetimeUTC, y = value_std)) +
  geom_point() +
  ylab("Turbidity") +
  facet_wrap(~usgs_site, scales = "free_y")
# much less turbidity data than solutes

# filter out catchments w/ all NA for any of the selected constituents in select.chems
select.chems <- select.chems %>% group_by(usgs_site, CharacteristicName) %>%
                                 filter(all(!is.na(value_std)))

dir.create(here("USGS_data"))

write.csv(select.chems, here("USGS_data", "USGS_select_chems.csv"), row.names = FALSE)

## Check whether sensor data are in the USGS chems database here

# CA nitrate sensor site: 0381142122015801
select.chems %>% filter(CharacteristicName == "NO3_NO2") %>%
  filter(grepl("0381142122015801", usgs_site)) %>%
  ggplot(aes(x = datetimeUTC, y = value_std)) +
  geom_point() +
  ylab("Nitrate") +
  facet_wrap(~usgs_site, scales = "free_y")

# AZ turbidity sensor site: 09497830
select.chems %>% filter(CharacteristicName == "Turbidity") %>%
  filter(grepl("09497830", usgs_site)) %>%
  ggplot(aes(x = datetimeUTC, y = value_std)) +
  geom_point() +
  ylab("Turbidity")

# AZ SPC sensor site: 09380000
select.chems %>% filter(CharacteristicName == "SPC") %>%
  filter(grepl("09380000", usgs_site)) %>%
  ggplot(aes(x = datetimeUTC, y = value_std)) +
  geom_point() +
  ylab("SPC")
# None of these are present in the dataset