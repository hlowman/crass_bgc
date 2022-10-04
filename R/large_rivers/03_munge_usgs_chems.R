### Exploratory plotting of USGS database ###

## 
## Join Q, chem, & spatial data queried from firearea_db in 02_query_firearea_database.R

library(tidyverse)
library(here)

### Visualize and pare down chems ###
chems <- catchments_water_chemistry

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

select.chems <- rbind(chems.SPC, chems.NO3, chems.NO3_NO2, chems.TSS, chems.Turb)

## How much data are there? 
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
  ylab("TSS") +
  facet_wrap(~usgs_site, scales = "free_y")
# much less turbidity

dir.create(here("USGS_data"))

write.csv(select.chems, here("USGS_data", "USGS_select_chems.csv"), row.names = FALSE)

### Next scripts:
### Pull Q for all catchments that have any chems ###
# filter out catchments w/ NA for all constituents in select.chems

### Merge w/ fire ###
# fire: merge w/ catchment chars and calc % catchment burned

## Plot C-Q

## Now add LULC filter & replot C-Q


