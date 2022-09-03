### Exploratory plotting of USGS database ###
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
                                                                ifelse(CharacteristicName == "Nitrogen, mixed forms (NH3), (NH4), organic, (NO2) and (NO3)", "Tmix",
                                                                  ifelse(CharacteristicName == "Inorganic nitrogen (nitrate and nitrite)", "NO3_NO2",
                                                                    ifelse(CharacteristicName == "Kjeldahl nitrogen", "KTDN",
                                                                      ifelse(CharacteristicName == "Total solids", "TS",
                                                                        ifelse(CharacteristicName == "UV 254", "abs254",
                                                                          ifelse(CharacteristicName == "Organic phosphorus", "orgP",
                                                                                 CharacteristicName)))))))))))))))
                          

# Working with ActivityStartDateTime as sampledatetime
chems <- chems %>% mutate(ActivityStartDateTime = as.POSIXct(ActivityStartDateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

### Harmonize units ###

# new col for unit-converted values "conc_conv" & new column for harmonized units "unit_conv"
# Everything to uM
# select only SPC (@25)
# check into TSS

# 
chems %>% filter(CharacteristicName == "SPC" & ResultMeasure.MeasureUnitCode == "uS/cm @25C") %>%
              filter(grepl("092", usgs_site)) %>%
          ggplot(aes(x = ActivityStartDateTime, y = ResultMeasureValue)) +
            geom_point() +
            ylab("SPC") +
            facet_wrap(~usgs_site, scales = "free_y")
                 
chems %>% filter(CharacteristicName == "Nitrate") %>%
  filter(grepl("092", usgs_site)) %>%
  ggplot(aes(x = ActivityStartDateTime, y = ResultMeasureValue)) +
  geom_point() +
  ylab("Nitrate") +
  facet_wrap(~usgs_site, scales = "free_y")

### Pull Q for all catchments that have any chems ###

### Merge w/ fire ###
## Plot C-Q

## Now add LULC filter & replot C-Q


