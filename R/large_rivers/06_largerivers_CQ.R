### Pre/post fire C-Q relationships for large rivers ###
## Demonstration of fire effects in large rivers: figure for conceptual paper

library(here)
library(tidyverse)

## Data ##
chemQfire <- read.csv(here("USGS_data", "chemQfire.csv"))

# format date
chemQfire$Date_UTC <- as.Date(chemQfire$Date_UTC, format = "%Y-%m-%d")
  
########################
### C-Q large rivers ###
########################
large <- chemQfire %>% filter(catchment_area > 18000)
# 78 sites

# How many sites with sufficient data for C-Q analysis?
large %>% group_by(usgs_site, post_lg) %>% 
          filter(CharacteristicName == "SPC" & 
                 !is.na(mn_value_std) & 
                 !is.na(post_lg) &
                 n() > 10) %>%
          ungroup() %>%
          distinct(usgs_site)
# 73

large %>% group_by(usgs_site, post_lg) %>% 
  filter(CharacteristicName == "Turbidity" & 
           !is.na(mn_value_std) & 
           !is.na(post_lg) &
           n() > 10) %>%
  ungroup() %>%
  distinct(usgs_site)
# 42

large %>% group_by(usgs_site, post_lg) %>% 
  filter(CharacteristicName == "Nitrate" & 
           !is.na(mn_value_std) & 
           !is.na(post_lg) &
           n() > 10) %>%
  ungroup() %>%
  distinct(usgs_site)
# 53

large %>% group_by(usgs_site, post_lg) %>% 
  filter(CharacteristicName == "NO3_NO2" & 
           !is.na(mn_value_std) & 
           !is.na(post_lg) &
           n() > 10) %>%
  ungroup() %>%
  distinct(usgs_site)
# 45

## Visualize pre/post C-Q ##
## Identify catchments with sufficient post-fire observations.
large %>% group_by(usgs_site, CharacteristicName, post_lg) %>% 
                      filter(!is.na(mn_value_std) & 
                             !is.na(Flow) &
                             !is.na(post_lg) &
                             #yrs_lg > -5 & 
                             yrs_lg <= 3 &
                             n() >= 8 &
                             CharacteristicName == "NO3_NO2") %>%
  ggplot(aes(x = log(Flow), y = log(mn_value_std), group = post_lg)) +
  geom_point(aes(color = post_lg)) +
  facet_wrap(~usgs_site, scales = "free")

large %>% group_by(usgs_site, CharacteristicName, post_lg) %>% 
  filter(!is.na(mn_value_std) & 
           !is.na(Flow) &
           !is.na(post_rct) &
           #yrs_recent > -5 & 
           yrs_recent <= 3 &
           n() >= 8 &
           CharacteristicName == "Turbidity") %>%
  ggplot(aes(x = log(Flow), y = log(mn_value_std), group = post_lg)) +
  geom_point(aes(color = post_lg)) +
  facet_wrap(~usgs_site, scales = "free")

large %>% group_by(usgs_site, CharacteristicName, post_lg) %>% 
  filter(!is.na(mn_value_std) & 
           !is.na(Flow) &
           !is.na(post_rct) &
           #yrs_recent > -5 & 
           yrs_recent <= 3 &
           n() >= 8 &
           CharacteristicName == "Potassium") %>%
  ggplot(aes(x = log(Flow), y = log(mn_value_std), group = post_lg)) +
  geom_point(aes(color = post_lg)) +
  facet_wrap(~usgs_site, scales = "free")

#################
## Time series ##
#################
large %>% group_by(usgs_site, CharacteristicName, post_lg) %>% 
  filter(!is.na(mn_value_std) & 
           #!is.na(Flow) &
           !is.na(post_lg) &
           n() > 10 &
           CharacteristicName == "NO3_NO2") %>%
  ggplot(aes(x = Date_UTC, y = mn_value_std, group = post_lg)) +
  geom_point(aes(color = post_lg)) +
  facet_wrap(~usgs_site, scales = "free")





### Plot C-Q pre/post burn ###
# cumulative fire vs. recent mean chem

## Then apply filters: catchment size, Now add LULC filter & replot C-Q

######
# C-Q pre-post burn for burned catchments (how long post-burn window?): same year, 1 y, 2, y, 3 y, etc.
# time series for burned catchments


