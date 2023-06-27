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
# 75

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
# 54

large %>% group_by(usgs_site, post_lg) %>% 
  filter(CharacteristicName == "NO3_NO2" & 
           !is.na(mn_value_std) & 
           !is.na(post_lg) &
           n() > 10) %>%
  ungroup() %>%
  distinct(usgs_site)
# 46


##############################
### Visualize pre/post C-Q ###
##############################

### Plot C-Q pre/post burn ###
# cumulative fire vs. recent mean chem

## Then apply filters: catchment size, Now add LULC filter & replot C-Q

######
# C-Q pre-post burn for burned catchments (how long post-burn window?): same year, 1 y, 2, y, 3 y, etc.
# time series for burned catchments

## Identify catchments with sufficient post-fire observations.
NO3NO2 <- large %>% group_by(usgs_site, CharacteristicName, post_lg) %>% 
                    filter(!is.na(mn_value_std) & 
                           !is.na(Flow) &
                           !is.na(post_lg) &
                           #yrs_lg > -5 & 
                           yrs_lg <= 3 &
                           n() >= 8 &
                           CharacteristicName == "NO3_NO2") %>%
                    ggplot(aes(x = log(Flow), y = log(mn_value_std), group = post_lg)) +
                      geom_point(aes(color = post_lg)) +
                      ylab("log(NO3_NO2)") +
                      facet_wrap(~usgs_site, scales = "free")

ggsave(NO3NO2, path = here("USGS_data", "plots"), file = "lgCQ_NO3NO2.pdf", width = 15, height = 14, units = "in")

turb <- large %>% group_by(usgs_site, CharacteristicName, post_lg) %>% 
                  filter(!is.na(mn_value_std) & 
                         !is.na(Flow) &
                         !is.na(post_rct) &
                         #yrs_recent > -5 & 
                         yrs_recent <= 3 &
                         n() >= 8 &
                         CharacteristicName == "Turbidity") %>%
                  ggplot(aes(x = log(Flow), y = log(mn_value_std), group = post_lg)) +
                    geom_point(aes(color = post_lg)) +
                    ylab("log(Turbidity)") +
                    facet_wrap(~usgs_site, scales = "free")

ggsave(turb, path = here("USGS_data", "plots"), file = "lgCQ_turb.pdf", width = 16, height = 14, units = "in")

Kp <- large %>% group_by(usgs_site, CharacteristicName, post_lg) %>% 
                filter(!is.na(mn_value_std) & 
                       !is.na(Flow) &
                       !is.na(post_rct) &
                       #yrs_recent > -5 & 
                       yrs_recent <= 3 &
                       n() >= 8 &
                       CharacteristicName == "Potassium") %>%
                ggplot(aes(x = log(Flow), y = log(mn_value_std), group = post_lg)) +
                  geom_point(aes(color = post_lg)) +
                  ylab("log(Potassium") +
                  facet_wrap(~usgs_site, scales = "free")

ggsave(Kp, path = here("USGS_data", "plots"), file = "lgCQ_K.pdf", width = 16, height = 14, units = "in")

large %>% group_by(usgs_site, CharacteristicName, post_lg) %>% 
  filter(!is.na(mn_value_std) & 
           !is.na(Flow) &
           !is.na(post_rct) &
           #yrs_recent > -5 & 
           yrs_recent <= 3 &
           n() >= 8 &
           CharacteristicName == "SPC") %>%
  ggplot(aes(x = log(Flow), y = log(mn_value_std), group = post_lg)) +
  geom_point(aes(color = post_lg)) +
  facet_wrap(~usgs_site, scales = "free")

## Export data for Bayesian regression test
NCQ <- chemQfire %>% group_by(usgs_site, CharacteristicName, post_lg) %>% 
  filter(!is.na(mn_value_std) & 
           !is.na(Flow) &
           !is.na(post_lg) &
           CharacteristicName == "NO3_NO2" &
           yrs_lg >= -3 & 
           yrs_lg <= 3 &
           n() >= 8) %>%
  ungroup() %>%
  group_by(usgs_site) %>%
  filter(all(c("post", "pre") %in% post_lg)) %>%
  select(c(usgs_site:catchment_area, pctburn_lg, yrs_lg:post_lg))

write.csv(NCQ, here("USGS_data", "data_summaries", "NO3CQ.csv"), row.names = FALSE)

NCQ.pl <- chemQfire %>% group_by(usgs_site, CharacteristicName, post_lg) %>% 
  filter(!is.na(mn_value_std) & 
           !is.na(Flow) &
           !is.na(post_lg) &
           CharacteristicName == "NO3_NO2" &
           yrs_lg >= -3 & 
           yrs_lg <= 3 &
           n() >= 8) %>%
  ungroup() %>%
  group_by(usgs_site) %>%
  filter(all(c("post", "pre") %in% post_lg)) %>% 
  ggplot(aes(x = log(Flow), y = log(mn_value_std), group = post_lg)) +
  geom_point(aes(color = post_lg)) +
  ylab("log(NO3_NO2)") +
  facet_wrap(~usgs_site, scales = "free")

group_by(ID) %>%
  filter(all(c("up", "down") %in% DIR) )


#################
## Time series ##
#################
NO3NO2 <- large %>% group_by(usgs_site, CharacteristicName, post_lg) %>% 
                    filter(CharacteristicName == "NO3_NO2" &
                           !is.na(post_lg) &
                           n() > 10) %>%
                    ggplot(aes(x = Date_UTC, y = mn_value_std, group = post_lg)) +
                      geom_point(aes(color = post_lg)) +
                      ylab("NO3_NO2") +
                      facet_wrap(~usgs_site, scales = "free")

ggsave(NO3NO2, path = here("USGS_data", "plots"), file = "lg_ts_NO3NO2.pdf", width = 15, height = 14, units = "in")

turb <- large %>% group_by(usgs_site, CharacteristicName, post_lg) %>% 
                  filter(!is.na(mn_value_std) & 
                        #!is.na(Flow) &
                        !is.na(post_lg) &
                        n() > 10 &
                        CharacteristicName == "Turbidity") %>%
                  ggplot(aes(x = Date_UTC, y = mn_value_std, group = post_lg)) +
                    geom_point(aes(color = post_lg)) +
                    ylab("Turbidity") +
                    facet_wrap(~usgs_site, scales = "free")

ggsave(turb, path = here("USGS_data", "plots"), file = "lg_ts_turb.pdf", width = 16, height = 14, units = "in")

Kp <- large %>% group_by(usgs_site, CharacteristicName, post_lg) %>% 
                filter(!is.na(mn_value_std) & 
                      #!is.na(Flow) &
                       !is.na(post_lg) &
                       n() > 10 &
                       CharacteristicName == "Potassium") %>%
                ggplot(aes(x = Date_UTC, y = mn_value_std, group = post_lg)) +
                  geom_point(aes(color = post_lg)) +
                  ylab("Potassium") +
                  facet_wrap(~usgs_site, scales = "free")

ggsave(Kp, path = here("USGS_data", "plots"), file = "lg_ts_Potassium.pdf", width = 16, height = 14, units = "in")





