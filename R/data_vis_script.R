# SBC LTER data overview figures
# 2/11/2021
# Heili Lowman

# Load packages
library(tidyverse)
library(naniar)
library(viridis)

# Load data
dat <- read_csv("data/sbclter_stream_chemistry_allyears_registered_stations_20190628.csv") %>% # Load in LTER dataset
  replace_with_na_all(condition = ~.x == -999) # Replace all -999s with NAs.

dat_long <- dat %>%
  pivot_longer(
    cols = nh4_uM:spec_cond_uSpercm,
    names_to = "analyte",
    values_to = "value")

# Generate overall plot

# New facet label names for regional board variable
alabs <- c("NH4 (μM)", "NO3 (μM)", "PO4 (μM)", "TDN (μM)", "TDP (μM)", "Spec. Cond. (μS/cm)")
names(alabs) <- c("nh4_uM", "no3_uM", "po4_uM", "tdn_uM", "tdp_uM", "spec_cond_uSpercm")

fig_panel <- dat_long %>% # use pivoted dataset
  filter(analyte == "nh4_uM" |
           analyte == "no3_uM" |
           analyte == "po4_uM" |
           analyte == "tdn_uM" |
           analyte == "tdp_uM" |
           analyte == "spec_cond_uSpercm") %>% # filter for more available analytes
  mutate(analyte_f = factor(as.character(analyte))) %>% # recreate column
  ggplot(aes(x = timestamp_local, y = value)) +
  geom_point(aes(color = site_code)) +
  scale_color_viridis(discrete = TRUE) +
  labs(x = "Date",
       y = "Units",
       color = "Site") +
  facet_wrap(facets = vars(analyte_f), 
             scales = "free",
             labeller = labeller(analyte_f = alabs)) +
  theme_bw()

fig_panel

ggsave("figures/crass_testplot.png", fig_panel,
       height = 4, width = 10, units = "in")

# End of script.
