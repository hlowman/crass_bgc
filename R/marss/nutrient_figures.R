# SB MARSS models for other solutes
# Script started February 7, 2023
# Heili Lowman, Alex Webster

# The following script will be used to create the timeseries figure
# for NH4.

# Setup
library(tidyverse)
library(lubridate)
library(MARSS)
library(naniar)
library(emoGG)
library(calecopal)

# Load in data used for MARSS models.
df <- readRDS("data_working/marss_data_sb_080823.rds")

# Need to make new columns for designating pre- and post-fire dates
# for plotting purposes.
df <- df %>%
  mutate(day = 1) %>%
  mutate(date = make_date(year, month, day)) %>%
  # if multiple fires in a watershed occurred, I am coloring based
  # on the first of the fires to occur
  mutate(post_fire = case_when(site == "AB00" & date > as.Date("2009-05-01") |
                                 site == "GV01" & date > as.Date("2004-06-01") |
                                 site == "HO00" & date > as.Date("2004-06-01") |
                                 site == "RS02" & date > as.Date("2008-11-01") ~ "A",
                               TRUE ~ "B")) %>%
  # and re-order sites for better plotting of sequential fires
  mutate(site_factor = factor(site, levels = c("GV01", "HO00", "RS02", "AB00")))

# Create figure.
(nh4_fig <- ggplot(df, aes(x = date, y = vwm_nh4, color = site_factor, alpha = post_fire)) +
  geom_point(size = 5) +
  geom_line() +
  scale_y_log10() +
  scale_color_manual(values = cal_palette("oak"), guide = "none") +
  scale_alpha_manual(values = c(1, 0.3), guide = "none") +
  labs(x = "Date", 
       y = "Volume-Weighted Mean Monthly Ammonium Concentration",
       color = "Site") +
  facet_grid(site_factor~.) +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 20),
        strip.text.x = element_text(size = 24),
        strip.text.y = element_text(size = 24),
        strip.background = element_rect(colour="white", fill="white")))

# Export plot.
# ggsave(("TS_SB_NH4_081123.png"),
#        path = "figures",
#        width = 30,
#        height = 30,
#        units = "cm"
# )

# End of script.
