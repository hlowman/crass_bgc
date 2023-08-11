# MARSS models 
# Script started February 7, 2023
# Heili Lowman, Alex Webster

# The following script will be used to create the timeseries figures
# for conductivity NH4.

# Setup
library(tidyverse)
library(lubridate)
library(MARSS)
library(naniar)
library(emoGG)
library(calecopal)
library(patchwork)

# Load in data used for MARSS models.
df <- readRDS("data_working/marss_data_sb_080823.rds")
df_cond <- readRDS("data_working/marss_data_sb_vc_nolegacies_081123.rds")

#### NH4 ####

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
  scale_color_manual(values = cal_palette("sbchannel"), guide = "none") +
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

#### Sp. Conductivity ####

# Need to make new columns for designating pre- and post-fire dates
# for plotting purposes.
df_cond <- df_cond %>%
  mutate(day = 1) %>%
  mutate(date = make_date(year, month, day)) %>%
  # if multiple fires in a watershed occurred, I am coloring based
  # on the first of the fires to occur
  mutate(post_fire = case_when(site == "AB00" & date > as.Date("2009-05-01") |
                                 site == "GV01" & date > as.Date("2004-06-01") |
                                 site == "HO00" & date > as.Date("2004-06-01") |
                                 site == "RS02" & date > as.Date("2008-11-01") |
                                 site == "EFJ" & date > as.Date("2011-06-26") |
                                 site == "RED" & date > as.Date("2013-05-31") |
                                 site == "RSA" & date > as.Date("2011-06-26") |
                                 site == "RSAW" & date > as.Date("2011-06-26") ~ "A",
                               TRUE ~ "B")) %>%
  # and re-order sites for better plotting of sequential fires
  mutate(site_factor = factor(site, levels = c("EFJ", "RSA", "RSAW", "RED", 
                                               "GV01", "HO00", "RS02", "AB00")))

# Create figure.
(cond_fig_vc <- ggplot(df_cond %>%
                         filter(region == "VC"), 
                       aes(x = date, y = mean_cond_uScm, 
                           color = site_factor, alpha = post_fire)) +
    geom_point(size = 5) +
    geom_line() +
    scale_color_manual(values = c("#BED6B3", "#92A587", "#4A5438", "#2F3525"), 
                       guide = "none") +
    scale_alpha_manual(values = c(1, 0.3), guide = "none") +
    labs(x = "Date", 
         y = "Mean Monthly Specific Conductivity") +
    facet_grid(site_factor~.) +
    theme_bw() +
    theme(axis.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          strip.text.x = element_text(size = 24),
          strip.text.y = element_text(size = 24),
          strip.background = element_rect(colour="white", fill="white")))

(cond_fig_sb <- ggplot(df_cond %>%
                         filter(region == "SB"), 
                       aes(x = date, y = mean_cond_uScm, 
                           color = site_factor, alpha = post_fire)) +
    geom_point(size = 5) +
    geom_line() +
    scale_color_manual(values = cal_palette("sbchannel"), guide = "none") +
    scale_alpha_manual(values = c(1, 0.3), guide = "none") +
    labs(x = "Date") +
    facet_grid(site_factor~.) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          strip.text.x = element_text(size = 24),
          strip.text.y = element_text(size = 24),
          strip.background = element_rect(colour="white", fill="white")))

(cond_fig <- cond_fig_vc + cond_fig_sb)

# Export plot.
ggsave(("TS_SB_VC_SpCond_081123.png"),
       path = "figures",
       width = 60,
       height = 30,
       units = "cm"
)

# End of script.
