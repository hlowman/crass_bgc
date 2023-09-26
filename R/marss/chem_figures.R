# MARSS models 
# Script started February 7, 2023
# Heili Lowman, Alex Webster

# The following script will be used to create the timeseries figures
# for conductivity NH4.

#### Setup ####

# Load packages
library(tidyverse)
library(lubridate)
library(MARSS)
library(naniar)
library(emoGG)
library(calecopal)
library(patchwork)

# Load in data used for MARSS models.
df <- readRDS("data_working/marss_data_sb_092123.rds")
df_cond <- readRDS("data_working/marss_data_sb_vc_091123.rds")

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
(nh4_fig <- ggplot(df, aes(x = date, y = vwm_nh4, alpha = post_fire)) +
  geom_point(size = 5, color = "#6592D6") +
  geom_line(color = "#6592D6") +
  scale_y_log10() +
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
# ggsave(("TS_SB_NH4_092123.png"),
#        path = "figures",
#        width = 60,
#        height = 30,
#        units = "cm"
# )

#### Sp. Conductivity ####

# Need to make new columns for designating pre- and post-fire dates
# for plotting purposes.
df_cond <- df_cond %>%
  filter(site %in% c("AB00", "GV01", "HO00", "RS02",
                     "EFJ", "RED", "RSAW")) %>%
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
                                 site == "RSAW" & date > as.Date("2011-06-26") ~ "A",
                               TRUE ~ "B")) %>%
  # and re-order sites for better plotting of sequential fires
  mutate(site_factor = factor(site, levels = c("EFJ", "RSAW", "RED", 
                                               "GV01", "HO00", "RS02", "AB00")))

# Create figure.
(cond_fig_vc <- ggplot(df_cond %>%
                         filter(region == "VC"), 
                       aes(x = date, y = mean_cond_uScm, alpha = post_fire)) +
    geom_point(size = 5, color = "#92A587") +
    geom_line(color = "#92A587") +
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

(cond_fig_sb <- ggplot(df_cond %>%
                         filter(region == "SB"), 
                       aes(x = date, y = mean_cond_uScm, alpha = post_fire)) +
    geom_point(size = 5, color = "#6592D6") +
    geom_line(color = "#6592D6") +
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

design <- "AAABBB
           AAABBB
           AAABBB
           AAABBB
           AAA###"

(cond_fig <- cond_fig_sb + cond_fig_vc +
    plot_layout(design = design))

# Export plot.
# ggsave(("TS_SB_VC_SpCond_092623.png"),
#        path = "figures",
#        width = 60,
#        height = 30,
#        units = "cm")

#### Summary Stats ####

# Running a few additional summary statistics to include in the results section.

df_sb <- df_cond %>%
  filter(region == "SB")

df_sb_summary <- df_sb %>%
  summarize(mean_cond = mean(mean_cond_uScm, na.rm = TRUE),
            sd_cond = sd(mean_cond_uScm, na.rm = TRUE),
            min_cond = min(mean_cond_uScm, na.rm = TRUE),
            max_cond = max(mean_cond_uScm, na.rm = TRUE),
            mean_ppt = mean(cumulative_precip_mm, na.rm = TRUE),
            sd_ppt = sd(cumulative_precip_mm, na.rm = TRUE),
            min_ppt = min(cumulative_precip_mm, na.rm = TRUE),
            max_ppt = max(cumulative_precip_mm, na.rm = TRUE))

unique(df_sb$fire_perc_ws) #  0.0 45.5  8.8 17.1  9.3 68.4

df_vc <- df_cond %>%
  filter(region == "VC")

df_vc_summary <- df_vc %>%
  summarize(mean_cond = mean(mean_cond_uScm, na.rm = TRUE),
            sd_cond = sd(mean_cond_uScm, na.rm = TRUE),
            min_cond = min(mean_cond_uScm, na.rm = TRUE),
            max_cond = max(mean_cond_uScm, na.rm = TRUE),
            mean_ppt = mean(cumulative_precip_mm, na.rm = TRUE),
            sd_ppt = sd(cumulative_precip_mm, na.rm = TRUE),
            min_ppt = min(cumulative_precip_mm, na.rm = TRUE),
            max_ppt = max(cumulative_precip_mm, na.rm = TRUE))

unique(df_vc$fire_perc_ws) #  0.0 37.0 21.7 42.8  7.8 79.9 93.5

df_summary <- df %>%
  summarize(mean_NH4 = mean(vwm_nh4, na.rm = TRUE),
            sd_NH4 = sd(vwm_nh4, na.rm = TRUE),
            min_NH4 = min(vwm_nh4, na.rm = TRUE),
            max_NH4 = max(vwm_nh4, na.rm = TRUE),
            mean_NO3 = mean(vwm_no3, na.rm = TRUE),
            sd_NO3 = sd(vwm_no3, na.rm = TRUE),
            min_NO3 = min(vwm_no3, na.rm = TRUE),
            max_NO3 = max(vwm_no3, na.rm = TRUE),
            mean_PO4 = mean(vwm_po4, na.rm = TRUE),
            sd_PO4 = sd(vwm_po4, na.rm = TRUE),
            min_PO4 = min(vwm_po4, na.rm = TRUE),
            max_PO4 = max(vwm_po4, na.rm = TRUE))

# End of script.
