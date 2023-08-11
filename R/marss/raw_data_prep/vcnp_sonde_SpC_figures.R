# Valles Caldera National Preserve (vcnp)
# Sonde Data Review
# August 8, 2023
# Alex Webster, Heili Lowman

#### Read Me ####
# the purpose of this script is to put together some timeseries of sonde data from VCNP. 
# Data was provided by Mic from the Crossey lab off of their Aquarius data logs. 
# Data has been appended (stitched together) and needs to be QAQCed.

#### Setup ####

# load libraries
library(lubridate)
library(tidyverse)
library(gridExtra)
library(zoo)
library(xts)

#### Heili's code ####

# Load sonde data.
sonde_dat <- readRDS("data_raw/VCNP_sonde.rds")

# Load grab sample data.
grab_dat <- readRDS("data_working/VCNPchem_edited_120521.rds")

# Trim both datasets as appropriate.

# Make a dataframe.
sonde_dat_trim <- map_df(sonde_dat,
                         ~as.data.frame(.x), .id="site_name") %>%
  # Filter for sites of interest.
  filter(site_name %in% c("EFJ", "RED", "RSA", "RSAW")) %>%
  # Select only conductivity data.
  select(site_name, datetime_NM, SPC_uS_Value) %>%
  # And add column to denote data source
  mutate(source = "Sonde") %>%
  # And rename to match the grab sample dataset
  rename(DateTime = datetime_NM,
         cond_uScm = SPC_uS_Value)

# Also, need to remove QAQC to remove outliers from sonde data.
# Determine mean and sd conductivity values for each site.
summary <- sonde_dat_trim %>%
  group_by(site_name) %>%
  summarize(mean = mean(cond_uScm, na.rm = TRUE),
            sd = sd(cond_uScm, na.rm = TRUE)) %>%
  ungroup()

meanEFJ <- 96.51694
sdEFJ <- 37.80697
meanRED <- 61.84900
sdRED <- 37.91663
meanRSA <- 88.43305
sdRSA <- 20.54047
meanRSAW <- 117.55648
sdRSAW <- 24.41510

# And QAQC original dataset accordingly.
sonde_dat_trim_ed <- sonde_dat_trim %>%
  # First, removing all values equal to or less than 25 (not reasonable)
  filter(cond_uScm > 25) %>%
  # And, filter out outliers.
  mutate(cond_uScm_ed = case_when(site_name == "EFJ" & cond_uScm >= meanEFJ+(4*sdEFJ) ~ NA,
                                  site_name == "EFJ" & cond_uScm <= meanEFJ-(4*sdEFJ) ~ NA,
                                  site_name == "RED" & cond_uScm >= meanRED+(4*sdRED) ~ NA,
                                  site_name == "RED" & cond_uScm <= meanRED-(4*sdRED) ~ NA,
                                  site_name == "RSA" & cond_uScm >= meanRSA+(4*sdRSA) ~ NA,
                                  site_name == "RSA" & cond_uScm <= meanRSA-(4*sdRSA) ~ NA,
                                  site_name == "RSAW" & cond_uScm >= meanRSAW+(4*sdRSAW) ~ NA,
                                  site_name == "RSAW" & cond_uScm <= meanRSAW-(4*sdRSAW) ~ NA,
                                  TRUE ~ cond_uScm)) %>%
  select(site_name, DateTime, cond_uScm_ed, source) %>%
  rename(cond_uScm = cond_uScm_ed)

# Repeat with grab sample data.
grab_dat_trim <- grab_dat %>%
  mutate(site_name = case_when(site_code == "East Fork Jemez River" ~ "EFJ",
                               site_code == "Redondo Creek" ~ "RED",
                               site_code %in% c("San Antonio Creek - Toledo",
                                                "San Antonio Creek- Toledo",
                                                "San Antonio Creek -Toledo") ~ "RSA",
                               site_code == "San Antonio - West" ~ "RSAW")) %>%
  filter(site_name %in% c("EFJ", "RED", "RSA", "RSAW")) %>%
  select(site_name, DateTime, mean_cond_uScm) %>%
  rename(cond_uScm = mean_cond_uScm) %>%
  mutate(source = "Grab")

# Combine datasets.
all_dat <- full_join(grab_dat_trim, sonde_dat_trim_ed)

# Plot datasets to compare to one another.
# East Fork Jemez River
(EFJ_fig <- ggplot(all_dat %>%
                    filter(site_name == "EFJ"), 
                  aes(x = DateTime, y = cond_uScm, color = source)) +
  geom_point(alpha = 0.7) +
  geom_line(alpha = 0.6) +
  scale_color_manual(values = c("#0FB2D3", "#026779")) +
  geom_vline(xintercept = ymd_hms("2011-06-26 00:00:00"), color = "red") + # Las Conchas fire
  geom_vline(xintercept = ymd_hms("2016-05-11 00:00:00"), color = "red") + # prescribed fire
  labs(x = "Date-Time", y = "Spec. Conductivity (uS)", 
       title = "East Fork Jemez River (EFJ)") +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(rows = vars(source)))

# Export.
# ggsave(EFJ_fig,
#        filename = "figures/EFJ_grab_vs_sonde_081123.jpg",
#        width = 20,
#        height = 10,
#        units = "cm")

# Redondo Creek
(RED_fig <- ggplot(all_dat %>%
                     filter(site_name == "RED"), 
                   aes(x = DateTime, y = cond_uScm, color = source)) +
    geom_point(alpha = 0.7) +
    geom_line(alpha = 0.6) +
    scale_color_manual(values = c("#0FB2D3", "#026779")) +
    geom_vline(xintercept = ymd_hms("2013-05-31 00:00:00"), 
               color = "red") + # Thompson Ridge fire
    labs(x = "Date-Time", y = "Spec. Conductivity (uS)", 
         title = "Redondo Creek (RED)") +
    theme_bw() +
    theme(legend.position = "none") +
    facet_grid(rows = vars(source)))

# Export.
# ggsave(RED_fig,
#        filename = "figures/RED_grab_vs_sonde_081123.jpg",
#        width = 20,
#        height = 10,
#        units = "cm")

# San Antonio Creek - Toledo
(RSA_fig <- ggplot(all_dat %>%
                     filter(site_name == "RSA"), 
                   aes(x = DateTime, y = cond_uScm, color = source)) +
    geom_point(alpha = 0.7) +
    geom_line(alpha = 0.6) +
    scale_color_manual(values = c("#0FB2D3", "#026779")) +
    geom_vline(xintercept = ymd_hms("2011-06-26 00:00:00"), color = "red") + # Las Conchas fire
    labs(x = "Date-Time", y = "Spec. Conductivity (uS)", 
         title = "San Antonio Creek - Toledo (RSA)") +
    theme_bw() +
    theme(legend.position = "none") +
    facet_grid(rows = vars(source)))

# Export.
# ggsave(RSA_fig,
#        filename = "figures/RSA_grab_vs_sonde_081123.jpg",
#        width = 20,
#        height = 10,
#        units = "cm")

# San Antonio - West
(RSAW_fig <- ggplot(all_dat %>%
                     filter(site_name == "RSAW"), 
                   aes(x = DateTime, y = cond_uScm, color = source)) +
    geom_point(alpha = 0.7) +
    geom_line(alpha = 0.6) +
    scale_color_manual(values = c("#0FB2D3", "#026779")) +
    geom_vline(xintercept = ymd_hms("2011-06-26 00:00:00"), color = "red") + # Las Conchas fire
    geom_vline(xintercept = ymd_hms("2013-05-31 00:00:00"), 
               color = "red") + # Thompson Ridge fire
    labs(x = "Date-Time", y = "Spec. Conductivity (uS)", 
         title = "San Antonio - West (RSAW)") +
    theme_bw() +
    theme(legend.position = "none") +
    facet_grid(rows = vars(source)))

# Export.
# ggsave(RSAW_fig,
#        filename = "figures/RSAW_grab_vs_sonde_081123.jpg",
#        width = 20,
#        height = 10,
#        units = "cm")

#### Alex's code ####

# load data
sondez = readRDS("data_raw/VCNP_sonde.rds")

# extract and smooth SpC data

# going to start with just EFJ
EFJ_SpC = sondez[["EFJ"]] %>% 
  select(datetime_NM, SPC_uS_Value)

summary(EFJ_SpC)
# NA's   :0

# remove zeros, negative values, and one extreme value
EFJ_SpC$SPC_uS_Value[EFJ_SpC$SPC_uS_Value<1] = NA
EFJ_SpC$SPC_uS_Value[EFJ_SpC$SPC_uS_Value==0] = NA
EFJ_SpC$SPC_uS_Value[EFJ_SpC$SPC_uS_Value>1000] = NA

# Make univariate zoo time series#
EFJ_SpC_ts<-read.zoo(EFJ_SpC, index.column=1, format="%Y-%m-%d %H:%M:%S", tz="America/Denver")
# try smoothing to reduce noise
plot(EFJ_SpC_ts, type="l",main='Simple Moving Average (SMA)', ylim=c(-1,1000))
lines(rollmean(EFJ_SpC_ts,8), type="l",col='red')
# looks better

# turn smoothed ts back into df
EFJ_SpC_ts_sm = rollmean(EFJ_SpC_ts,8)
# revert back to df
EFJ_SpC_sm = as.data.frame(EFJ_SpC_ts_sm)
EFJ_SpC_sm$datetime_NM = as.POSIXct(row.names(EFJ_SpC_sm), tz="America/Denver")
names(EFJ_SpC_sm) = c("SPC_uS_Value","datetime_NM")

# plot by year
EFJ_SpC_sm$yr = year(EFJ_SpC_sm$datetime_NM)
EFJ_SpC_sm$doy = yday(EFJ_SpC_sm$datetime_NM)
ggplot(data= EFJ_SpC_sm, aes(y=SPC_uS_Value, x=doy))+
  geom_line()+
  facet_wrap(~yr)

###  JCP

JCP_SpC = sondez[["JCP"]] %>% 
  select(datetime_NM, SPC_uS_Value)

summary(JCP_SpC)
# NA's   :10134

# remove zeros, negative values, and one extreme value
JCP_SpC$SPC_uS_Value[JCP_SpC$SPC_uS_Value<1] = NA
JCP_SpC$SPC_uS_Value[JCP_SpC$SPC_uS_Value==0] = NA
JCP_SpC$SPC_uS_Value[JCP_SpC$SPC_uS_Value>1000] = NA

# Make univariate zoo time series#
JCP_SpC_ts<-read.zoo(JCP_SpC, index.column=1, format="%Y-%m-%d %H:%M:%S", tz="America/Denver")
# try smoothing to reduce noise
plot(JCP_SpC_ts, type="l",main='Simple Moving Average (SMA)', ylim=c(-1,1000))
lines(rollmean(JCP_SpC_ts,8), type="l",col='red')
# looks better

# turn smoothed ts back into df
JCP_SpC_ts_sm = rollmean(JCP_SpC_ts,8)
# revert back to df
JCP_SpC_sm = as.data.frame(JCP_SpC_ts_sm)
JCP_SpC_sm$datetime_NM = as.POSIXct(row.names(JCP_SpC_sm), tz="America/Denver")
names(JCP_SpC_sm) = c("SPC_uS_Value","datetime_NM")

# plot by year
JCP_SpC_sm$yr = year(JCP_SpC_sm$datetime_NM)
JCP_SpC_sm$doy = yday(JCP_SpC_sm$datetime_NM)
ggplot(data= JCP_SpC_sm, aes(y=SPC_uS_Value, x=doy))+
  geom_line()+
  facet_wrap(~yr)

###  RED

RED_SpC = sondez[["RED"]] %>% 
  select(datetime_NM, SPC_uS_Value)

summary(RED_SpC)
# NA's   :2

# remove zeros, negative values, and one extreme value
RED_SpC$SPC_uS_Value[RED_SpC$SPC_uS_Value<1] = NA
RED_SpC$SPC_uS_Value[RED_SpC$SPC_uS_Value==0] = NA
RED_SpC$SPC_uS_Value[RED_SpC$SPC_uS_Value>1000] = NA

# Make univariate zoo time series#
RED_SpC_ts<-read.zoo(RED_SpC, index.column=1, format="%Y-%m-%d %H:%M:%S", tz="America/Denver")
# try smoothing to reduce noise
plot(RED_SpC_ts, type="l",main='Simple Moving Average (SMA)', ylim=c(-1,1000))
lines(rollmean(RED_SpC_ts,8), type="l",col='red')
# looks better

# turn smoothed ts back into df
RED_SpC_ts_sm = rollmean(RED_SpC_ts,8)
# revert back to df
RED_SpC_sm = as.data.frame(RED_SpC_ts_sm)
RED_SpC_sm$datetime_NM = as.POSIXct(row.names(RED_SpC_sm), tz="America/Denver")
names(RED_SpC_sm) = c("SPC_uS_Value","datetime_NM")

# plot by year
RED_SpC_sm$yr = year(RED_SpC_sm$datetime_NM)
RED_SpC_sm$doy = yday(RED_SpC_sm$datetime_NM)
ggplot(data= RED_SpC_sm, aes(y=SPC_uS_Value, x=doy))+
  geom_line()+
  facet_wrap(~yr)

###  RSAW

RSAW_SpC = sondez[["RSAW"]] %>% 
  select(datetime_NM, SPC_uS_Value)

summary(RSAW_SpC)
# NA's   :0

# remove zeros, negative values, and one extreme value
RSAW_SpC$SPC_uS_Value[RSAW_SpC$SPC_uS_Value<1] = NA
RSAW_SpC$SPC_uS_Value[RSAW_SpC$SPC_uS_Value==0] = NA
RSAW_SpC$SPC_uS_Value[RSAW_SpC$SPC_uS_Value>1000] = NA

# Make univariate zoo time series#
RSAW_SpC_ts<-read.zoo(RSAW_SpC, index.column=1, format="%Y-%m-%d %H:%M:%S", tz="America/Denver")
# try smoothing to reduce noise
plot(RSAW_SpC_ts, type="l",main='Simple Moving Average (SMA)', ylim=c(-1,1000))
lines(rollmean(RSAW_SpC_ts,8), type="l",col='red')
# looks better

# turn smoothed ts back into df
RSAW_SpC_ts_sm = rollmean(RSAW_SpC_ts,8)
# revert back to df
RSAW_SpC_sm = as.data.frame(RSAW_SpC_ts_sm)
RSAW_SpC_sm$datetime_NM = as.POSIXct(row.names(RSAW_SpC_sm), tz="America/Denver")
names(RSAW_SpC_sm) = c("SPC_uS_Value","datetime_NM")

# plot by year
RSAW_SpC_sm$yr = year(RSAW_SpC_sm$datetime_NM)
RSAW_SpC_sm$doy = yday(RSAW_SpC_sm$datetime_NM)
ggplot(data= RSAW_SpC_sm, aes(y=SPC_uS_Value, x=doy))+
  geom_line()+
  facet_wrap(~yr)

###  RSA

RSA_SpC = sondez[["RSA"]] %>% 
  select(datetime_NM, SPC_uS_Value)

summary(RSA_SpC)
# NA's   :10819

# remove zeros, negative values, and one extreme value
RSA_SpC$SPC_uS_Value[RSA_SpC$SPC_uS_Value<1] = NA
RSA_SpC$SPC_uS_Value[RSA_SpC$SPC_uS_Value==0] = NA
RSA_SpC$SPC_uS_Value[RSA_SpC$SPC_uS_Value>1000] = NA

# Make univariate zoo time series#
RSA_SpC_ts<-read.zoo(RSA_SpC, index.column=1, format="%Y-%m-%d %H:%M:%S", tz="America/Denver")
# try smoothing to reduce noise
plot(RSA_SpC_ts, type="l",main='Simple Moving Average (SMA)', ylim=c(-1,1000))
lines(rollmean(RSA_SpC_ts,8), type="l",col='red')
# looks better

# turn smoothed ts back into df
RSA_SpC_ts_sm = rollmean(RSA_SpC_ts,8)
# revert back to df
RSA_SpC_sm = as.data.frame(RSA_SpC_ts_sm)
RSA_SpC_sm$datetime_NM = as.POSIXct(row.names(RSA_SpC_sm), tz="America/Denver")
names(RSA_SpC_sm) = c("SPC_uS_Value","datetime_NM")

# plot by year
RSA_SpC_sm$yr = year(RSA_SpC_sm$datetime_NM)
RSA_SpC_sm$doy = yday(RSA_SpC_sm$datetime_NM)
ggplot(data= RSA_SpC_sm, aes(y=SPC_uS_Value, x=doy))+
  geom_line()+
  facet_wrap(~yr)

# End of script.
