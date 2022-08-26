# SB MARSS models for other solutes
# Script started August 26, 2022
# Heili Lowman, Alex Webster

#### READ ME ####

# The following script will prepare data and run MARSS analyses at SBC sites for the CRASS project.
# MARSS modeling sections each contain 1 model or 1 model comparision and restart the script and import data anew each time.
# Be sure to load libraries before starting modeling sections. 
# Much of this code has been copied over from the "sbc_vcnp_firexppt&legacies.R" script written by Alex, but removing any of the NM site processing/modeling.

#### Load packages - ***do this first for all sections!*** ####
library(tidyverse)
library(lubridate)
library(MARSS)
library(naniar) 
library(beepr)

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))

#### Compile data ####

# Load datasets
# Stream Chemistry - all sites
# see script "sbc_chem_compilation.R" for the code used to tidy/generate this dataset
chem <- readRDS("data_working/SBchem_edited_120321.rds")

# Precipitation - all sites
precip <- readRDS("data_working/SBprecip_edited_120121.rds")

# Fire Events - all sites
fire <- readRDS("data_working/SBfire_edited_072922.rds")

# Site Location information
location <- read_csv("data_raw/sbc_sites_stream_hydro.csv")

# I first need to identify the matching stream sites for the precip data
precip_ed <- precip %>%
  mutate(sitecode_match = factor(case_when(
    sitecode == "GV202" ~ "GV01",
    sitecode == "BARA" ~ "HO00", # HO201 doesn't start until 2002
    sitecode == "RG202" ~ "RG01",
    sitecode == "SMPA" ~ "SP02",
    sitecode == "GORY" ~ "AT07",
    sitecode == "CAWTP" ~ "AB00",
    sitecode == "STFS" ~ "MC06", # BOGA doesn't start until 2005
    sitecode == "ELDE" ~ "RS02"))) %>%
  dplyr::rename(cumulative_precip_mm = c_precip_mm,
                sitecode_precip = sitecode) %>%
  mutate(Day = 1) %>% # new column of "days"
  mutate(Date = make_date(Year, Month, Day))

# check edited fire data
sitez = c("AB00", "GV01", "MC06", "RG01", "RS02", "HO00")

ggplot(fire %>% filter(site %in% sitez), 
       aes(x = date, y = fire_pa)) +
  geom_point() +
  labs(x = "Date") +
  facet_wrap(.~site, scales = "free",
             ncol = 1) +
  theme_bw()

ggplot(fire %>% filter(site %in% sitez), 
       aes(x = date, y = fire_perc_ws)) +
  geom_point() +
  labs(x = "Date") +
  facet_wrap(.~site, scales = "free",
             ncol = 1) +
  theme_bw()

## Timeframe Selection

# Examine precip data coverage for MARSS timescale delineation
precip_ed %>%
  drop_na(sitecode_match) %>%
  ggplot(aes(x = Year, y = sitecode_match, color = sitecode_match)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none")

# So, covariate data, which cannot be missing, can run from a maximum of 9/2002 to 7/2016.

# To avoid strange gaps in data, I'm going to start with the fire data,
# since I know it extends from 9/1/2002 to 7/1/2016. This will ensure all other datasets are
# joined to these dates in full (which was causing problems earlier).
fire_precip <- left_join(fire, precip_ed, by = c("site" = "sitecode_match", "date" = "Date"))

fire_precip <- fire_precip %>%
  mutate(year = year(date),
         month = month(date)) # and the Year/Month didn't populate, so adding in new columns

# Then, left join with chemistry so as not to lose any data.
dat <- left_join(fire_precip, chem, by = c("site", "year" = "Year", "month" = "Month"))

# Adding in dummy covariates by season
n_months <- dat$date %>%
  unique() %>%
  length()

seas_1 <- sin(2 * pi * seq(n_months) / 12)
seas_2 <- cos(2 * pi * seq(n_months) / 12)

dat <- dat %>%
  mutate(Season1 = rep(seas_1, 8),
         Season2 = rep(seas_2, 8),
         index = rep(seq(1,166), 8))

# AJW: replace NaNs with NAs
dat[is.nan(dat)] = NA

# And, inspect dataset for missing covariate data.
sum(is.na(dat$fire)) # 0
sum(is.na(dat$cumulative_precip_mm)) # 0
# Great!

# Trim down to columns of interest
dat_select <- dat %>%
  select(year, month, index, site, ws_area_m2,
         cumulative_precip_mm, 
         fire_ID, ig_date, fire_pa, ws_fire_area_m2, fire_perc_ws,
         mean_nh4_uM, mean_no3_uM, mean_po4_uM, 
         Season1, Season2) %>%
  mutate(region = "SB")

# And export to save progress
saveRDS(dat_select, "data_working/marss_data_sb_N_P_082622.rds")

#### Import compiled data and add fire x ppt interactions and legacy effects ####

dat1 <- readRDS("data_working/marss_data_sb_N_P_082622.rds")

# reorganize
names(dat1)
dat2 = dat1[,c(1,2,3,17,4,5, #"year" "month" "index" "region" "site" "ws_area_m2"
               6, # cumulative_precip_mm
               7:11, # fire info
               15,16, # seasonal effects
               12:14)] # solutes

qplot(index, fire_pa, data=dat2, colour=site, geom="path", facets = "region")

# create df of sites x fire info
firez = unique(dat1[,c(4,7,8,10,11)])
firez = firez[complete.cases(firez),] # removes NAs
firez$ig_date = gsub("/","-",firez$ig_date) # reformat date
firez$year = year(as.Date(firez$ig_date)) # create year column
firez$month = month(as.Date(firez$ig_date)) # create month column
firez$date = as.Date(paste(firez$year, firez$month, "01", sep="-")) # reformat date again

#### Add 6 m window to fire effect ####

# fire_pa_6m: 1=fire ignited in prior 6 m, 0=no fire ignitions in 6 m
# creating 6 months of dates
firedates_6m = c(firez$date,
                 firez$date + 31*1,
                 firez$date + 31*2,
                 firez$date + 31*3,
                 firez$date + 31*4,
                 firez$date + 31*5)

# creating dataset with fire data repeated for 6 months
firedates_6m = data.frame(year = year(firedates_6m),
                          month = month(firedates_6m),
                          site = rep(firez$site, 6),
                          fire_pa_6m = 1,
                          ws_fire_area_m2_6m = rep(firez$ws_fire_area_m2, 6),
                          fire_perc_ws_6m = rep(firez$fire_perc_ws, 6))

# join with original data
dat3 = left_join(dat2, firedates_6m, by=c("year","month","site"))

# remove NAs to be zeros in fire info columns so they can be used as covariates
dat3[,18:20][is.na(dat3[,18:20])] = 0

qplot(index, fire_pa_6m, data=dat3, colour=site, geom="path", facets = "region")
qplot(index, fire_perc_ws_6m, data=dat3, colour=site, geom="point", facets = "region")

# And export to save progress
saveRDS(dat3, "data_working/marss_data_sb_N_P_6mon_082622.rds")

#### Add fire 6 m window legacy effects ####

# Waiting for Alex's code re: persistent legacy effects

#### Add fire x ppt legacy effects ####

# for fire pa
dat3$fire_pa_6m_ppt = dat3$fire_pa_6m*dat3$cumulative_precip_mm

# for fire area
dat3$ws_fire_area_m2_6m_ppt = dat3$ws_fire_area_m2_6m*dat3$cumulative_precip_mm

# for fire perc burn ws
dat3$fire_perc_ws_6m_ppt = dat3$fire_perc_ws_6m*dat3$cumulative_precip_mm

#### Add fire 2 m window legacy effects ####

# Also skipping for now.

#### Examine data closely ####

# precip
ggplot(dat3, aes(x = index, y = cumulative_precip_mm, color = site)) +
  geom_point() +
  facet_wrap(.~site) +
  theme(legend.position = "none")

# NH4 (0-60 uM)
ggplot(dat3, aes(x = index, y = mean_nh4_uM, color = site)) +
  geom_point() +
  facet_wrap(.~site) +
  theme(legend.position = "none")

# NO3 (0-600 uM)
ggplot(dat3, aes(x = index, y = mean_no3_uM, color = site)) +
  geom_point() +
  facet_wrap(.~site) +
  theme(legend.position = "none")

# PO4 (0-15 uM)
ggplot(dat3, aes(x = index, y = mean_po4_uM, color = site)) +
  geom_point() +
  facet_wrap(.~site) +
  theme(legend.position = "none")

#### Remove outliers ####

# Keeping all data points for now.

#### Export data ready for modeling ####

saveRDS(dat3, "data_working/marss_data_sb_082622.rds")

#### Plot data ready for modeling ####

# remove data
rm(list=ls())

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

# load data with fire x ppt interactions and legacy effects
dat <- readRDS("data_working/marss_data_sb_082622.rds")

# add date
dat$date = as.Date(paste(dat$year, dat$month, "01", sep="-"))

# examine data to determine which sites to use for modeling

# NH4 in all sites
ggplot(dat, aes(x = date, y = mean_nh4_uM)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="AT07" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RG01" & fire_pa==1), aes(xintercept=date),
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="SP02" & fire_pa==1), aes(xintercept=date), 
             colour="red")

# NO3 in all sites
ggplot(dat, aes(x = date, y = mean_no3_uM)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="AT07" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RG01" & fire_pa==1), aes(xintercept=date),
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="SP02" & fire_pa==1), aes(xintercept=date), 
             colour="red")

# PO4 in all sites
ggplot(dat, aes(x = date, y = mean_po4_uM)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="AT07" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="RG01" & fire_pa==1), aes(xintercept=date),
             colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="SP02" & fire_pa==1), aes(xintercept=date), 
             colour="red")

# select sites for modeling
# these have the longest most complete ts for NH4/NO3/PO4 and have coverage before and after fires): 

# RG01 has a fire on the last data point
# SP02 has only 6 pre-fire data points

# AB00, AT07, GV01, HO00, MC06, RS02
sitez = c("AB00", "AT07", "GV01", "HO00", "MC06", "RS02")
dat = dat[dat$site %in% sitez,]

# Examine remaining covariates with the filtered dataset

# Fire area in all sites - check to be sure lines up with fire dates
ggplot(dat, aes(x = date, y = ws_fire_area_m2)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="AT07" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date),
             colour="red")+
  # geom_vline(data=filter(dat, site=="RG01" & fire_pa==1), aes(xintercept=date), 
  #            colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")
  # geom_vline(data=filter(dat, site=="SP02" & fire_pa==1), aes(xintercept=date), 
  #            colour="red")

# Fire % of watershed burned in all sites - check to be sure lines up with fire dates
ggplot(dat, aes(x = date, y = fire_perc_ws)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="AT07" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  # geom_vline(data=filter(dat, site=="RG01" & fire_pa==1), aes(xintercept=date), 
  #            colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")
  # geom_vline(data=filter(dat, site=="SP02" & fire_pa==1), aes(xintercept=date), 
  #            colour="red")+

# ppt in all sites - check for outliers
ggplot(dat, aes(x = date, y = cumulative_precip_mm)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="AT07" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  # geom_vline(data=filter(dat, site=="RG01" & fire_pa==1), aes(xintercept=date), 
  #            colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")
  # geom_vline(data=filter(dat, site=="SP02" & fire_pa==1), aes(xintercept=date), 
  #            colour="red")+

# Fire p/a X ppt in all sites - check if ppt interacts w/ fire in 6 mo timeframe
ggplot(dat, aes(x = date, y = fire_pa_6m_ppt)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="AT07" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  # geom_vline(data=filter(dat, site=="RG01" & fire_pa==1), aes(xintercept=date), 
  #            colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")
  # geom_vline(data=filter(dat, site=="SP02" & fire_pa==1), aes(xintercept=date), 
  #            colour="red")+

# Fire area X ppt in all sites - again check for interaction with precip
ggplot(dat, aes(x = date, y = ws_fire_area_m2_6m_ppt)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="AT07" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  # geom_vline(data=filter(dat, site=="RG01" & fire_pa==1), aes(xintercept=date), 
  #            colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")
  # geom_vline(data=filter(dat, site=="SP02" & fire_pa==1), aes(xintercept=date), 
  #            colour="red")+

# Fire % watershed burned X ppt in all sites - again check for interaction
ggplot(dat, aes(x = date, y = fire_perc_ws_6m_ppt)) +
  geom_point() +
  geom_line() +
  labs(x = "Date")+
  facet_wrap(region~site, scales = "free") +
  theme_bw()+
  geom_vline(data=filter(dat, site=="AB00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="AT07" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="GV01" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="HO00" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  geom_vline(data=filter(dat, site=="MC06" & fire_pa==1), aes(xintercept=date), 
             colour="red")+
  # geom_vline(data=filter(dat, site=="RG01" & fire_pa==1), aes(xintercept=date), 
  #            colour="red")+
  geom_vline(data=filter(dat, site=="RS02" & fire_pa==1), aes(xintercept=date), 
             colour="red")
  # geom_vline(data=filter(dat, site=="SP02" & fire_pa==1), aes(xintercept=date), 
  #            colour="red")+



# End of script.
