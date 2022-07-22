# MARSS model @ Santa Barbara Coastal LTER Sites
# October 24, 2021
# Heili Lowman, Alex Webster

# The following script will run a MARSS analysis at SBC LTER & NM Valles Caldera sites for the CRASS project.
# M - Multi-variate
# AR - Auto-regressive
# S - State
# S - Space

# A few notes about the process below:
# First, if you would like to skip over the data tidying portions
# you can load in the dataset on line 368 with all SB and VC data.

# Following a discussion with John Melack, the following sites have both enough chemistry
# and precip data, as well as fires that occurred within that timeframe to be included in the
# MARSS analysis performed below.

# Arroyo Burro - AB00
# Atascadero - AT07
# Gaviota - GV01
# Arroyo Hondo - HO00
# Mission Creek (at Rocky Nook) - MC06
# Refugio - RG01
# Rattlesnake - RS02
# San Pedro - SP02

# The following site may be added later, but a longer precipitation record needs to be
# identified for it - Bell Canyon (BC02).

# See "data_raw/VCNP_sonde_site_codes_names.csv" for other possible site names used across other files (E.g., sonde data, GIS, etc.)
# "Redondo Creek" ~ "RED",
# "East Fork Jemez River" ~ "EFJ",
# "San Antonio - West" ~ "RSAW",
# "San Antonio Creek - Toledo" ~ "RSA",
# "Indios Creek" ~ "IND",
# "Indios Creek - Post Fire (Below Burn)" ~ "IND_BB",
# "Indios Creek - above burn" ~ "IND_AB",
# "Sulfur Creek" ~ "SULF",

#### Setup ####

# Load packages
library(tidyverse)
library(lubridate)
library(here)
library(MARSS)
library(naniar)

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))

# Load datasets
# Stream Chemistry - all sites
chem <- readRDS("data_working/SBchem_edited_120321.rds")
chem_nm <- readRDS("data_working/VCNPchem_edited_120521.rds")
# Precipitation - all sites
precip <- readRDS("data_working/SBprecip_edited_120121.rds")
precip_nm <- readRDS("data_working/VCNPprecip_m_cum_edited_20220721.rds")
# Fire Events - all sites
fire <- readRDS("data_working/SBfire_edited_060622.rds")
fire_nm <- readRDS("data_working/VCNPfire_edited_060622.rds")
# Site Location information
location <- read_csv("data_raw/sbc_sites_stream_hydro.csv")
location_nm <- read_csv("data_raw/VCNP_sonde_site_codes_names.csv")
  
#### Filter & Joining Data ####

#### Santa Barbara ####

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

#### Valles Caldera ####

# Now to do the same for NM data
precip_nm_ed <- precip_nm %>%
  dplyr::rename(sitecode_precip = ID) %>%
  mutate(Day = 1) %>% # new column of "days"
  mutate(Date = as.Date(paste(year, month, "01", sep="-")))%>%
  mutate(Year = year,
         Month = month)

## Timeframe Selection
# Examine precip data coverage for MARSS timescale delineation
precip_nm_ed %>%
  ggplot(aes(x = Year, y = sitecode_precip, color = sitecode_precip)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none")

# So, covariate data, which cannot be missing, can run from a maximum of 2005 to 2021.
# To avoid strange gaps in data, I'm going to start with the fire data,
# since I know it extends from 6/1/2005 to 3/1/2019. This will ensure all other datasets are 
# joined to these dates in full (which was causing problems earlier).
fire_precip_nm <- left_join(fire_nm, precip_nm_ed, by = c("site" = "sitecode_precip", "date" = "Date"))

# Edit chemistry dataset to permit joining
namelist <- c("San Antonio Creek - Toledo", "San Antonio Creek- Toledo", "San Antonio Creek -Toledo")

# Make abbreviations appear at appropriate sites
chem_nm_ed <- chem_nm %>%
  mutate(site = factor(case_when(
    site_code == "Redondo Creek" ~ "RED",
    site_code == "East Fork Jemez River" ~ "EFJ",
    site_code == "San Antonio - West" ~ "RSAW",
    site_code %in% namelist ~ "RSA",
    site_code == "Indios Creek" ~ "IND",
    site_code == "Indios Creek - Post Fire (Below Burn)" ~ "IND_BB",
    site_code == "Indios Creek - above burn" ~ "IND_AB",
    site_code == "Sulfur Creek" ~ "SULF",
    TRUE ~ NA_character_))) %>%
  mutate(NH4_mgL = gsub("<", "", nh4_mgL),
         NO3_mgL = gsub("<", "", nO2_nO3_mgL),
         PO4_mgL = gsub("<", "", po4_mgL)) # remove all "<" symbols

# need to also add in LODs and convert analytes appropriately
mylist <- c("Contaminated")

chem_nm_ed2 <- chem_nm_ed %>%
  replace_with_na_at(.vars = c("NH4_mgL", "NO3_mgL", "PO4_mgL"),
                     condition = ~.x %in% mylist) %>% # Remove "contaminated" and replace with NA
  mutate(NH4_mgL_lod = as.numeric(ifelse(NH4_mgL <= 0.1, 0.05, NH4_mgL)),
         NO3_mgL_lod = as.numeric(ifelse(NO3_mgL <= 0.1, 0.05, NO3_mgL)),
         PO4_mgL_lod = as.numeric(ifelse(PO4_mgL <= 0.1, 0.05, PO4_mgL))) # Report low values at 1/2 LOD

# convert to uM and summarize by month
chem_nm_ed3 <- chem_nm_ed2 %>%
  mutate(nh4_uM = (NH4_mgL_lod/80.043)*1000, # convert to uM
         no3_uM = (NO3_mgL_lod/62.0049)*1000,
         po4_uM = (PO4_mgL_lod/94.9714)*1000) 

chem_nm_monthly <- chem_nm_ed3 %>%
  dplyr::group_by(site, Year, Month) %>%
  dplyr::summarize(mean_nh4_uM = mean(nh4_uM, na.rm = TRUE),
            mean_no3_uM = mean(no3_uM, na.rm = TRUE),
            mean_po4_uM = mean(po4_uM, na.rm = TRUE),
            mean_cond_uScm = mean(mean_cond_uScm, na.rm = TRUE)) %>%
  dplyr::ungroup()

# Then, left join with chemistry so as not to lose any data.
dat_nm <- left_join(fire_precip_nm, chem_nm_monthly, by = c("site", "Year", "Month"))

# Adding a plot to examine analyte availability
chem_nm_monthly[is.nan(chem_nm_monthly)] = NA

# plot of chem data availability
chem_nm_monthly %>%
  ggplot(aes(x = Year, y = site, color = site)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none")

temp = chem_nm_monthly[!is.na(chem_nm_monthly$mean_cond_uScm),]
temp %>%
  ggplot(aes(x = Year, y = site, color = site)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none")

dat_nm_trim <- dat_nm %>%
  filter(date >= "2005-06-01") 
# should be 166 records pre site to match "dat" above
# 166 records x 7 sites = 1162 records
# dates between SB and NM datasets don't need to match
# but lengths of time series do
table(dat_nm_trim$site)

# Adding in dummy covariates by season
n_months_nm <- dat_nm_trim$date %>%
  unique() %>%
  length()

seas_1_nm <- sin(2 * pi * seq(n_months_nm) / 12)
seas_2_nm <- cos(2 * pi * seq(n_months_nm) / 12)

dat_nm_trim <- dat_nm_trim %>%
  mutate(Season1 = rep(seas_1_nm, 7),
         Season2 = rep(seas_2_nm, 7),
         index = rep(seq(1,166), 7))

# AJW: replace NaNs with NAs
dat_nm_trim[is.nan(dat_nm_trim)] = NA

# Note for future me - be VERY careful with the joining above. Something weird was happening previously where precip data that IS present was simply dropping off.

# Also, because the Valles Caldera data isn't always monthly, some fire 1s were dropping off when they shouldn't, so Heili made edits 2/24/22 to rectify the situation.

## Added extra figures John requested to this code chunk. ##
# list of sites included in cond MARSS model below
sites_used <- c("EFJ", "RED", "RSA", "RSAW", "SULF")
# new facet labels for sites
site.labs <- c("East Fork Jemez River", "Redondo Creek", "San Antonio Creek - Toledo", "San Antonio Creek - West", "Sulfur Creek")
names(site.labs) <- sites_used

(no3_fig <- ggplot(chem_nm_ed3 %>% 
                    filter(site %in% sites_used), 
                  aes(x = DateTime, y = no3_uM)) +
  geom_point() +
  labs(x = "Date") +
  ylab(expression(paste(NO[3]^{"-"}," (μM)"))) +
  facet_wrap(.~site, scales = "free",
             ncol = 1,
             labeller = labeller(site = site.labs)) +
  theme_bw())

# ggsave(("figures/nm_sites_no3.png"),
#        width = 12,
#        height = 24,
#        units = "cm"
# )

(nh4_fig <- ggplot(chem_nm_ed3 %>% 
                     filter(site %in% sites_used), 
                   aes(x = DateTime, y = nh4_uM)) +
    geom_point() +
    labs(x = "Date") +
    ylab(expression(paste(NH[4]^{"+"}," (μM)"))) +
    facet_wrap(.~site, scales = "free",
               ncol = 1,
               labeller = labeller(site = site.labs)) +
    theme_bw())

# ggsave(("figures/nm_sites_nh4.png"),
#        width = 12,
#        height = 24,
#        units = "cm"
# )

(po4_fig <- ggplot(chem_nm_ed3 %>% 
                     filter(site %in% sites_used), 
                   aes(x = DateTime, y = po4_uM)) +
    geom_point() +
    labs(x = "Date") +
    ylab(expression(paste(PO[4]^{"3-"}," (μM)"))) +
    facet_wrap(.~site, scales = "free",
               ncol = 1,
               labeller = labeller(site = site.labs)) +
    theme_bw())

# ggsave(("figures/nm_sites_po4.png"),
#        width = 12,
#        height = 24,
#        units = "cm"
# )

(cond_fig <- ggplot(chem_nm_ed3 %>% 
                     filter(site %in% sites_used), 
                   aes(x = DateTime, y = mean_cond_uScm)) +
    geom_point() +
    labs(x = "Date") +
    ylab(expression(paste(Conductivity," (μS/cm)"))) +
    facet_wrap(.~site, scales = "free",
               ncol = 1,
               labeller = labeller(site = site.labs)) +
    theme_bw())

# ggsave(("figures/nm_sites_cond.png"),
#        width = 12,
#        height = 24,
#        units = "cm"
# )

## End of additional code making NM/VC figures for nutrients.

#### Join Full Dataset ####

# And finally, join the Santa Barbara and New Mexico datasets

# Created NM dataset for joining with SB data
dat_nm_select <- dat_nm_trim %>%
  # dplyr::rename(year = Year,
  #        month = Month) %>%
  select(year, month, site, 
         cumulative_precip_mm, 
         RED_Thompson, EFJ_Thompson, EFJ_Conchas, RSAW_Thompson, RSAW_Conchas, RSA_Conchas, IND_Conchas, IND_BB_Conchas, SULF_Thompson,
         RED_Thompson_d, EFJ_Thompson_d, EFJ_Conchas_d, RSAW_Thompson_d, RSAW_Conchas_d, RSA_Conchas_d, IND_Conchas_d, IND_BB_Conchas_d, SULF_Thompson_d,
         mean_nh4_uM, mean_no3_uM, mean_po4_uM, mean_cond_uScm, 
         Season1, Season2, index) %>%
  mutate(region = "VC") %>%
  mutate(AB00_Tea = 0,
         AB00_Jesusita = 0,
         AT07_Jesusita = 0,
         GV01_Gaviota = 0,
         HO00_Gaviota = 0,
         HO00_Sherpa = 0,
         MC06_Tea = 0,
         MC06_Jesusita = 0,
         RG01_Gaviota = 0,
         RG01_Sherpa = 0,
         RS02_Tea = 0,
         RS02_Jesusita = 0,
         SP02_Gap = 0,
         AB00_Tea_d = 0,
         AB00_Jesusita_d = 0,
         AT07_Jesusita_d = 0,
         GV01_Gaviota_d = 0,
         HO00_Gaviota_d = 0,
         HO00_Sherpa_d = 0,
         MC06_Tea_d = 0,
         MC06_Jesusita_d = 0,
         RG01_Gaviota_d = 0,
         RG01_Sherpa_d = 0,
         RS02_Tea_d = 0,
         RS02_Jesusita_d = 0,
         SP02_Gap_d = 0) # need to add in empty values for SB fire columns so the rbind below works

dat_select <- dat %>%
  select(year, month, site, 
         cumulative_precip_mm, 
         AB00_Tea,AB00_Jesusita,AT07_Jesusita,GV01_Gaviota,HO00_Gaviota,HO00_Sherpa,MC06_Tea,MC06_Jesusita,RG01_Gaviota,RG01_Sherpa,RS02_Tea,RS02_Jesusita,SP02_Gap,
         AB00_Tea_d,AB00_Jesusita_d,AT07_Jesusita_d,GV01_Gaviota_d,HO00_Gaviota_d,HO00_Sherpa_d,MC06_Tea_d,MC06_Jesusita_d,RG01_Gaviota_d,RG01_Sherpa_d,RS02_Tea_d,RS02_Jesusita_d,SP02_Gap_d,
         mean_nh4_uM, mean_no3_uM, mean_po4_uM, mean_cond_uScm, 
         Season1, Season2, index) %>%
  mutate(region = "SB") %>%
  mutate(RED_Thompson = 0, 
         EFJ_Thompson = 0, 
         EFJ_Conchas = 0, 
         RSAW_Thompson = 0, 
         RSAW_Conchas = 0, 
         RSA_Conchas = 0, 
         IND_Conchas = 0,
         IND_BB_Conchas = 0, 
         SULF_Thompson = 0,
         RED_Thompson_d = 0, 
         EFJ_Thompson_d = 0, 
         EFJ_Conchas_d = 0, 
         RSAW_Thompson_d = 0, 
         RSAW_Conchas_d = 0, 
         RSA_Conchas_d = 0, 
         IND_Conchas_d = 0,
         IND_BB_Conchas_d = 0, 
         SULF_Thompson_d = 0) # need to add in empty values for NM fire columns so the rbind below works

#dat_agu <- rbind(dat_select, dat_nm_select)
dat_new22 <- rbind(dat_select, dat_nm_select)


# replace fire NAs with zeros
dat_new22[,c(5:30,39:56)][is.na(dat_new22[,c(5:30,39:56)])] = 0

# And export to save progress
#saveRDS(dat_new22, "data_working/marss_data_sb_vc_022422.rds")
#saveRDS(dat_new22, "data_working/marss_data_sb_vc_060622.rds")
# with fixed NM ppt data (AJW):
saveRDS(dat_new22, "data_working/marss_data_sb_vc_072222.rds")

#### Model fitting ####

# Data: Stream Chemistry analytes (NH4, NO3, TDN, TPN, PO4, TDP, TP, TPP, TPC, TDS/TSS, SCond)
# Covariates: Month, Precip, Fire

# Note, we cannot analyze TSS, because we only have TSS data
# from the CA sites and TDS data from the NM sites.

#### NH4 ####
# HL - edited the following script on 3/23/22 to
# try and match Alex's successful modeling of conductivity.

# first, subset and pivot the NH4 data
dat_nh4 <- dat_new22 %>%
  select(site, index, Season1, Season2, 
         mean_nh4_uM, cumulative_precip_mm,
         AB00_Tea, AB00_Jesusita, AT07_Jesusita, GV01_Gaviota, HO00_Gaviota, HO00_Sherpa, MC06_Tea, MC06_Jesusita, RG01_Gaviota, RG01_Sherpa, RS02_Tea, RS02_Jesusita, RED_Thompson, EFJ_Thompson, EFJ_Conchas, RSAW_Thompson, RSAW_Conchas, RSA_Conchas) %>% 
  pivot_wider(names_from = site, values_from = c(mean_nh4_uM, cumulative_precip_mm, AB00_Tea, AB00_Jesusita, AT07_Jesusita, GV01_Gaviota, HO00_Gaviota, HO00_Sherpa, MC06_Tea, MC06_Jesusita, RG01_Gaviota, RG01_Sherpa, RS02_Tea, RS02_Jesusita, RED_Thompson, EFJ_Thompson, EFJ_Conchas, RSAW_Thompson, RSAW_Conchas, RSA_Conchas)) %>%
  select(index, Season1, Season2, 
         mean_nh4_uM_AB00, mean_nh4_uM_AT07, mean_nh4_uM_GV01, mean_nh4_uM_HO00, mean_nh4_uM_MC06, mean_nh4_uM_RG01, mean_nh4_uM_RS02, mean_nh4_uM_EFJ, mean_nh4_uM_RED, mean_nh4_uM_RSA, mean_nh4_uM_RSAW,
         cumulative_precip_mm_AB00, cumulative_precip_mm_AT07, cumulative_precip_mm_GV01, cumulative_precip_mm_HO00, cumulative_precip_mm_MC06, cumulative_precip_mm_RG01, cumulative_precip_mm_RS02, cumulative_precip_mm_EFJ, cumulative_precip_mm_RED, cumulative_precip_mm_RSA, cumulative_precip_mm_RSAW, 
         AB00_Tea_AB00, AB00_Jesusita_AB00, AT07_Jesusita_AT07, GV01_Gaviota_GV01, HO00_Gaviota_HO00, HO00_Sherpa_HO00, MC06_Tea_MC06, MC06_Jesusita_MC06, RG01_Gaviota_RG01, RG01_Sherpa_RG01, RS02_Tea_RS02, RS02_Jesusita_RS02, EFJ_Thompson_EFJ, EFJ_Conchas_EFJ, RED_Thompson_RED, RSA_Conchas_RSA, RSAW_Thompson_RSAW, RSAW_Conchas_RSAW)

# Remove the NaNs and replace with NAs.
dat_nh4[is.nan(dat_nh4)] <- NA

# log and scale transform response var
names(dat_nh4)
dat_nh4_log <- dat_nh4
dat_nh4_log[,4:14] <- log10(dat_nh4_log[,4:14])
dat_nh4_log[,4:14] <- scale(dat_nh4_log[,4:14])
sum(is.nan(dat_nh4_log[,4:14])) # 0 yay!!
sum(is.na(dat_nh4_log[,4:14])) # 724
range(dat_nh4_log[,4:14], na.rm = T)

### Plot response vars ###
# note - this is a HUGE plot
par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
plot(dat_nh4_log$mean_nh4_uM_AB00, type="o")
plot(dat_nh4_log$mean_nh4_uM_AT07, type="o")
plot(dat_nh4_log$mean_nh4_uM_GV01, type="o")
plot(dat_nh4_log$mean_nh4_uM_HO00, type="o")
plot(dat_nh4_log$mean_nh4_uM_MC06, type="o")
plot(dat_nh4_log$mean_nh4_uM_RG01, type="o")
plot(dat_nh4_log$mean_nh4_uM_RS02, type="o")
#plot(dat_nh4_log$mean_nh4_uM_SP02, type="o")
plot(dat_nh4_log$mean_nh4_uM_EFJ, type="o")
#plot(dat_nh4_log$mean_nh4_uM_IND, type="o")
#plot(dat_nh4_log$mean_nh4_uM_IND_AB, type="o")
#plot(dat_nh4_log$mean_nh4_uM_IND_BB, type="o")
plot(dat_nh4_log$mean_nh4_uM_RED, type="o")
plot(dat_nh4_log$mean_nh4_uM_RSA, type="o")
plot(dat_nh4_log$mean_nh4_uM_RSAW, type="o")
#plot(dat_nh4_log$mean_nh4_uM_SULF, type="o")
# SB sites clearly have better NH4 coverage than VC sites

#### Scenario 1 : all catchments are separate states #### 

# Pull out only NH4 data
names(dat_nh4_log)
# Removing shorter timeseries since this worked for conductivity
# Remaining sites: AB00, AT07, GV01, HO00, MC06, RG01, RS02
# EFJ, RED, RSA, RSAW, & SULF
# 12 sites total with 166 observations each
dat_dep <- t(dat_nh4_log[,c(4:10,12,16:19)])
row.names(dat_dep)

# Make covariate inputs
# Removing shorter timeseries since this worked for conductivity
# 33 covariates total with 166 observations each
dat_cov <- dat_nh4_log[,c(2:3, 
                          20:26, 28,32,33,34,35,
                          36:47, 50,51,49,54,52,53,56)]
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)

# make C matrix
# this matrix controls what covars predict what response vars; 
# in contrast to unconstrained: separately estimates ALL correlations of 
# predictor vars with x, i.e. how predictor vars drive ts dynamics. 
# This structure is ok given that covars are unique to each site, but it 
# would not be ok if there is a mix of shared and unique covars.
# 12 response variables - so 12 rows in the matrix
# 33 covariates in total - so 33 columns in the matrix

# without short ts sites
CC <- matrix(list(# season 1 covariate
                  "Season1", "Season1", "Season1", "Season1", 
                  "Season1", "Season1", "Season1", "Season1", 
                  "Season1", "Season1", "Season1", "Season1",
                  # season 2 covariate
                  "Season2", "Season2", "Season2", "Season2", 
                  "Season2", "Season2", "Season2", "Season2", 
                  "Season2", "Season2", "Season2", "Season2",
                  # precipitation covariates
                  "AB00_precip", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                  0, "AT07_precip", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                  0, 0, "GV01_precip", 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, "HO00_precip", 0, 0, 0, 0, 0, 0, 0, 0, 
                  0, 0, 0, 0, "MC06_precip", 0, 0, 0, 0, 0, 0, 0, 
                  0, 0, 0, 0, 0, "RG01_precip", 0, 0, 0, 0, 0, 0, 
                  0, 0, 0, 0, 0, 0, "RS02_precip", 0, 0, 0, 0, 0, 
                  0, 0, 0, 0, 0, 0, 0, "EFJ_precip", 0, 0, 0, 0, 
                  0, 0, 0, 0, 0, 0, 0, 0, "RED_precip", 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, "RSA_precip", 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "RSAW_precip", 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "SULF_precip",
                  # fire covariates
                  "AB00_Tea_fire", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  "AB00_Jesusita_fire", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, "AT07_Jesusita_fire", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, "GV01_Gaviota_fire", 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, "HO00_Gaviota_fire", 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, "HO00_Sherpa_fire", 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, "MC06_Tea_fire", 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, "MC06_Jesusita_fire", 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, "RG01_Gaviota_fire", 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, "RG01_Sherpa_fire", 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, "RS01_Tea_fire", 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, "RS01_Jesusita_fire", 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, "EFJ_Thompson_fire", 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, "EFJ_Conchas_fire", 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, "RED_Thompson_fire", 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, "RSA_Conchas_fire", 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "RSAW_Thompson_fire", 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "RSAW_Conchas_fire", 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "SULF_Thompson_fire"),
             12, 33)
                  
# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit = 100, allow.degen = TRUE, trace = 1),
                fit = TRUE) # used to get some priors

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", 
             inits=kemfit$par) # actual BFGS method model

fit_em <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 2000, allow.degen = TRUE, trace = 1),
             fit = TRUE) # actual EM method model

# see pg 5 in MARSS manual for notes on method BFGS vs method EM: EM algorithm gives more robust estimation for datasets replete with missing values and for high-dimensional models with various constraints. BFGS is faster and is good enough for some datasets. Typically, both should be tried.

# export model fit
saveRDS(fit, file = "data_working/marss_test_run/fit_020222_12state_nh4_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_RED_RSA_RSAW_SULF_mBFGS.rds")

### DIAGNOSES ###

## check for hidden errors
fit[["errors"]]
# none- Yay!!

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)

# What does this error message mean?
# Warning messages:
# 1: In MARSShessian(MLEobj, method = hessian.fun) :
#   MARSShessian: Hessian could not be inverted to compute the parameter var-cov matrix. parSigma set to NULL.  See MARSSinfo("HessianNA").
# 
# 2: In MARSSparamCIs(fit) :
#   MARSSparamCIs: No parSigma element returned by Hessian function.  See marssMLE object errors (MLEobj$errors)

saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_020222_12state_nh4_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_RED_RSA_RSAW_SULF_mBFGS.rds")

CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)

### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00", "AT07", "GV01", "HO00", 
             "MC06", "RG01", "RS02", "EFJ",
             "RED", "RSA", "RSAW", "SULF")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>% # rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "EFJ_precip", "RED_precip", 
                                             "RSA_precip", "RSAW_precip","SULF_precip",
                                             "AB00_Tea_fire","AB00_Jesusita_fire",
                                             "AT07_Jesusita_fire",
                                             "GV01_Gaviota_fire",
                                             "HO00_Gaviota_fire","HO00_Sherpa_fire",
                                             "MC06_Tea_fire", "MC06_Jesusita_fire",
                                             "RG01_Gaviota_fire", "RG01_Sherpa_fire",
                                             "RS02_Tea_fire", "RS02_Jesusita_fire",
                                             "EFJ_Thompson_fire", "EFJ_Conchas_fire",
                                             "RED_Thompson_fire",
                                             "RSA_Conchas_fire",
                                             "RSAW_Thompson_fire", "RSAW_Conchas_fire",
                                             "SULF_Thompson_fire")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "NH4 MARSS modeling results - 02/02/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))

# adding some formatting to match cond figure output
CIs_fit_ed2 = CIs_fit_ed[!(CIs_fit_ed$Site=="RSA" & CIs_fit_ed$Parameter=="RSAW_precip"),] 
CIs_fit_ed2 = CIs_fit_ed2[!(CIs_fit_ed2$Site=="RSA" & CIs_fit_ed2$Parameter=="RSAW_Thompson"),] 
CIs_fit_ed2 = CIs_fit_ed2[!(CIs_fit_ed2$Site=="RSA" & CIs_fit_ed2$Parameter=="RSAW_Conchas"),] 
CIs_fit_ed2$region = c(rep("Coastal California",33),rep("Subalpine New Mexico",22))

(RESULTS_ALL2 <- ggplot(CIs_fit_ed2, aes(Parameter, Est., color=region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=.7) +
    geom_point(position=position_dodge(width=0.3), size=5) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Ammonium MARSS modeling results - 2/2/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(vars(region, Site), scales = "free"))

# ggsave("figures/MARSS_12states_nh4_precip_fire_020222.png",
#        width = 40,
#        height = 20,
#        units = "cm")

## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- residuals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, 
                                    allow.degen=TRUE, trace=1), 
                     fit=TRUE) # default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), 
                  method = "BFGS", 
                  inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#           dAIC df
# mod.null  0.0  69
# fit       0.4  36
# RESULT: covar model is better than null

### Plot response vars ###
par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
plot(dat_dep[1,], type="o")
plot(dat_dep[2,], type="o")
plot(dat_dep[3,], type="o")
plot(dat_dep[4,], type="o")
plot(dat_dep[5,], type="o")
plot(dat_dep[6,], type="o")
plot(dat_dep[7,], type="o")
plot(dat_dep[8,], type="o")
plot(dat_dep[8,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  # forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}
# RESULT: ???
# Getting the same error message for both forecasts:
# Error in ts(x) : 'ts' object must have one or more observations

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  # qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"), 
  #        pch=16, 
  #        xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  # qqline(resids$model.residuals[i,])
  qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16, 
         xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# Also getting the following error here:
# Error in qqnorm.default(resids$state.residuals[i, ], main = paste(i, "state residuals"),  : y is empty or has only NAs

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

#### Scenario 1 : (SB only) all catchments are separate states #### 

# Heili revised the below section following a conversation with Tamara in 02/2022

# Pull out only NH4 data
names(dat_nh4_log)
# Remaining sites: AB00, AT07, GV01, HO00, MC06, RG01, RS02
# 7 sites total with 166 observations each
dat_dep <- t(dat_nh4_log[,c(4:10)])
row.names(dat_dep)

# Make covariate inputs
# 21 covariates total with 166 observations each
dat_cov <- dat_nh4_log[,c(2:3, # seasons
                          15:21, # precip
                          26:37)] # fires
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)

# make C matrix
# this matrix controls what covars predict what response vars; 
# in contrast to unconstrained: separately estimates ALL correlations of 
# predictor vars with x, i.e. how predictor vars drive ts dynamics. 
# This structure is ok given that covars are unique to each site, but it 
# would not be ok if there is a mix of shared and unique covars.
# 7 response variables - so 7 rows in the matrix
# 21 covariates in total - so 21 columns in the matrix

# without short ts sites
CC <- matrix(list(# season 1 covariate
  "Season1", "Season1", "Season1", "Season1", 
  "Season1", "Season1", "Season1",
  # season 2 covariate
  "Season2", "Season2", "Season2", "Season2", 
  "Season2", "Season2", "Season2",
  # precipitation covariates
  "AB00_precip", 0, 0, 0, 0, 0, 0,
  0, "AT07_precip", 0, 0, 0, 0, 0, 
  0, 0, "GV01_precip", 0, 0, 0, 0,
  0, 0, 0, "HO00_precip", 0, 0, 0,
  0, 0, 0, 0, "MC06_precip", 0, 0,
  0, 0, 0, 0, 0, "RG01_precip", 0,
  0, 0, 0, 0, 0, 0, "RS02_precip",
  # fire covariates
  "AB00_Tea_fire", 0, 0, 0, 0, 0, 0,
  "AB00_Jesusita_fire", 0, 0, 0, 0, 0, 0,
  0, "AT07_Jesusita_fire", 0, 0, 0, 0, 0,
  0, 0, "GV01_Gaviota_fire", 0, 0, 0, 0,
  0, 0, 0, "HO00_Gaviota_fire", 0, 0, 0,
  0, 0, 0, "HO00_Sherpa_fire", 0, 0, 0,
  0, 0, 0, 0, "MC06_Tea_fire", 0, 0,
  0, 0, 0, 0, "MC06_Jesusita_fire", 0, 0,
  0, 0, 0, 0, 0, "RG01_Gaviota_fire", 0,
  0, 0, 0, 0, 0, "RG01_Sherpa_fire", 0,
  0, 0, 0, 0, 0, 0, "RS02_Tea_fire",
  0, 0, 0, 0, 0, 0, "RS02_Jesusita_fire"),
  7, 21)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit = 100, allow.degen = TRUE, trace = 1),
                fit = TRUE) # used to get some priors

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", 
             inits=kemfit$par) # actual BFGS method model

# fit_em <- MARSS(y = dat_dep, model = mod_list,
#                 control = list(maxit = 2000, allow.degen = TRUE, trace = 1),
#                 fit = TRUE) # actual EM method model

# see pg 5 in MARSS manual for notes on method BFGS vs method EM: EM algorithm gives more robust estimation for datasets replete with missing values and for high-dimensional models with various constraints. BFGS is faster and is good enough for some datasets. Typically, both should be tried.

# export model fit
#saveRDS(fit, file = "data_working/marss_test_run/fit_032322_7state_nh4_AB00_AT07_GV01_HO00_MC06_RG01_RS02_mBFGS.rds")

### DIAGNOSES ###

## check for hidden errors
fit[["errors"]]
# NULL- Yay!!

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)

# No error messages!
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_032322_7state_nh4_AB00_AT07_GV01_HO00_MC06_RG01_RS02_mBFGS.rds")

CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)

### Plot Results for All Sites ###

# First, create dataset of all outputs
my_list <- c("AB00", "AT07", "GV01", "HO00", 
             "MC06", "RG01", "RS02")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>% # rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip",
                                             "AB00_Tea_fire","AB00_Jesusita_fire",
                                             "AT07_Jesusita_fire",
                                             "GV01_Gaviota_fire",
                                             "HO00_Gaviota_fire","HO00_Sherpa_fire",
                                             "MC06_Tea_fire", "MC06_Jesusita_fire",
                                             "RG01_Gaviota_fire", "RG01_Sherpa_fire",
                                             "RS02_Tea_fire", "RS02_Jesusita_fire")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "NH4 MARSS modeling results SB only - 03/23/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))

# adding some formatting to match cond figure output
CIs_fit_ed2 = CIs_fit_ed
CIs_fit_ed2$region = c(rep("Coastal California",33))

(RESULTS_ALL2 <- ggplot(CIs_fit_ed2, aes(Parameter, Est., color=region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=.7) +
    geom_point(position=position_dodge(width=0.3), size=5) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Ammonium (NH4) MARSS modeling results SB only - 3/23/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm"),# need to play with margins to make it all fit
          legend.position = "none") +
    facet_wrap(vars(region, Site), scales = "free"))

# ggsave("figures/MARSS_7states_SBonly_nh4_precip_fire_032322.png",
#        width = 25,
#        height = 15,
#        units = "cm")

## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- MARSSresiduals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, 
                                    allow.degen=TRUE, trace=1), 
                     fit=TRUE) # default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), 
                  method = "BFGS", 
                  inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

# dAIC - delta AIC
# 0.0 = always the value for the lowest model AIC
#          dAIC df
# fit       0.0 42
# null.fit 66.5 21
# RESULT: covar model is better than null - yay!!

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  # forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}

# Model residuals look alright - no major patterns, although need to look into results for model 6/RG01 (crossing threshold at 4, 9, and 15 month lags).
# State residuals look fine - no major patterns again, although should look at results for state 4/HO00 (crossing threshold at 1, 2, 11, and 12 month lags).

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  # qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"),
  #        pch=16,
  #        xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  # qqline(resids$model.residuals[i,])
  qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16,
         xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# State residuals look ok....state 2/AT07 does have a step in it.
# Model residuals look TERRIBLE...all flat lines.

# reset plotting window
dev.off()

#### Scenario 1 : (SB + SULF) all catchments are separate states #### 

# Pull out only NH4 data
names(dat_nh4_log)
# Removing shorter timeseries since this worked for conductivity
# Remaining sites: AB00, AT07, GV01, HO00, MC06, RG01, RS02
# & SULF
# 8 sites total with 166 observations each
dat_dep <- t(dat_nh4_log[,c(4:10,19)])
row.names(dat_dep)

# Make covariate inputs
# Removing shorter timeseries since this worked for conductivity
# 23 covariates total with 166 observations each
dat_cov <- dat_nh4_log[,c(2:3, 
                          20:26, 35,
                          36:47, 56)]
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)

# make C matrix
# this matrix controls what covars predict what response vars; 
# in contrast to unconstrained: separately estimates ALL correlations of 
# predictor vars with x, i.e. how predictor vars drive ts dynamics. 
# This structure is ok given that covars are unique to each site, but it 
# would not be ok if there is a mix of shared and unique covars.
# 8 response variables - so 8 rows in the matrix
# 23 covariates in total - so 23 columns in the matrix

# without short ts sites
CC <- matrix(list(# season 1 covariate
  "Season1", "Season1", "Season1", "Season1", 
  "Season1", "Season1", "Season1", "Season1",
  # season 2 covariate
  "Season2", "Season2", "Season2", "Season2", 
  "Season2", "Season2", "Season2", "Season2",
  # precipitation covariates
  "AB00_precip", 0, 0, 0, 0, 0, 0, 0,
  0, "AT07_precip", 0, 0, 0, 0, 0, 0, 
  0, 0, "GV01_precip", 0, 0, 0, 0, 0,
  0, 0, 0, "HO00_precip", 0, 0, 0, 0, 
  0, 0, 0, 0, "MC06_precip", 0, 0, 0,
  0, 0, 0, 0, 0, "RG01_precip", 0, 0, 
  0, 0, 0, 0, 0, 0, "RS02_precip", 0,
  0, 0, 0, 0, 0, 0, 0, "SULF_precip",
  # fire covariates
  "AB00_Tea_fire", 0, 0, 0, 0, 0, 0, 0, 
  "AB00_Jesusita_fire", 0, 0, 0, 0, 0, 0, 0, 
  0, "AT07_Jesusita_fire", 0, 0, 0, 0, 0, 0, 
  0, 0, "GV01_Gaviota_fire", 0, 0, 0, 0, 0, 
  0, 0, 0, "HO00_Gaviota_fire", 0, 0, 0, 0, 
  0, 0, 0, "HO00_Sherpa_fire", 0, 0, 0, 0, 
  0, 0, 0, 0, "MC06_Tea_fire", 0, 0, 0, 
  0, 0, 0, 0, "MC06_Jesusita_fire", 0, 0, 0, 
  0, 0, 0, 0, 0, "RG01_Gaviota_fire", 0, 0, 
  0, 0, 0, 0, 0, "RG01_Sherpa_fire", 0, 0, 
  0, 0, 0, 0, 0, 0, "RS01_Tea_fire", 0, 
  0, 0, 0, 0, 0, 0, "RS01_Jesusita_fire", 0, 
  0, 0, 0, 0, 0, 0, 0, "SULF_Thompson_fire"),
  8, 23)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit = 100, allow.degen = TRUE, trace = 1),
                fit = TRUE) # used to get some priors

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", 
             inits=kemfit$par) # actual BFGS method model

fit_em <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit = 2000, allow.degen = TRUE, trace = 1),
                fit = TRUE) # actual EM method model

# see pg 5 in MARSS manual for notes on method BFGS vs method EM: EM algorithm gives more robust estimation for datasets replete with missing values and for high-dimensional models with various constraints. BFGS is faster and is good enough for some datasets. Typically, both should be tried.

# export model fit
saveRDS(fit, file = "data_working/marss_test_run/fit_020322_8state_nh4_AB00_AT07_GV01_HO00_MC06_RG01_RS02_SULF_mBFGS.rds")

### DIAGNOSES ###

## check for hidden errors
fit[["errors"]]
# NULL - Yay!!

### Plot coef and coef estimates ###
## estimates
# hessian method is much faster but not ideal for final results
est_fit <- MARSSparamCIs(fit)

# The error message is back...
# What does this error message mean?
# Warning messages:
# 1: In MARSShessian(MLEobj, method = hessian.fun) :
#   MARSShessian: Hessian could not be inverted to compute the parameter var-cov matrix. parSigma set to NULL.  See MARSSinfo("HessianNA").
# 
# 2: In MARSSparamCIs(fit) :
#   MARSSparamCIs: No parSigma element returned by Hessian function.  See marssMLE object errors (MLEobj$errors)

saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_020322_8state_nh4_AB00_AT07_GV01_HO00_MC06_RG01_RS02_SULF_mBFGS.rds")

CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)

### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00", "AT07", "GV01", "HO00", 
             "MC06", "RG01", "RS02", "SULF")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>% # rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "SULF_precip",
                                             "AB00_Tea_fire","AB00_Jesusita_fire",
                                             "AT07_Jesusita_fire",
                                             "GV01_Gaviota_fire",
                                             "HO00_Gaviota_fire","HO00_Sherpa_fire",
                                             "MC06_Tea_fire", "MC06_Jesusita_fire",
                                             "RG01_Gaviota_fire", "RG01_Sherpa_fire",
                                             "RS02_Tea_fire", "RS02_Jesusita_fire",
                                             "SULF_Thompson_fire")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "NH4 MARSS modeling results - 02/02/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))

# adding some formatting to match cond figure output
CIs_fit_ed2 = CIs_fit_ed
CIs_fit_ed2$region = c(rep("Coastal California",31),rep("Subalpine New Mexico",4))

(RESULTS_ALL2 <- ggplot(CIs_fit_ed2, aes(Parameter, Est., color=region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=.7) +
    geom_point(position=position_dodge(width=0.3), size=5) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Ammonium MARSS modeling results SB + SULF - 2/3/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(vars(region, Site), scales = "free"))

# ggsave("figures/MARSS_8states_sn_sulf_nh4_precip_fire_020322.png",
#        width = 40,
#        height = 20,
#        units = "cm")

## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- residuals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, 
                                    allow.degen=TRUE, trace=1), 
                     fit=TRUE) # default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), 
                  method = "BFGS", 
                  inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#           dAIC df
# mod.null  0.0  47
# fit       4.4  24
# RESULT: covar model is better than null

### Plot response vars ###
# par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
# plot(dat_dep[1,], type="o")
# plot(dat_dep[2,], type="o")
# plot(dat_dep[3,], type="o")
# plot(dat_dep[4,], type="o")
# plot(dat_dep[5,], type="o")
# plot(dat_dep[6,], type="o")
# plot(dat_dep[7,], type="o")
# plot(dat_dep[8,], type="o")
# plot(dat_dep[8,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  # forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}
# RESULT: ???
# Getting the same error message for both forecasts:
# Error in ts(x) : 'ts' object must have one or more observations

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  # qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"), 
  #        pch=16, 
  #        xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  # qqline(resids$model.residuals[i,])
  qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16, 
         xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# Also getting the following error here:
# Error in qqnorm.default(resids$state.residuals[i, ], main = paste(i, "state residuals"),  : y is empty or has only NAs

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

#### Scenario 2 : catchments in two ecoregions #### 
# not using fire for now to simplify
# Make covariate inputs
dat_nh4 <- dat_nh4 %>%
  mutate(c_precip_SB = (cumulative_precip_mm_AB00 + cumulative_precip_mm_AT07 + cumulative_precip_mm_GV01 + cumulative_precip_mm_HO00 + cumulative_precip_mm_MC06 + cumulative_precip_mm_RG01 + cumulative_precip_mm_RS02 + cumulative_precip_mm_SP02)/8,
         c_precip_NM = (cumulative_precip_mm_EFJ + cumulative_precip_mm_IND + cumulative_precip_mm_IND_AB + cumulative_precip_mm_IND_BB + cumulative_precip_mm_RED + cumulative_precip_mm_RSA + cumulative_precip_mm_RSAW + cumulative_precip_mm_SULF)/8)

dat_cov2 <- dat_nh4[,c(2:3,36:37)]
dat_cov2 <- t(scale(dat_cov2))

CC2 <- matrix(list("Season1", "Season1",
                   "Season2", "Season2",
                   "c_precip_SB", 0,
                   0, "c_precip_NM"),2,4)

# Model setup
mod_list2 <- list(
  # tinitx = "zero", # setting initial state value to time = 0
  B = "diagonal and unequal",
  U = "zero", # zero: does NOT allow a drift term in process model to be estimated # removing due to lack of anticipated monotonic trend
  C = CC2, # see new matrix above
  c = dat_cov2, # we should probably de-mean and scale covariates to units of sd 
  Q = "diagonal and unequal", # diagonal and unequal: allows for and estimates the covariance matrix of process errors
  Z = factor(c("SantaBarbara", "SantaBarbara",
               "SantaBarbara", "SantaBarbara",
               "SantaBarbara", "SantaBarbara",
               "SantaBarbara", "SantaBarbara",
               "VallesCaldera", "VallesCaldera",
               "VallesCaldera", "VallesCaldera",
               "VallesCaldera", "VallesCaldera",
               "VallesCaldera", "VallesCaldera")),
  A = "zero",
  R = "zero" # diagonal and equal: allows for and estimates the covariance matrix of observations errors (may want to provide a number for this from method precision etc if possible) - changed to "zero" on 11/22 to "turn off" observation error
)

# Fit model
fit2 <- MARSS(y = dat_dep, model = mod_list2,
              control = list(maxit = 5000), method = "BFGS")

# fit2 <- MARSS(y = dat_dep, model = mod_list2,
#              control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

# export model fit
saveRDS(fit2, file = "data_working/marss_test_run/fit_120321_2states.rds")

#### Scenario 3 : all catchments in a single state #### 
# not using fire for now to simplify
# Make covariate inputs
dat_nh4 <- dat_nh4 %>%
  mutate(c_precip_avg = (cumulative_precip_mm_AB00 + cumulative_precip_mm_AT07 + cumulative_precip_mm_GV01 + cumulative_precip_mm_HO00 + cumulative_precip_mm_MC06 + cumulative_precip_mm_RG01 + cumulative_precip_mm_RS02 + cumulative_precip_mm_SP02 + cumulative_precip_mm_EFJ + cumulative_precip_mm_IND + cumulative_precip_mm_IND_AB + cumulative_precip_mm_IND_BB + cumulative_precip_mm_RED + cumulative_precip_mm_RSA + cumulative_precip_mm_RSAW + cumulative_precip_mm_SULF)/16)
         
dat_cov3 <- dat_nh4[,c(2:3,36)]
dat_cov3 <- t(scale(dat_cov3))

CC3 <- matrix(list("Season1", "Season2", "c_precip_avg"),1,3)

# Model setup
mod_list3 <- list(
  # tinitx = "zero", # setting initial state value to time = 0
  B = "diagonal and unequal",
  U = "zero", # zero: does NOT allow a drift term in process model to be estimated # removing due to lack of anticipated monotonic trend
  C = CC3, # see new matrix above
  c = dat_cov3, # we should probably de-mean and scale covariates to units of sd 
  Q = "diagonal and unequal", # diagonal and unequal: allows for and estimates the covariance matrix of process errors
  Z = matrix(1, nrow = 16, ncol = 1),
  A = "zero",
  R = "zero" # diagonal and equal: allows for and estimates the covariance matrix of observations errors (may want to provide a number for this from method precision etc if possible) - changed to "zero" on 11/22 to "turn off" observation error
)

# Fit model
fit3 <- MARSS(y = dat_dep, model = mod_list3,
                control = list(maxit = 5000), method = "BFGS")

# fit3 <- MARSS(y = dat_dep, model = mod_list3,
#              control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

# export model fit
saveRDS(fit3, file = "data_working/marss_test_run/fit_120321_1state.rds")

#### NO3 ####
# Going to use NO3 to add one catchment at a time to see what's going on with model fit
# And to simplify, for the AGU talk, taking out the fire variables for now

dat_no3 <- dat_agu %>%
  select(site, index, Season1, Season2, 
         mean_no3_uM, cumulative_precip_mm) %>% 
  pivot_wider(names_from = site, values_from = c(mean_no3_uM, cumulative_precip_mm)) %>%
  select(index, Season1, Season2, 
         mean_no3_uM_AB00, mean_no3_uM_AT07, mean_no3_uM_GV01, mean_no3_uM_HO00, mean_no3_uM_MC06, mean_no3_uM_RG01, mean_no3_uM_RS02, mean_no3_uM_SP02, mean_no3_uM_EFJ, mean_no3_uM_IND, mean_no3_uM_IND_AB, mean_no3_uM_IND_BB, mean_no3_uM_RED, mean_no3_uM_RSA, mean_no3_uM_RSAW, mean_no3_uM_SULF,
         cumulative_precip_mm_AB00, cumulative_precip_mm_AT07, cumulative_precip_mm_GV01, cumulative_precip_mm_HO00, cumulative_precip_mm_MC06, cumulative_precip_mm_RG01, cumulative_precip_mm_RS02, cumulative_precip_mm_SP02, cumulative_precip_mm_EFJ, cumulative_precip_mm_IND, cumulative_precip_mm_IND_AB, cumulative_precip_mm_IND_BB, cumulative_precip_mm_RED, cumulative_precip_mm_RSA, cumulative_precip_mm_RSAW, cumulative_precip_mm_SULF)

dat_no3[is.nan(dat_no3)] <- NA

# log and scale transform response var
names(dat_no3)
dat_no3_log = dat_no3
dat_no3_log[,4:19] = log10(dat_no3_log[,4:19])
dat_no3_log[,4:19] = scale(dat_no3_log[,4:19])
sum(is.nan(dat_no3_log[,4:19]))
sum(is.na(dat_no3_log[,4:19]))
range(dat_no3_log[,4:19], na.rm = T)

### Plot response vars ###
par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
plot(dat_no3_log$mean_no3_uM_AB00, type="o")
plot(dat_no3_log$mean_no3_uM_AT07, type="o")
plot(dat_no3_log$mean_no3_uM_GV01, type="o")
plot(dat_no3_log$mean_no3_uM_HO00, type="o")
plot(dat_no3_log$mean_no3_uM_MC06, type="o")
plot(dat_no3_log$mean_no3_uM_RG01, type="o")
plot(dat_no3_log$mean_no3_uM_RS02, type="o")
plot(dat_no3_log$mean_no3_uM_SP02, type="o")
plot(dat_no3_log$mean_no3_uM_EFJ, type="o")
plot(dat_no3_log$mean_no3_uM_IND, type="o")
plot(dat_no3_log$mean_no3_uM_IND_AB, type="o")
plot(dat_no3_log$mean_no3_uM_IND_BB, type="o")
plot(dat_no3_log$mean_no3_uM_RED, type="o")
plot(dat_no3_log$mean_no3_uM_RSA, type="o")
plot(dat_no3_log$mean_no3_uM_RSAW, type="o")
plot(dat_no3_log$mean_no3_uM_SULF, type="o")

par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
plot(dat_no3$mean_no3_uM_AB00, type="o")
plot(dat_no3$mean_no3_uM_AT07, type="o")
plot(dat_no3$mean_no3_uM_GV01, type="o")
plot(dat_no3$mean_no3_uM_HO00, type="o")
plot(dat_no3$mean_no3_uM_MC06, type="o")
plot(dat_no3$mean_no3_uM_RG01, type="o")
plot(dat_no3$mean_no3_uM_RS02, type="o")
plot(dat_no3$mean_no3_uM_SP02, type="o")
plot(dat_no3$mean_no3_uM_EFJ, type="o")
plot(dat_no3$mean_no3_uM_IND, type="o")
plot(dat_no3$mean_no3_uM_IND_AB, type="o")
plot(dat_no3$mean_no3_uM_IND_BB, type="o")
plot(dat_no3$mean_no3_uM_RED, type="o")
plot(dat_no3$mean_no3_uM_RSA, type="o")
plot(dat_no3$mean_no3_uM_RSAW, type="o")
plot(dat_no3$mean_no3_uM_SULF, type="o")


#### Scenario 1 : all catchments are separate states - 1 catchment #### 

# Pull out only response var
#dat_dep <- t(dat_nh4_log[,4:19])
# w/o IND sites:
#dat_dep <- t(dat_no3_log[,c(4:12,16:19)])
# starting with HO00
dat_dep <- t(dat_no3_log[,c(7)])

# Make covariate inputs
#dat_cov <- dat_nh4_log[,c(2:3,20:35)]
# w/o IND sites:
# dat_cov <- dat_nh4_log[,c(2:3,20:28,32:35)]
# starting with HO00
dat_cov <- dat_no3_log[,c(2:3,23)]
dat_cov <- t(scale(dat_cov))


#### make C matrix
# starting with HO00
CC <- matrix(list("Season1", 
                  "Season2", 
                  "HO00_precip"),1,3)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observtion model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model
# fit <- MARSS(y = dat_dep, model = mod_list,
#              control = list(maxit = 5000), method = "BFGS")
# see pg 5 in MARSS manual for notes on method BFGS vs method EM: EM algorithm gives more robust estimation for datasets replete with missing values and for high-dimensional models with various constraints. BFGS is faster and is good enough for some datasets. Typically, both should be tried.

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

# export model fit
saveRDS(fit, file = "data_working/marss_test_run/fit_120521_1state_HO00_EM.rds")


### DIAGNOSES ###
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## check for hidden errors
fit[["errors"]]
# lots of errors when IND sites were included, none when these were removed

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
#est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)

saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_1state_HO00_EM_hessian.rds")

CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)


### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00", "AT07", "GV01", "HO00", "MC06", "RG01", "RS02", "SP02","EFJ","RED","RSA","RSAW","SULF")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>% # rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "SP02_precip", 
                                             "EFJ_precip", "RED_precip", "RSA_precip", "RSAW_precip",
                                             "SULF_precip")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "NO3 MARSS modeling results - 12/5/2021") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))

## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- residuals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)
mod.null <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE)
# mod.null <- MARSS(y = dat_dep, model = mod_list_null,
#                    control = list(maxit = 5000), method = "BFGS")

bbmle::AICtab(fit, mod.null)

#           dAIC df
# fit      0.0  6 
# mod.null 0.1  3 
# RESULT: models are equivelent 

### Do resids have temporal autocorrelation? ###
par(mfrow=c(1,2),oma = c(0, 0, 2, 0))
forecast::Acf(resids$model.residuals[1,], main="1 model residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$state.residuals[1,], main="1 state residuals", na.action=na.pass, lag.max = 24)
mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)

### Do resids have temporal trend? ###
par(mfrow=c(1,2),oma = c(0, 0, 2, 0))
plot(resids$model.residuals[1,], ylab="model residual", xlab="", main="1 model residuals")
abline(h=0)
plot(resids$state.residuals[1,], ylab="state residual", xlab="", main="1 state residuals")
abline(h=0)
# plot(resids$model.residuals[2,], ylab="model residual", xlab="", main="RC01 model residuals")
# abline(h=0)
# plot(resids$state.residuals[2,], ylab="state residual", xlab="", main="RC01 state residuals")
# abline(h=0)
mtext("Do resids have temporal trend?", outer = TRUE, cex = 1.5)

### Are resids normal? ###
par(mfrow=c(1,2),oma = c(0, 0, 2, 0))
qqnorm(resids$model.residuals[1,], main="1 model residuals", pch=16, 
       xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[1,])[1]))
qqline(resids$model.residuals[1,])
qqnorm(resids$state.residuals[1,], main="1 state residuals", pch=16, 
       xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[1,])[1]))
qqline(resids$state.residuals[1,])
# qqnorm(resids$model.residuals[2,], main="2 model residuals", pch=16, 
#        xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[2,])[1]))
# qqline(resids$model.residuals[2,])
# qqnorm(resids$state.residuals[2,], main="3 state residuals", pch=16, 
#        xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[2,])[1]))
# qqline(resids$state.residuals[2,])
mtext("Are resids normal?", outer = TRUE, cex = 1.5)

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))


#### Scenario 1 : all catchments are separate states - 1 catchment #### 

# Pull out only response var
names(dat_no3_log)
# AB00
dat_dep <- t(dat_no3_log[,c(4)])
row.names(dat_dep)

# Make covariate inputs
dat_cov <- dat_no3_log[,c(2:3,
                          20)]
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)


#### make C matrix
# starting with HO00
CC <- matrix(list("Season1", 
                  "Season2", 
                  "AB00_precip"),1,3)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observtion model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model
# fit <- MARSS(y = dat_dep, model = mod_list,
#              control = list(maxit = 5000), method = "BFGS")
# see pg 5 in MARSS manual for notes on method BFGS vs method EM: EM algorithm gives more robust estimation for datasets replete with missing values and for high-dimensional models with various constraints. BFGS is faster and is good enough for some datasets. Typically, both should be tried.

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

# export model fit
saveRDS(fit, file = "data_working/marss_test_run/fit_120521_1state_AB00_EM.rds")


### DIAGNOSES ###
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## check for hidden errors
fit[["errors"]]
# lots of errors when IND sites were included, none when these were removed

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
#est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)

saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_1state_AB00_EM_hessian.rds")

CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)


### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00", "AT07", "GV01", "HO00", "MC06", "RG01", "RS02", "SP02","EFJ","RED","RSA","RSAW","SULF")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>% # rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "SP02_precip", 
                                             "EFJ_precip", "RED_precip", "RSA_precip", "RSAW_precip",
                                             "SULF_precip")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "NO3 MARSS modeling results - 12/5/2021") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))


## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- residuals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)
mod.null <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE)
# mod.null <- MARSS(y = dat_dep, model = mod_list_null,
#                    control = list(maxit = 5000), method = "BFGS")

bbmle::AICtab(fit, mod.null)

#           dAIC df
# fit      0.0  6 
# mod.null 4.4  3 
# RESULT: covar model is better than null

### Do resids have temporal autocorrelation? ###
par(mfrow=c(1,2),oma = c(0, 0, 2, 0))
forecast::Acf(resids$model.residuals[1,], main="1 model residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$state.residuals[1,], main="1 state residuals", na.action=na.pass, lag.max = 24)
mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)

### Do resids have temporal trend? ###
par(mfrow=c(1,2),oma = c(0, 0, 2, 0))
plot(resids$model.residuals[1,], ylab="model residual", xlab="", main="1 model residuals")
abline(h=0)
plot(resids$state.residuals[1,], ylab="state residual", xlab="", main="1 state residuals")
abline(h=0)
# plot(resids$model.residuals[2,], ylab="model residual", xlab="", main="RC01 model residuals")
# abline(h=0)
# plot(resids$state.residuals[2,], ylab="state residual", xlab="", main="RC01 state residuals")
# abline(h=0)
mtext("Do resids have temporal trend?", outer = TRUE, cex = 1.5)

### Are resids normal? ###
par(mfrow=c(1,2),oma = c(0, 0, 2, 0))
qqnorm(resids$model.residuals[1,], main="1 model residuals", pch=16, 
       xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[1,])[1]))
qqline(resids$model.residuals[1,])
qqnorm(resids$state.residuals[1,], main="1 state residuals", pch=16, 
       xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[1,])[1]))
qqline(resids$state.residuals[1,])
# qqnorm(resids$model.residuals[2,], main="2 model residuals", pch=16, 
#        xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[2,])[1]))
# qqline(resids$model.residuals[2,])
# qqnorm(resids$state.residuals[2,], main="3 state residuals", pch=16, 
#        xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[2,])[1]))
# qqline(resids$state.residuals[2,])
mtext("Are resids normal?", outer = TRUE, cex = 1.5)

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))


#### Scenario 1 : all catchments are separate states - 2 catchments #### 

# Pull out only response var
names(dat_no3_log)
# AB00, AT07
dat_dep <- t(dat_no3_log[,c(4,5)])
row.names(dat_dep)

# Make covariate inputs
dat_cov <- dat_no3_log[,c(2:3,
                          20,21)]
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)


#### make C matrix
# starting with HO00
CC <- matrix(list("Season1", "Season1", 
                  "Season2", "Season2", 
                  "AB00_precip",0,
                  0,"AT07_precip"),2,4)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observtion model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model
# fit <- MARSS(y = dat_dep, model = mod_list,
#              control = list(maxit = 5000), method = "BFGS")
# see pg 5 in MARSS manual for notes on method BFGS vs method EM: EM algorithm gives more robust estimation for datasets replete with missing values and for high-dimensional models with various constraints. BFGS is faster and is good enough for some datasets. Typically, both should be tried.

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

# export model fit
saveRDS(fit, file = "data_working/marss_test_run/fit_120521_1state_AB00_EM.rds")
saveRDS(fit, file = "data_working/marss_test_run/fit_120521_2state_AB00_AT07_EM.rds")


### DIAGNOSES ###
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## check for hidden errors
fit[["errors"]]
# lots of errors when IND sites were included, none when these were removed

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
#est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)

#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_1state_AB00_EM_hessian.rds")
saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_2state_AB00_AT07_EM_hessian.rds")

CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)


### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00", "AT07", "GV01", "HO00", "MC06", "RG01", "RS02", "SP02","EFJ","RED","RSA","RSAW","SULF")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>% # rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "SP02_precip", 
                                             "EFJ_precip", "RED_precip", "RSA_precip", "RSAW_precip",
                                             "SULF_precip")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "NO3 MARSS modeling results - 12/5/2021") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))


## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- residuals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)
mod.null <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE)
# mod.null <- MARSS(y = dat_dep, model = mod_list_null,
#                    control = list(maxit = 5000), method = "BFGS")

bbmle::AICtab(fit, mod.null)

#           dAIC df
# fit       0.0 10
# mod.null 14.7 6 
# RESULT: covar model is better than null

### Plot response vars ###
par(mfrow=c(2,1),oma = c(0, 0, 2, 0))
plot(dat_dep[1,], type="o")
plot(dat_dep[2,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
forecast::Acf(resids$model.residuals[1,], main="1 model residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$state.residuals[1,], main="1 state residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$model.residuals[2,], main="2 model residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$state.residuals[2,], main="2 state residuals", na.action=na.pass, lag.max = 24)
mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)

### Do resids have temporal trend? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
plot(resids$model.residuals[1,], ylab="model residual", xlab="", main="1 model residuals")
abline(h=0)
plot(resids$state.residuals[1,], ylab="state residual", xlab="", main="1 state residuals")
abline(h=0)
plot(resids$model.residuals[2,], ylab="model residual", xlab="", main="2 model residuals")
abline(h=0)
plot(resids$state.residuals[2,], ylab="state residual", xlab="", main="2 state residuals")
abline(h=0)
mtext("Do resids have temporal trend?", outer = TRUE, cex = 1.5)

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
qqnorm(resids$model.residuals[1,], main="1 model residuals", pch=16, 
       xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[1,])[1]))
qqline(resids$model.residuals[1,])
qqnorm(resids$state.residuals[1,], main="1 state residuals", pch=16, 
       xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[1,])[1]))
qqline(resids$state.residuals[1,])
qqnorm(resids$model.residuals[2,], main="2 model residuals", pch=16,
       xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[2,])[1]))
qqline(resids$model.residuals[2,])
qqnorm(resids$state.residuals[2,], main="2 state residuals", pch=16,
       xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[2,])[1]))
qqline(resids$state.residuals[2,])
mtext("Are resids normal?", outer = TRUE, cex = 1.5)

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))


#### Scenario 1 : all catchments are separate states - 3 catchments #### 

# Pull out only response var
names(dat_no3_log)
# AB00, AT07
dat_dep <- t(dat_no3_log[,c(4,5,6)])
row.names(dat_dep)

# Make covariate inputs
dat_cov <- dat_no3_log[,c(2:3,
                          20,21,22)]
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)


#### make C matrix
# starting with HO00
CC <- matrix(list("Season1", "Season1", "Season1", 
                  "Season2", "Season2", "Season2", 
                  "AB00_precip",0, 0,
                  0,"AT07_precip",0,
                  0,0,"GV01_precip"),3,5)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observtion model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model
# fit <- MARSS(y = dat_dep, model = mod_list,
#              control = list(maxit = 5000), method = "BFGS")
# see pg 5 in MARSS manual for notes on method BFGS vs method EM: EM algorithm gives more robust estimation for datasets replete with missing values and for high-dimensional models with various constraints. BFGS is faster and is good enough for some datasets. Typically, both should be tried.

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

# export model fit
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_1state_AB00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_2state_AB00_AT07_EM.rds")
saveRDS(fit, file = "data_working/marss_test_run/fit_120521_3state_AB00_AT07_GV01_EM.rds")


### DIAGNOSES ###
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## check for hidden errors
fit[["errors"]]
# lots of errors when IND sites were included, none when these were removed

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
#est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)

#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_1state_AB00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_2state_AB00_AT07_EM_hessian.rds")
saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_3state_AB00_AT07_GV01_EM_hessian.rds")


CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)


### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00", "AT07", "GV01", "HO00", "MC06", "RG01", "RS02", "SP02","EFJ","RED","RSA","RSAW","SULF")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>% # rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "SP02_precip", 
                                             "EFJ_precip", "RED_precip", "RSA_precip", "RSAW_precip",
                                             "SULF_precip")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "NO3 MARSS modeling results - 12/5/2021") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))


## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- residuals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)
mod.null <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE)
# mod.null <- MARSS(y = dat_dep, model = mod_list_null,
#                    control = list(maxit = 5000), method = "BFGS")

bbmle::AICtab(fit, mod.null)

#           dAIC df
# fit       0.0 14
# mod.null 74.9 9 
# RESULT: covar model is better than null

### Plot response vars ###
par(mfrow=c(3,1),oma = c(0, 0, 2, 0))
plot(dat_dep[1,], type="o")
plot(dat_dep[2,], type="o")
plot(dat_dep[3,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(3,2),oma = c(0, 0, 2, 0))
forecast::Acf(resids$model.residuals[1,], main="1 model residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$state.residuals[1,], main="1 state residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$model.residuals[2,], main="2 model residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$state.residuals[2,], main="2 state residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$model.residuals[2,], main="3 model residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$state.residuals[2,], main="3 state residuals", na.action=na.pass, lag.max = 24)
mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)

### Are resids normal? ###
par(mfrow=c(3,2),oma = c(0, 0, 2, 0))
qqnorm(resids$model.residuals[1,], main="1 model residuals", pch=16, 
       xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[1,])[1]))
qqline(resids$model.residuals[1,])
qqnorm(resids$state.residuals[1,], main="1 state residuals", pch=16, 
       xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[1,])[1]))
qqline(resids$state.residuals[1,])
#
qqnorm(resids$model.residuals[2,], main="2 model residuals", pch=16,
       xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[2,])[1]))
qqline(resids$model.residuals[2,])
qqnorm(resids$state.residuals[2,], main="2 state residuals", pch=16,
       xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[2,])[1]))
qqline(resids$state.residuals[2,])
#
qqnorm(resids$model.residuals[3,], main="2 model residuals", pch=16,
       xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[2,])[1]))
qqline(resids$model.residuals[3,])
qqnorm(resids$state.residuals[3,], main="2 state residuals", pch=16,
       xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[2,])[1]))
qqline(resids$state.residuals[3,])
mtext("Are resids normal?", outer = TRUE, cex = 1.5)

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))


#### Scenario 1 : all catchments are separate states - 4 catchments #### 

# Pull out only response var
names(dat_no3_log)
# AB00, AT07
dat_dep <- t(dat_no3_log[,c(4,5,6,7)])
row.names(dat_dep)

# Make covariate inputs
dat_cov <- dat_no3_log[,c(2:3,
                          20,21,22,23)]
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)


#### make C matrix
# starting with HO00
CC <- matrix(list("Season1", "Season1", "Season1",  "Season1", 
                  "Season2", "Season2", "Season2",  "Season2", 
                  "AB00_precip",0, 0, 0,
                  0,"AT07_precip",0,0,
                  0,0,"GV01_precip",0,
                  0,0,0,"HO00_precip"),4,6)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observtion model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model
# fit <- MARSS(y = dat_dep, model = mod_list,
#              control = list(maxit = 5000), method = "BFGS")
# see pg 5 in MARSS manual for notes on method BFGS vs method EM: EM algorithm gives more robust estimation for datasets replete with missing values and for high-dimensional models with various constraints. BFGS is faster and is good enough for some datasets. Typically, both should be tried.

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

# export model fit
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_1state_AB00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_2state_AB00_AT07_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_3state_AB00_AT07_GV01_EM.rds")
saveRDS(fit, file = "data_working/marss_test_run/fit_120521_4state_AB00_AT07_GV01_HO00_EM.rds")



### DIAGNOSES ###
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## check for hidden errors
fit[["errors"]]
# lots of errors when IND sites were included, none when these were removed

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
#est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)

#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_1state_AB00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_2state_AB00_AT07_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_3state_AB00_AT07_GV01_EM_hessian.rds")
saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_4state_AB00_AT07_GV01_HO00_EM_hessian.rds")


CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)


### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00", "AT07", "GV01", "HO00", "MC06", "RG01", "RS02", "SP02","EFJ","RED","RSA","RSAW","SULF")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>% # rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "SP02_precip", 
                                             "EFJ_precip", "RED_precip", "RSA_precip", "RSAW_precip",
                                             "SULF_precip")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "NO3 MARSS modeling results - 12/5/2021") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))


## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- residuals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)
mod.null <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE)
# mod.null <- MARSS(y = dat_dep, model = mod_list_null,
#                    control = list(maxit = 5000), method = "BFGS")

bbmle::AICtab(fit, mod.null)

#           dAIC df
# fit       0   18
# mod.null 78   12
# RESULT: covar model is better than null

### Plot response vars ###
par(mfrow=c(4,1),oma = c(0, 0, 2, 0))
plot(dat_dep[1,], type="o")
plot(dat_dep[2,], type="o")
plot(dat_dep[3,], type="o")
plot(dat_dep[4,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
forecast::Acf(resids$model.residuals[1,], main="1 model residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$state.residuals[1,], main="1 state residuals", na.action=na.pass, lag.max = 24)
#
forecast::Acf(resids$model.residuals[2,], main="2 model residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$state.residuals[2,], main="2 state residuals", na.action=na.pass, lag.max = 24)
#
forecast::Acf(resids$model.residuals[3,], main="3 model residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$state.residuals[3,], main="3 state residuals", na.action=na.pass, lag.max = 24)
#
forecast::Acf(resids$model.residuals[4,], main="4 model residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$state.residuals[4,], main="4 state residuals", na.action=na.pass, lag.max = 24)
#
mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)

### Are resids normal? ###
par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
qqnorm(resids$model.residuals[1,], main="1 model residuals", pch=16, 
       xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[1,])[1]))
qqline(resids$model.residuals[1,])
qqnorm(resids$state.residuals[1,], main="1 state residuals", pch=16, 
       xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[1,])[1]))
qqline(resids$state.residuals[1,])
#
qqnorm(resids$model.residuals[2,], main="2 model residuals", pch=16,
       xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[2,])[1]))
qqline(resids$model.residuals[2,])
qqnorm(resids$state.residuals[2,], main="2 state residuals", pch=16,
       xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[2,])[1]))
qqline(resids$state.residuals[2,])
#
qqnorm(resids$model.residuals[3,], main="3 model residuals", pch=16,
       xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[2,])[1]))
qqline(resids$model.residuals[3,])
qqnorm(resids$state.residuals[3,], main="3 state residuals", pch=16,
       xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[2,])[1]))
qqline(resids$state.residuals[3,])
#
qqnorm(resids$model.residuals[4,], main="4 model residuals", pch=16,
       xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[2,])[1]))
qqline(resids$model.residuals[4,])
qqnorm(resids$state.residuals[4,], main="4 state residuals", pch=16,
       xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[2,])[1]))
qqline(resids$state.residuals[4,])
mtext("Are resids normal?", outer = TRUE, cex = 1.5)

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))


#### Scenario 1 : all catchments are separate states - 5 catchments #### 

# Pull out only response var
names(dat_no3_log)
# AB00, AT07
dat_dep <- t(dat_no3_log[,c(4,5,6,7,8)])
row.names(dat_dep)

# Make covariate inputs
dat_cov <- dat_no3_log[,c(2:3,
                          20,21,22,23,24)]
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)


#### make C matrix
# starting with HO00
CC <- matrix(list("Season1", "Season1", "Season1",  "Season1", "Season1", 
                  "Season2", "Season2", "Season2",  "Season2", "Season2", 
                  "AB00_precip",0, 0, 0,0,
                  0,"AT07_precip",0,0,0,
                  0,0,"GV01_precip",0,0,
                  0,0,0,"HO00_precip",0,
                  0,0,0,0,"MC06_precip"),5,7)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observtion model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model
# fit <- MARSS(y = dat_dep, model = mod_list,
#              control = list(maxit = 5000), method = "BFGS")
# see pg 5 in MARSS manual for notes on method BFGS vs method EM: EM algorithm gives more robust estimation for datasets replete with missing values and for high-dimensional models with various constraints. BFGS is faster and is good enough for some datasets. Typically, both should be tried.

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

# export model fit
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_1state_AB00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_2state_AB00_AT07_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_3state_AB00_AT07_GV01_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_4state_AB00_AT07_GV01_HO00_EM.rds")
saveRDS(fit, file = "data_working/marss_test_run/fit_120521_5state_AB00_AT07_GV01_HO00_MC06_EM.rds")



### DIAGNOSES ###
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## check for hidden errors
fit[["errors"]]
# lots of errors when IND sites were included, none when these were removed

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
#est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)

#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_1state_AB00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_2state_AB00_AT07_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_3state_AB00_AT07_GV01_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_4state_AB00_AT07_GV01_HO00_EM_hessian.rds")
saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_5state_AB00_AT07_GV01_HO00_MC06_EM_hessian.rds")

CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)


### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00", "AT07", "GV01", "HO00", "MC06", "RG01", "RS02", "SP02","EFJ","RED","RSA","RSAW","SULF")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>% # rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "SP02_precip", 
                                             "EFJ_precip", "RED_precip", "RSA_precip", "RSAW_precip",
                                             "SULF_precip")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "NO3 MARSS modeling results - 12/5/2021") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))


## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- residuals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)
mod.null <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE)
# mod.null <- MARSS(y = dat_dep, model = mod_list_null,
#                    control = list(maxit = 5000), method = "BFGS")

bbmle::AICtab(fit, mod.null)

#           dAIC df
# fit       0.0 22
# mod.null 80.3 15
# RESULT: covar model is better than null

### Plot response vars ###
par(mfrow=c(5,1),oma = c(0, 0, 2, 0))
plot(dat_dep[1,], type="o")
plot(dat_dep[2,], type="o")
plot(dat_dep[3,], type="o")
plot(dat_dep[4,], type="o")
plot(dat_dep[5,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:5)){
  forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:5)){
  qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"), 
         pch=16, 
         xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  qqline(resids$model.residuals[i,])
  qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16, 
         xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))


#### Scenario 1 : all catchments are separate states - 6 catchments #### 

# Pull out only response var
names(dat_no3_log)
# AB00, AT07
dat_dep <- t(dat_no3_log[,c(4,5,6,7,8,9)])
row.names(dat_dep)

# Make covariate inputs
dat_cov <- dat_no3_log[,c(2:3,
                          20,21,22,23,24,25)]
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)


#### make C matrix
# starting with HO00
CC <- matrix(list("Season1", "Season1", "Season1",  "Season1", "Season1","Season1", 
                  "Season2", "Season2", "Season2",  "Season2", "Season2", "Season2", 
                  "AB00_precip",0,0,0,0,0,
                  0,"AT07_precip",0,0,0,0,
                  0,0,"GV01_precip",0,0,0,
                  0,0,0,"HO00_precip",0,0,
                  0,0,0,0,"MC06_precip",0,
                  0,0,0,0,0,"RG01_precip"),6,8)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observtion model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model
# fit <- MARSS(y = dat_dep, model = mod_list,
#              control = list(maxit = 5000), method = "BFGS")
# see pg 5 in MARSS manual for notes on method BFGS vs method EM: EM algorithm gives more robust estimation for datasets replete with missing values and for high-dimensional models with various constraints. BFGS is faster and is good enough for some datasets. Typically, both should be tried.

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

# export model fit
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_1state_AB00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_2state_AB00_AT07_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_3state_AB00_AT07_GV01_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_4state_AB00_AT07_GV01_HO00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_5state_AB00_AT07_GV01_HO00_MC06_EM.rds")
saveRDS(fit, file = "data_working/marss_test_run/fit_120521_6state_AB00_AT07_GV01_HO00_MC06_RG01_EM.rds")



### DIAGNOSES ###
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## check for hidden errors
fit[["errors"]]
# lots of errors when IND sites were included, none when these were removed

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
#est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)

#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_1state_AB00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_2state_AB00_AT07_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_3state_AB00_AT07_GV01_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_4state_AB00_AT07_GV01_HO00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_5state_AB00_AT07_GV01_HO00_MC06_EM_hessian.rds")
saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_6state_AB00_AT07_GV01_HO00_MC06_RG01_EM_hessian.rds")

CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)


### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00", "AT07", "GV01", "HO00", "MC06", "RG01", "RS02", "SP02","EFJ","RED","RSA","RSAW","SULF")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>% # rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "SP02_precip", 
                                             "EFJ_precip", "RED_precip", "RSA_precip", "RSAW_precip",
                                             "SULF_precip")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "NO3 MARSS modeling results - 12/5/2021") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))


## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- residuals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)
mod.null <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE)
# mod.null <- MARSS(y = dat_dep, model = mod_list_null,
#                    control = list(maxit = 5000), method = "BFGS")

bbmle::AICtab(fit, mod.null)

#           dAIC df
# fit        0.0 26
# mod.null 102.3 18
# RESULT: covar model is better than null

### Plot response vars ###
par(mfrow=c(3,2),oma = c(0, 0, 2, 0))
plot(dat_dep[1,], type="o")
plot(dat_dep[2,], type="o")
plot(dat_dep[3,], type="o")
plot(dat_dep[4,], type="o")
plot(dat_dep[5,], type="o")
plot(dat_dep[6,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:6)){
  forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:6)){
  qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"), 
         pch=16, 
         xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  qqline(resids$model.residuals[i,])
  qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16, 
         xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))


#### Scenario 1 : all catchments are separate states - 7 catchments #### 

# Pull out only response var
names(dat_no3_log)
# AB00, AT07
dat_dep <- t(dat_no3_log[,c(4,5,6,7,8,9,10)])
row.names(dat_dep)

# Make covariate inputs
dat_cov <- dat_no3_log[,c(2:3,
                          20,21,22,23,24,25,26)]
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)


#### make C matrix
# starting with HO00
CC <- matrix(list("Season1", "Season1", "Season1",  "Season1", "Season1","Season1","Season1", 
                  "Season2", "Season2", "Season2",  "Season2", "Season2", "Season2", "Season2", 
                  "AB00_precip",0,0,0,0,0,0,
                  0,"AT07_precip",0,0,0,0,0,
                  0,0,"GV01_precip",0,0,0,0,
                  0,0,0,"HO00_precip",0,0,0,
                  0,0,0,0,"MC06_precip",0,0,
                  0,0,0,0,0,"RG01_precip",0,
                  0,0,0,0,0,0,"RS02_precip"),7,9)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observtion model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model
# fit <- MARSS(y = dat_dep, model = mod_list,
#              control = list(maxit = 5000), method = "BFGS")
# see pg 5 in MARSS manual for notes on method BFGS vs method EM: EM algorithm gives more robust estimation for datasets replete with missing values and for high-dimensional models with various constraints. BFGS is faster and is good enough for some datasets. Typically, both should be tried.

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

# export model fit
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_1state_AB00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_2state_AB00_AT07_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_3state_AB00_AT07_GV01_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_4state_AB00_AT07_GV01_HO00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_5state_AB00_AT07_GV01_HO00_MC06_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_6state_AB00_AT07_GV01_HO00_MC06_RG01_EM.rds")
saveRDS(fit, file = "data_working/marss_test_run/fit_120521_7state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EM.rds")



### DIAGNOSES ###
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## check for hidden errors
fit[["errors"]]
# lots of errors when IND sites were included, none when these were removed

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
#est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)

#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_1state_AB00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_2state_AB00_AT07_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_3state_AB00_AT07_GV01_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_4state_AB00_AT07_GV01_HO00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_5state_AB00_AT07_GV01_HO00_MC06_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_6state_AB00_AT07_GV01_HO00_MC06_RG01_EM_hessian.rds")
saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_7state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EM_hessian.rds")

CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)


### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00", "AT07", "GV01", "HO00", "MC06", "RG01", "RS02", "SP02","EFJ","RED","RSA","RSAW","SULF")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>% # rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "SP02_precip", 
                                             "EFJ_precip", "RED_precip", "RSA_precip", "RSAW_precip",
                                             "SULF_precip")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "NO3 MARSS modeling results - 12/5/2021") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))


## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- residuals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)
mod.null <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE)
# mod.null <- MARSS(y = dat_dep, model = mod_list_null,
#                    control = list(maxit = 5000), method = "BFGS")

bbmle::AICtab(fit, mod.null)

#           dAIC df
# fit        0.0 30
# mod.null 155.3 21
# RESULT: covar model is better than null

### Plot response vars ###
par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
plot(dat_dep[1,], type="o")
plot(dat_dep[2,], type="o")
plot(dat_dep[3,], type="o")
plot(dat_dep[4,], type="o")
plot(dat_dep[5,], type="o")
plot(dat_dep[6,], type="o")
plot(dat_dep[7,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:7)){
  forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:7)){
  qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"), 
         pch=16, 
         xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  qqline(resids$model.residuals[i,])
  qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16, 
         xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))


#### Scenario 1 : all catchments are separate states - 8 catchments ***ADDING SP02 IS WHERE MODEL FIT BREAKS DOWN*** #### 

# Pull out only response var
names(dat_no3_log)
# AB00, AT07
dat_dep <- t(dat_no3_log[,c(4,5,6,7,8,9,10,11)])
row.names(dat_dep)

# Make covariate inputs
dat_cov <- dat_no3_log[,c(2:3,
                          20,21,22,23,24,25,26,27)]
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)


#### make C matrix
# starting with HO00
CC <- matrix(list("Season1", "Season1", "Season1",  "Season1", "Season1","Season1","Season1", "Season1", 
                  "Season2", "Season2", "Season2",  "Season2", "Season2", "Season2", "Season2", "Season2", 
                  "AB00_precip",0,0,0,0,0,0,0,
                  0,"AT07_precip",0,0,0,0,0,0,
                  0,0,"GV01_precip",0,0,0,0,0,
                  0,0,0,"HO00_precip",0,0,0,0,
                  0,0,0,0,"MC06_precip",0,0,0,
                  0,0,0,0,0,"RG01_precip",0,0,
                  0,0,0,0,0,0,"RS02_precip",0,
                  0,0,0,0,0,0,0,"SP02_precip"),8,10)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observtion model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model
# fit <- MARSS(y = dat_dep, model = mod_list,
#              control = list(maxit = 5000), method = "BFGS")
# see pg 5 in MARSS manual for notes on method BFGS vs method EM: EM algorithm gives more robust estimation for datasets replete with missing values and for high-dimensional models with various constraints. BFGS is faster and is good enough for some datasets. Typically, both should be tried.

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

# export model fit
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_1state_AB00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_2state_AB00_AT07_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_3state_AB00_AT07_GV01_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_4state_AB00_AT07_GV01_HO00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_5state_AB00_AT07_GV01_HO00_MC06_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_6state_AB00_AT07_GV01_HO00_MC06_RG01_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_7state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EM.rds")
saveRDS(fit, file = "data_working/marss_test_run/fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_SP02_EM.rds")



### DIAGNOSES ###
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## check for hidden errors
fit[["errors"]]
# lots of errors when IND sites were included, none when these were removed

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
#est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)
#***************** SP02 IS WHERE THE HESSIAN METHOD BREAKS DOWN!!! *********************

#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_1state_AB00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_2state_AB00_AT07_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_3state_AB00_AT07_GV01_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_4state_AB00_AT07_GV01_HO00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_5state_AB00_AT07_GV01_HO00_MC06_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_6state_AB00_AT07_GV01_HO00_MC06_RG01_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_7state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_SP02_EM_hessian.rds") ****** NOT SAVED *******


#### Scenario 1 : all catchments are separate states - 8 catchments ***ADDING EFJ WORKS BUT THE DATA & RESIDS NOT GREAT #### 

# Pull out only response var
names(dat_no3_log)
# AB00, AT07
dat_dep <- t(dat_no3_log[,c(4,5,6,7,8,9,10,12)])
row.names(dat_dep)

# Make covariate inputs
dat_cov <- dat_no3_log[,c(2:3,
                          20,21,22,23,24,25,26,28)]
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)


#### make C matrix
# starting with HO00
CC <- matrix(list("Season1", "Season1", "Season1",  "Season1", "Season1","Season1","Season1", "Season1", 
                  "Season2", "Season2", "Season2",  "Season2", "Season2", "Season2", "Season2", "Season2", 
                  "AB00_precip",0,0,0,0,0,0,0,
                  0,"AT07_precip",0,0,0,0,0,0,
                  0,0,"GV01_precip",0,0,0,0,0,
                  0,0,0,"HO00_precip",0,0,0,0,
                  0,0,0,0,"MC06_precip",0,0,0,
                  0,0,0,0,0,"RG01_precip",0,0,
                  0,0,0,0,0,0,"RS02_precip",0,
                  0,0,0,0,0,0,0,"EFJ_precip"),8,10)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observtion model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model
# fit <- MARSS(y = dat_dep, model = mod_list,
#              control = list(maxit = 5000), method = "BFGS")
# see pg 5 in MARSS manual for notes on method BFGS vs method EM: EM algorithm gives more robust estimation for datasets replete with missing values and for high-dimensional models with various constraints. BFGS is faster and is good enough for some datasets. Typically, both should be tried.

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

# export model fit
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_1state_AB00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_2state_AB00_AT07_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_3state_AB00_AT07_GV01_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_4state_AB00_AT07_GV01_HO00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_5state_AB00_AT07_GV01_HO00_MC06_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_6state_AB00_AT07_GV01_HO00_MC06_RG01_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_7state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_SP02_EM.rds")
saveRDS(fit, file = "data_working/marss_test_run/fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_EM.rds")



### DIAGNOSES ###
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## check for hidden errors
fit[["errors"]]
# lots of errors when IND sites were included, none when these were removed

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
#est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)
#***************** SP02 IS WHERE THE HESSIAN METHOD BREAKS DOWN!!! *********************

#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_1state_AB00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_2state_AB00_AT07_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_3state_AB00_AT07_GV01_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_4state_AB00_AT07_GV01_HO00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_5state_AB00_AT07_GV01_HO00_MC06_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_6state_AB00_AT07_GV01_HO00_MC06_RG01_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_7state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_SP02_EM_hessian.rds") ****** NOT SAVED *******
saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_EM_hessian.rds")

CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)


### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00", "AT07", "GV01", "HO00", "MC06", "RG01", "RS02", "SP02","EFJ","RED","RSA","RSAW","SULF")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>% # rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "SP02_precip", 
                                             "EFJ_precip", "RED_precip", "RSA_precip", "RSAW_precip",
                                             "SULF_precip")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "NO3 MARSS modeling results - 12/5/2021") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))


## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- residuals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)
mod.null <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE)
# mod.null <- MARSS(y = dat_dep, model = mod_list_null,
#                    control = list(maxit = 5000), method = "BFGS")

bbmle::AICtab(fit, mod.null)

#           dAIC df
# fit        0.0 34
# mod.null 154.2 24
# RESULT: covar model is better than null

### Plot response vars ###
par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
plot(dat_dep[1,], type="o")
plot(dat_dep[2,], type="o")
plot(dat_dep[3,], type="o")
plot(dat_dep[4,], type="o")
plot(dat_dep[5,], type="o")
plot(dat_dep[6,], type="o")
plot(dat_dep[7,], type="o")
plot(dat_dep[8,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:8)){
  forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:8)){
  qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"), 
         pch=16, 
         xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  qqline(resids$model.residuals[i,])
  qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16, 
         xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

#************************** EFJ has worst normality in residuals than most - the data has very little variation in it

#### Scenario 1 : all catchments are separate states - 8 catchments ***ADDING RED WORKS BUT THE DATA & RESIDS NOT GREAT #### 

# Pull out only response var
names(dat_no3_log)
# AB00, AT07
dat_dep <- t(dat_no3_log[,c(4,5,6,7,8,9,10,16)])
row.names(dat_dep)

# Make covariate inputs
dat_cov <- dat_no3_log[,c(2:3,
                          20,21,22,23,24,25,26,32)]
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)


#### make C matrix
# starting with HO00
CC <- matrix(list("Season1", "Season1", "Season1",  "Season1", "Season1","Season1","Season1", "Season1", 
                  "Season2", "Season2", "Season2",  "Season2", "Season2", "Season2", "Season2", "Season2", 
                  "AB00_precip",0,0,0,0,0,0,0,
                  0,"AT07_precip",0,0,0,0,0,0,
                  0,0,"GV01_precip",0,0,0,0,0,
                  0,0,0,"HO00_precip",0,0,0,0,
                  0,0,0,0,"MC06_precip",0,0,0,
                  0,0,0,0,0,"RG01_precip",0,0,
                  0,0,0,0,0,0,"RS02_precip",0,
                  0,0,0,0,0,0,0,"RED_precip"),8,10)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observtion model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model
# fit <- MARSS(y = dat_dep, model = mod_list,
#              control = list(maxit = 5000), method = "BFGS")
# see pg 5 in MARSS manual for notes on method BFGS vs method EM: EM algorithm gives more robust estimation for datasets replete with missing values and for high-dimensional models with various constraints. BFGS is faster and is good enough for some datasets. Typically, both should be tried.

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

# export model fit
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_1state_AB00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_2state_AB00_AT07_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_3state_AB00_AT07_GV01_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_4state_AB00_AT07_GV01_HO00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_5state_AB00_AT07_GV01_HO00_MC06_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_6state_AB00_AT07_GV01_HO00_MC06_RG01_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_7state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_SP02_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_EM.rds")
saveRDS(fit, file = "data_working/marss_test_run/fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_RED_EM.rds")


### DIAGNOSES ###
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## check for hidden errors
fit[["errors"]]
# lots of errors when IND sites were included, none when these were removed

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
#est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)
#***************** SP02 IS WHERE THE HESSIAN METHOD BREAKS DOWN!!! *********************

#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_1state_AB00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_2state_AB00_AT07_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_3state_AB00_AT07_GV01_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_4state_AB00_AT07_GV01_HO00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_5state_AB00_AT07_GV01_HO00_MC06_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_6state_AB00_AT07_GV01_HO00_MC06_RG01_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_7state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_SP02_EM_hessian.rds") ****** NOT SAVED *******
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_EM_hessian.rds")
saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_RED_EM_hessian.rds")

CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)


### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00", "AT07", "GV01", "HO00", "MC06", "RG01", "RS02", "SP02","EFJ","RED","RSA","RSAW","SULF")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>% # rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "SP02_precip", 
                                             "EFJ_precip", "RED_precip", "RSA_precip", "RSAW_precip",
                                             "SULF_precip")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "NO3 MARSS modeling results - 12/5/2021") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))


## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- residuals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)
mod.null <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE)
# mod.null <- MARSS(y = dat_dep, model = mod_list_null,
#                    control = list(maxit = 5000), method = "BFGS")

bbmle::AICtab(fit, mod.null)

#           dAIC df
# fit        0.0 34
# mod.null 156.3 24
# RESULT: covar model is better than null

### Plot response vars ###
par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
plot(dat_dep[1,], type="o")
plot(dat_dep[2,], type="o")
plot(dat_dep[3,], type="o")
plot(dat_dep[4,], type="o")
plot(dat_dep[5,], type="o")
plot(dat_dep[6,], type="o")
plot(dat_dep[7,], type="o")
plot(dat_dep[8,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:8)){
  forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:8)){
  qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"), 
         pch=16, 
         xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  qqline(resids$model.residuals[i,])
  qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16, 
         xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

#************************** RED has worst normality in residuals than most - the data has very little variation in it

#### Scenario 1 : all catchments are separate states - 8 catchments ***ADDING SULF WORKS BUT THE DATA & RESIDS NOT GREAT #### 

# Pull out only response var
names(dat_no3_log)
# AB00, AT07
dat_dep <- t(dat_no3_log[,c(4,5,6,7,8,9,10,19)])
row.names(dat_dep)

# Make covariate inputs
dat_cov <- dat_no3_log[,c(2:3,
                          20,21,22,23,24,25,26,
                          35)]
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)


#### make C matrix
# starting with HO00
CC <- matrix(list("Season1", "Season1", "Season1",  "Season1", "Season1","Season1","Season1", "Season1", 
                  "Season2", "Season2", "Season2",  "Season2", "Season2", "Season2", "Season2", "Season2", 
                  "AB00_precip",0,0,0,0,0,0,0,
                  0,"AT07_precip",0,0,0,0,0,0,
                  0,0,"GV01_precip",0,0,0,0,0,
                  0,0,0,"HO00_precip",0,0,0,0,
                  0,0,0,0,"MC06_precip",0,0,0,
                  0,0,0,0,0,"RG01_precip",0,0,
                  0,0,0,0,0,0,"RS02_precip",0,
                  0,0,0,0,0,0,0,"SULF_precip"),8,10)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observtion model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model
# fit <- MARSS(y = dat_dep, model = mod_list,
#              control = list(maxit = 5000), method = "BFGS")
# see pg 5 in MARSS manual for notes on method BFGS vs method EM: EM algorithm gives more robust estimation for datasets replete with missing values and for high-dimensional models with various constraints. BFGS is faster and is good enough for some datasets. Typically, both should be tried.

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

# export model fit
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_1state_AB00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_2state_AB00_AT07_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_3state_AB00_AT07_GV01_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_4state_AB00_AT07_GV01_HO00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_5state_AB00_AT07_GV01_HO00_MC06_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_6state_AB00_AT07_GV01_HO00_MC06_RG01_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_7state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_SP02_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_EM.rds")
saveRDS(fit, file = "data_working/marss_test_run/fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_SULF_EM.rds")


### DIAGNOSES ###
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## check for hidden errors
fit[["errors"]]
# lots of errors when IND sites were included, none when these were removed

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
#est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)
#***************** SP02 IS WHERE THE HESSIAN METHOD BREAKS DOWN!!! *********************

#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_1state_AB00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_2state_AB00_AT07_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_3state_AB00_AT07_GV01_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_4state_AB00_AT07_GV01_HO00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_5state_AB00_AT07_GV01_HO00_MC06_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_6state_AB00_AT07_GV01_HO00_MC06_RG01_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_7state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_SP02_EM_hessian.rds") ****** NOT SAVED *******
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_EM_hessian.rds")
saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_SULF_EM_hessian.rds")

CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)


### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00", "AT07", "GV01", "HO00", "MC06", "RG01", "RS02", "SP02","EFJ","RED","RSA","RSAW","SULF")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>% # rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "SP02_precip", 
                                             "EFJ_precip", "RED_precip", "RSA_precip", "RSAW_precip",
                                             "SULF_precip")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "NO3 MARSS modeling results - 12/5/2021") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))


## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- residuals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)
mod.null <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE)
# mod.null <- MARSS(y = dat_dep, model = mod_list_null,
#                    control = list(maxit = 5000), method = "BFGS")

bbmle::AICtab(fit, mod.null)

#           dAIC df
# fit        0.0 34
# mod.null 157.1 24
# RESULT: covar model is better than null

### Plot response vars ###
par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
plot(dat_dep[1,], type="o")
plot(dat_dep[2,], type="o")
plot(dat_dep[3,], type="o")
plot(dat_dep[4,], type="o")
plot(dat_dep[5,], type="o")
plot(dat_dep[6,], type="o")
plot(dat_dep[7,], type="o")
plot(dat_dep[8,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:8)){
  forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:3,5:8)){
  qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"), 
         pch=16, 
         xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  qqline(resids$model.residuals[i,])
  qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16, 
         xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

#************************** SULF has worst normality in residuals than most - the data has very little variation in it

#### Scenario 1 : all catchments are separate states - 12 catchments ***Runs with BFGS+priors, but NM data&resids are shit  #### 

# Pull out only response var
names(dat_no3_log)
# AB00,AT07,GV01,HO00,MC06,RG01,RS02,
#  EFJ, RED, RSA,RSAW,SULF
dat_dep <- t(dat_no3_log[,c(4,5,6,7,8,9,10,12,16,17,18,19)])
row.names(dat_dep)

# Make covariate inputs
dat_cov <- dat_no3_log[,c(2:3,
                          20,21,22,23,24,25,26,28,32,33,34,35)]
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)


#### make C matrix
# starting with HO00
CC <- matrix(list("Season1", "Season1", "Season1",  "Season1", "Season1","Season1","Season1", "Season1", "Season1", "Season1", "Season1", "Season1", 
                  "Season2", "Season2", "Season2",  "Season2", "Season2", "Season2", "Season2", "Season2","Season2","Season2","Season2","Season2", 
                  "AB00_precip",0,0,0,0,0,0,0,0,0,0,0,
                  0,"AT07_precip",0,0,0,0,0,0,0,0,0,0,
                  0,0,"GV01_precip",0,0,0,0,0,0,0,0,0,
                  0,0,0,"HO00_precip",0,0,0,0,0,0,0,0,
                  0,0,0,0,"MC06_precip",0,0,0,0,0,0,0,
                  0,0,0,0,0,"RG01_precip",0,0,0,0,0,0,
                  0,0,0,0,0,0,"RS02_precip",0,0,0,0,0,
                  0,0,0,0,0,0,0,"EFJ_precip",0,0,0,0,
                  0,0,0,0,0,0,0,0,"RED_precip",0,0,0,
                  0,0,0,0,0,0,0,0,0,"RSA_precip",0,0,
                  0,0,0,0,0,0,0,0,0,0,"RSAW_precip",0,
                  0,0,0,0,0,0,0,0,0,0,0,"SULF_precip"),12,14)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observtion model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model
kemfit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_1state_AB00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_2state_AB00_AT07_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_3state_AB00_AT07_GV01_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_4state_AB00_AT07_GV01_HO00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_5state_AB00_AT07_GV01_HO00_MC06_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_6state_AB00_AT07_GV01_HO00_MC06_RG01_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_7state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_SP02_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_RED_EM.rds")
saveRDS(fit, file = "data_working/marss_test_run/fit_120521_12state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_RED_RSA_RSAW_SULF_BFGS.rds")


### DIAGNOSES ###
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## check for hidden errors
fit[["errors"]]
# lots of errors when IND sites were included, none when these were removed

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
#est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)
#***************** SP02 IS WHERE THE HESSIAN METHOD BREAKS DOWN!!! *********************

#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_1state_AB00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_2state_AB00_AT07_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_3state_AB00_AT07_GV01_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_4state_AB00_AT07_GV01_HO00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_5state_AB00_AT07_GV01_HO00_MC06_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_6state_AB00_AT07_GV01_HO00_MC06_RG01_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_7state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_SP02_EM_hessian.rds") ****** NOT SAVED *******
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_RED_EM_hessian.rds")
saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_12state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_RED_RSA_RSAW_SULF_BFGS.rds")

CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)


### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00", "AT07", "GV01", "HO00", "MC06", "RG01", "RS02", "SP02","EFJ","RED","RSA","RSAW","SULF")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>% # rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "SP02_precip", 
                                             "EFJ_precip", "RED_precip", "RSA_precip", "RSAW_precip",
                                             "SULF_precip")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "NO3 MARSS modeling results - 12/5/2021") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))


## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- residuals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)
null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"
null.fit <- MARSS(y = dat_dep, model = mod_list_null,
             control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#           dAIC df
# fit        0.0 50
# null.fit 162.7 36
# RESULT: covar model is better than null

### Plot response vars ###
par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
plot(dat_dep[1,], type="o")
plot(dat_dep[2,], type="o")
plot(dat_dep[3,], type="o")
plot(dat_dep[4,], type="o")
plot(dat_dep[5,], type="o")
plot(dat_dep[6,], type="o")
plot(dat_dep[7,], type="o")
plot(dat_dep[8,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  #forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  # qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"), 
  #        pch=16, 
  #        xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  # qqline(resids$model.residuals[i,])
  qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16, 
         xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

#************************** all the NM sites have terrible residuals, likely bc there is so little variation in the data and so many below detection values

#### Scenario 1 : (SB only) all catchments are separate states  #### 

# Mimicking NH4 model above, and redoing the dataset because all NO3 models above have no fire variables

# first, subset and pivot the NO3 data
dat_no3 <- dat_new22 %>%
  select(site, index, Season1, Season2, 
         mean_no3_uM, cumulative_precip_mm,
         AB00_Tea, AB00_Jesusita, AT07_Jesusita, GV01_Gaviota, HO00_Gaviota, HO00_Sherpa, MC06_Tea, MC06_Jesusita, RG01_Gaviota, RG01_Sherpa, RS02_Tea, RS02_Jesusita, RED_Thompson, EFJ_Thompson, EFJ_Conchas, RSAW_Thompson, RSAW_Conchas, RSA_Conchas) %>% 
  pivot_wider(names_from = site, values_from = c(mean_no3_uM, cumulative_precip_mm, AB00_Tea, AB00_Jesusita, AT07_Jesusita, GV01_Gaviota, HO00_Gaviota, HO00_Sherpa, MC06_Tea, MC06_Jesusita, RG01_Gaviota, RG01_Sherpa, RS02_Tea, RS02_Jesusita, RED_Thompson, EFJ_Thompson, EFJ_Conchas, RSAW_Thompson, RSAW_Conchas, RSA_Conchas)) %>%
  select(index, Season1, Season2, 
         mean_no3_uM_AB00, mean_no3_uM_AT07, mean_no3_uM_GV01, mean_no3_uM_HO00, mean_no3_uM_MC06, mean_no3_uM_RG01, mean_no3_uM_RS02, mean_no3_uM_EFJ, mean_no3_uM_RED, mean_no3_uM_RSA, mean_no3_uM_RSAW,
         cumulative_precip_mm_AB00, cumulative_precip_mm_AT07, cumulative_precip_mm_GV01, cumulative_precip_mm_HO00, cumulative_precip_mm_MC06, cumulative_precip_mm_RG01, cumulative_precip_mm_RS02, cumulative_precip_mm_EFJ, cumulative_precip_mm_RED, cumulative_precip_mm_RSA, cumulative_precip_mm_RSAW, 
         AB00_Tea_AB00, AB00_Jesusita_AB00, AT07_Jesusita_AT07, GV01_Gaviota_GV01, HO00_Gaviota_HO00, HO00_Sherpa_HO00, MC06_Tea_MC06, MC06_Jesusita_MC06, RG01_Gaviota_RG01, RG01_Sherpa_RG01, RS02_Tea_RS02, RS02_Jesusita_RS02, EFJ_Thompson_EFJ, EFJ_Conchas_EFJ, RED_Thompson_RED, RSA_Conchas_RSA, RSAW_Thompson_RSAW, RSAW_Conchas_RSAW)

# Remove the NaNs and replace with NAs.
dat_no3[is.nan(dat_no3)] <- NA

# log and scale transform response var
names(dat_no3)
dat_no3_log <- dat_no3
dat_no3_log[,4:14] <- log10(dat_no3_log[,4:14])
dat_no3_log[,4:14] <- scale(dat_no3_log[,4:14])
sum(is.nan(dat_no3_log[,4:14])) # 0 yay!!
sum(is.na(dat_no3_log[,4:14])) # 711
range(dat_no3_log[,4:14], na.rm = T)

# Pull out only NO3 data
names(dat_no3_log)
# Remaining sites: AB00, AT07, GV01, HO00, MC06, RG01, RS02
# 7 sites total with 166 observations each
dat_dep <- t(dat_no3_log[,c(4:10)])
row.names(dat_dep)

# Make covariate inputs
# 21 covariates total with 166 observations each
dat_cov <- dat_no3_log[,c(2:3, # seasons
                          15:21, # precip
                          26:37)] # fires
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)

# make C matrix
# this matrix controls what covars predict what response vars; 
# in contrast to unconstrained: separately estimates ALL correlations of 
# predictor vars with x, i.e. how predictor vars drive ts dynamics. 
# This structure is ok given that covars are unique to each site, but it 
# would not be ok if there is a mix of shared and unique covars.
# 7 response variables - so 7 rows in the matrix
# 21 covariates in total - so 21 columns in the matrix

# without short ts sites
CC <- matrix(list(# season 1 covariate
  "Season1", "Season1", "Season1", "Season1", 
  "Season1", "Season1", "Season1",
  # season 2 covariate
  "Season2", "Season2", "Season2", "Season2", 
  "Season2", "Season2", "Season2",
  # precipitation covariates
  "AB00_precip", 0, 0, 0, 0, 0, 0,
  0, "AT07_precip", 0, 0, 0, 0, 0, 
  0, 0, "GV01_precip", 0, 0, 0, 0,
  0, 0, 0, "HO00_precip", 0, 0, 0,
  0, 0, 0, 0, "MC06_precip", 0, 0,
  0, 0, 0, 0, 0, "RG01_precip", 0,
  0, 0, 0, 0, 0, 0, "RS02_precip",
  # fire covariates
  "AB00_Tea_fire", 0, 0, 0, 0, 0, 0,
  "AB00_Jesusita_fire", 0, 0, 0, 0, 0, 0,
  0, "AT07_Jesusita_fire", 0, 0, 0, 0, 0,
  0, 0, "GV01_Gaviota_fire", 0, 0, 0, 0,
  0, 0, 0, "HO00_Gaviota_fire", 0, 0, 0,
  0, 0, 0, "HO00_Sherpa_fire", 0, 0, 0,
  0, 0, 0, 0, "MC06_Tea_fire", 0, 0,
  0, 0, 0, 0, "MC06_Jesusita_fire", 0, 0,
  0, 0, 0, 0, 0, "RG01_Gaviota_fire", 0,
  0, 0, 0, 0, 0, "RG01_Sherpa_fire", 0,
  0, 0, 0, 0, 0, 0, "RS02_Tea_fire",
  0, 0, 0, 0, 0, 0, "RS02_Jesusita_fire"),
  7, 21)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit = 100, allow.degen = TRUE, trace = 1),
                fit = TRUE) # used to get some priors

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", 
             inits=kemfit$par) # actual BFGS method model

# fit_em <- MARSS(y = dat_dep, model = mod_list,
#                 control = list(maxit = 2000, allow.degen = TRUE, trace = 1),
#                 fit = TRUE) # actual EM method model

# see pg 5 in MARSS manual for notes on method BFGS vs method EM: EM algorithm gives more robust estimation for datasets replete with missing values and for high-dimensional models with various constraints. BFGS is faster and is good enough for some datasets. Typically, both should be tried.

# export model fit
#saveRDS(fit, file = "data_working/marss_test_run/fit_032322_7state_no3_AB00_AT07_GV01_HO00_MC06_RG01_RS02_mBFGS.rds")

### DIAGNOSES ###

## check for hidden errors
fit[["errors"]]
# NULL- Yay!!

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)

# No error messages!
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_032322_7state_no3_AB00_AT07_GV01_HO00_MC06_RG01_RS02_mBFGS.rds")

CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)

### Plot Results for All Sites ###

# First, create dataset of all outputs
my_list <- c("AB00", "AT07", "GV01", "HO00", 
             "MC06", "RG01", "RS02")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>% # rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip",
                                             "AB00_Tea_fire","AB00_Jesusita_fire",
                                             "AT07_Jesusita_fire",
                                             "GV01_Gaviota_fire",
                                             "HO00_Gaviota_fire","HO00_Sherpa_fire",
                                             "MC06_Tea_fire", "MC06_Jesusita_fire",
                                             "RG01_Gaviota_fire", "RG01_Sherpa_fire",
                                             "RS02_Tea_fire", "RS02_Jesusita_fire")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "NO3 MARSS modeling results SB only - 03/23/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))

# adding some formatting to match cond figure output
CIs_fit_ed2 = CIs_fit_ed
CIs_fit_ed2$region = c(rep("Coastal California",33))

(RESULTS_ALL2 <- ggplot(CIs_fit_ed2, aes(Parameter, Est., color=region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=.7) +
    geom_point(position=position_dodge(width=0.3), size=5) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Nitrate (NO3) MARSS modeling results SB only - 3/23/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm"),# need to play with margins to make it all fit
          legend.position = "none") +
    facet_wrap(vars(region, Site), scales = "free"))

# ggsave("figures/MARSS_7states_SBonly_no3_precip_fire_032322.png",
#        width = 25,
#        height = 15,
#        units = "cm")

## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- MARSSresiduals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)

null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, 
                                    allow.degen=TRUE, trace=1), 
                     fit=TRUE) # default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), 
                  method = "BFGS", 
                  inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

# dAIC - delta AIC
# 0.0 = always the value for the lowest model AIC
#          dAIC  df
# fit       0.0  42
# null.fit 150.1 21
# RESULT: covar model is better than null - yay!!

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  # forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}

# Model residuals look alright - only errant lags at certain sites.
# State residuals look fine - no major patterns again, although should look at results for state 5/MC06 (crossing threshold at 3, 4, 11, and 12 month lags).

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"),
         pch=16,
         xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  qqline(resids$model.residuals[i,])
  # qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16,
  #        xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  # qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# State residuals look ok....state 2/AT07 does have a step in it again.
# Model residuals look TERRIBLE...all flat lines again.

# reset plotting window
dev.off()

#### po4 ####

dat_po4 <- dat_agu %>%
  select(site, index, Season1, Season2, 
         mean_po4_uM, cumulative_precip_mm) %>% 
  pivot_wider(names_from = site, values_from = c(mean_po4_uM, cumulative_precip_mm)) %>%
  select(index, Season1, Season2, 
         mean_po4_uM_AB00, mean_po4_uM_AT07, mean_po4_uM_GV01, mean_po4_uM_HO00, mean_po4_uM_MC06, mean_po4_uM_RG01, mean_po4_uM_RS02, mean_po4_uM_SP02, mean_po4_uM_EFJ, mean_po4_uM_IND, mean_po4_uM_IND_AB, mean_po4_uM_IND_BB, mean_po4_uM_RED, mean_po4_uM_RSA, mean_po4_uM_RSAW, mean_po4_uM_SULF,
         cumulative_precip_mm_AB00, cumulative_precip_mm_AT07, cumulative_precip_mm_GV01, cumulative_precip_mm_HO00, cumulative_precip_mm_MC06, cumulative_precip_mm_RG01, cumulative_precip_mm_RS02, cumulative_precip_mm_SP02, cumulative_precip_mm_EFJ, cumulative_precip_mm_IND, cumulative_precip_mm_IND_AB, cumulative_precip_mm_IND_BB, cumulative_precip_mm_RED, cumulative_precip_mm_RSA, cumulative_precip_mm_RSAW, cumulative_precip_mm_SULF)

dat_po4[is.nan(dat_po4)] <- NA

# log and scale transform response var
names(dat_po4)
dat_po4_log = dat_po4
dat_po4_log[,4:19] = log10(dat_po4_log[,4:19])
dat_po4_log[,4:19] = scale(dat_po4_log[,4:19])
sum(is.nan(dat_po4_log[,4:19]))
sum(is.na(dat_po4_log[,4:19]))
range(dat_po4_log[,4:19], na.rm = T)

### Plot response vars ###
par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
plot(dat_po4_log$mean_po4_uM_AB00, type="o")
plot(dat_po4_log$mean_po4_uM_AT07, type="o")
plot(dat_po4_log$mean_po4_uM_GV01, type="o")
plot(dat_po4_log$mean_po4_uM_HO00, type="o")
plot(dat_po4_log$mean_po4_uM_MC06, type="o")
plot(dat_po4_log$mean_po4_uM_RG01, type="o")
plot(dat_po4_log$mean_po4_uM_RS02, type="o")
plot(dat_po4_log$mean_po4_uM_SP02, type="o")
plot(dat_po4_log$mean_po4_uM_EFJ, type="o")
plot(dat_po4_log$mean_po4_uM_IND, type="o")
plot(dat_po4_log$mean_po4_uM_IND_AB, type="o")
plot(dat_po4_log$mean_po4_uM_IND_BB, type="o")
plot(dat_po4_log$mean_po4_uM_RED, type="o")
plot(dat_po4_log$mean_po4_uM_RSA, type="o")
plot(dat_po4_log$mean_po4_uM_RSAW, type="o")
plot(dat_po4_log$mean_po4_uM_SULF, type="o")
# there are only a few data points in most of the NM sites

#### mean_cond_uScm ####

dat_cond <- dat_new22 %>%
  select(site, index, Season1, Season2, 
         mean_cond_uScm, cumulative_precip_mm,
         AB00_Tea, AB00_Jesusita, AT07_Jesusita, GV01_Gaviota, HO00_Gaviota, HO00_Sherpa, MC06_Tea, MC06_Jesusita, RG01_Gaviota, RG01_Sherpa, RS02_Tea, RS02_Jesusita, SP02_Gap, RED_Thompson, EFJ_Thompson, EFJ_Conchas, RSAW_Thompson, RSAW_Conchas, RSA_Conchas, IND_BB_Conchas, SULF_Thompson) %>% 
  pivot_wider(names_from = site, values_from = c(mean_cond_uScm, cumulative_precip_mm, AB00_Tea, AB00_Jesusita, AT07_Jesusita, GV01_Gaviota, HO00_Gaviota, HO00_Sherpa, MC06_Tea, MC06_Jesusita, RG01_Gaviota, RG01_Sherpa, RS02_Tea, RS02_Jesusita, SP02_Gap, RED_Thompson, EFJ_Thompson, EFJ_Conchas, RSAW_Thompson, RSAW_Conchas, RSA_Conchas, IND_BB_Conchas, SULF_Thompson)) %>%
  select(index, Season1, Season2, 
         mean_cond_uScm_AB00, mean_cond_uScm_AT07, mean_cond_uScm_GV01, mean_cond_uScm_HO00, mean_cond_uScm_MC06, mean_cond_uScm_RG01, mean_cond_uScm_RS02, mean_cond_uScm_SP02, mean_cond_uScm_EFJ, mean_cond_uScm_IND, mean_cond_uScm_IND_BB, mean_cond_uScm_RED, mean_cond_uScm_RSA, mean_cond_uScm_RSAW, mean_cond_uScm_SULF,
         cumulative_precip_mm_AB00, cumulative_precip_mm_AT07, cumulative_precip_mm_GV01, cumulative_precip_mm_HO00, cumulative_precip_mm_MC06, cumulative_precip_mm_RG01, cumulative_precip_mm_RS02, cumulative_precip_mm_SP02, cumulative_precip_mm_EFJ, cumulative_precip_mm_IND, cumulative_precip_mm_IND_BB, cumulative_precip_mm_RED, cumulative_precip_mm_RSA, cumulative_precip_mm_RSAW, cumulative_precip_mm_SULF,
         AB00_Tea_AB00, AB00_Jesusita_AB00, AT07_Jesusita_AT07, GV01_Gaviota_GV01, HO00_Gaviota_HO00, HO00_Sherpa_HO00, MC06_Tea_MC06, MC06_Jesusita_MC06, RG01_Gaviota_RG01, RG01_Sherpa_RG01, RS02_Tea_RS02, RS02_Jesusita_RS02, SP02_Gap_SP02, RED_Thompson_RED, EFJ_Thompson_EFJ, EFJ_Conchas_EFJ, RSAW_Thompson_RSAW, RSAW_Conchas_RSAW, RSA_Conchas_RSA, IND_BB_Conchas_IND_BB, SULF_Thompson_SULF)

dat_cond[is.nan(dat_cond)] <- NA

# log and scale transform response var
names(dat_cond)
dat_cond_log = dat_cond
dat_cond_log[,4:18] = log10(dat_cond_log[,4:18])
dat_cond_log[,4:18] = scale(dat_cond_log[,4:18])
sum(is.nan(dat_cond_log[,4:18]))
sum(is.na(dat_cond_log[,4:18]))
range(dat_cond_log[,4:18], na.rm = T)

### Plot response vars ###
par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
plot(dat_cond_log$mean_cond_uScm_AB00, type="o")
plot(dat_cond_log$mean_cond_uScm_AT07, type="o")
plot(dat_cond_log$mean_cond_uScm_GV01, type="o")
plot(dat_cond_log$mean_cond_uScm_HO00, type="o")
plot(dat_cond_log$mean_cond_uScm_MC06, type="o")
plot(dat_cond_log$mean_cond_uScm_RG01, type="o")
plot(dat_cond_log$mean_cond_uScm_RS02, type="o")
plot(dat_cond_log$mean_cond_uScm_SP02, type="o")
plot(dat_cond_log$mean_cond_uScm_EFJ, type="o")
plot(dat_cond_log$mean_cond_uScm_IND, type="o")
#plot(dat_cond_log$mean_cond_uScm_IND_AB, type="o")
plot(dat_cond_log$mean_cond_uScm_IND_BB, type="o")
plot(dat_cond_log$mean_cond_uScm_RED, type="o")
plot(dat_cond_log$mean_cond_uScm_RSA, type="o")
plot(dat_cond_log$mean_cond_uScm_RSAW, type="o")
plot(dat_cond_log$mean_cond_uScm_SULF, type="o")

### Plot response vars ###
par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
plot(dat_cond$mean_cond_uScm_AB00, type="o")
plot(dat_cond$mean_cond_uScm_AT07, type="o")
plot(dat_cond$mean_cond_uScm_GV01, type="o")
plot(dat_cond$mean_cond_uScm_HO00, type="o")
plot(dat_cond$mean_cond_uScm_MC06, type="o")
plot(dat_cond$mean_cond_uScm_RG01, type="o")
plot(dat_cond$mean_cond_uScm_RS02, type="o")
plot(dat_cond$mean_cond_uScm_SP02, type="o")
plot(dat_cond$mean_cond_uScm_EFJ, type="o")
plot(dat_cond$mean_cond_uScm_IND, type="o")
#plot(dat_cond$mean_cond_uScm_IND_AB, type="o")
plot(dat_cond$mean_cond_uScm_IND_BB, type="o")
plot(dat_cond$mean_cond_uScm_RED, type="o")
plot(dat_cond$mean_cond_uScm_RSA, type="o")
plot(dat_cond$mean_cond_uScm_RSAW, type="o")
plot(dat_cond$mean_cond_uScm_SULF, type="o")

#### Scenario 1 : all catchments are separate states - 12 catchments  #### 
# Pull out only response var
names(dat_cond_log)
# AB00,AT07,GV01,HO00,MC06,RG01,RS02,
#  EFJ, RED, RSA,RSAW,SULF
dat_dep <- t(dat_cond_log[,c(4,5,6,7,8,9,10,12,16,17,18,19)])
row.names(dat_dep)

# Make covariate inputs
# dat_cov <- dat_cond_log[,c(2:3, 20:56)]
# dat_cov <- t(scale(dat_cov))
# row.names(dat_cov)

# without short ts sites:
dat_cov <- dat_cond_log[,c(2:3, 
                           20:26, 28,32,33,34,35,
                           36:47, 50,51,49,54,52,53,56)]
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)

#### make C matrix
# CC <- matrix(list( # season 1
#                   "Season1", "Season1", "Season1", "Season1", 
#                   "Season1", "Season1", "Season1", "Season1", 
#                   "Season1", "Season1", "Season1", "Season1",
#                   "Season1", "Season1", "Season1", "Season1", 
#                   "Season1", "Season1", "Season1", "Season1", 
#                   "Season1", "Season1", "Season1", "Season1",
#                   "Season1", "Season1", "Season1", "Season1", 
#                   "Season1", "Season1", "Season1", "Season1",
#                   "Season1", "Season1", "Season1", "Season1",
#                   "Season1",
#                   # season 2
#                   "Season2", "Season2", "Season2", "Season2", 
#                   "Season2", "Season2", "Season2", "Season2",
#                   "Season2", "Season2", "Season2", "Season2",
#                   "Season2", "Season2", "Season2", "Season2", 
#                   "Season2", "Season2", "Season2", "Season2",
#                   "Season2", "Season2", "Season2", "Season2",
#                   "Season2", "Season2", "Season2", "Season2",
#                   "Season2", "Season2", "Season2", "Season2",
#                   "Season2", "Season2", "Season2", "Season2",
#                   "Season2",
#                   # precip by site
#                   "AB00_precip",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,"AT07_precip",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,"GV01_precip",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,"HO00_precip",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,"MC06_precip",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,"RG01_precip",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,"RS02_precip",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,"SP02_precip",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,"EFJ_precip",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,"IND_precip",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,"IND_AB_precip",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,"IND_BB_precip",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,"RED_precip",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,"RSA_precip",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,"RSAW_precip",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"SULF_precip",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   # fires by site -37
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"AB00_Tea",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"AB00_Jesusita",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"AT07_Jesusita",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"GV01_Gaviota",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"HO00_Gaviota",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"HO00_Sherpa",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"MC06_Tea",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"MC06_Jesusita",0,0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"RG01_Gaviota",0,0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"RG01_Sherpa",0,0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"RS02_Tea",0,0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"RS02_Jesusita",0,0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"SP02_Gap",0,0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"RED_Thompson",0,0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"EFJ_Thompson",0,0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"EFJ_Conchas",0,0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"RSAW_Thompson",0,0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"RSAW_Conchas",0,0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"RSA_Conchas",0,0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"IND_BB_Conchas",0,
#                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"SULF_Thompson"
#                   ),37,39)

# without short ts sites:
CC <- matrix(list( # season 1
  "Season1", "Season1", "Season1", "Season1", "Season1", "Season1", "Season1", "Season1", "Season1", "Season1", "Season1", "Season1",
  # season 2
  "Season2", "Season2", "Season2", "Season2", "Season2", "Season2", "Season2", "Season2","Season2", "Season2", "Season2", "Season2",
  # precip by site
  "AB00_precip",0,0,0,0,0,0,0,0,0,0,0,
  0,"AT07_precip",0,0,0,0,0,0,0,0,0,0,
  0,0,"GV01_precip",0,0,0,0,0,0,0,0,0,
  0,0,0,"HO00_precip",0,0,0,0,0,0,0,0,
  0,0,0,0,"MC06_precip",0,0,0,0,0,0,0,
  0,0,0,0,0,"RG01_precip",0,0,0,0,0,0,
  0,0,0,0,0,0,"RS02_precip",0,0,0,0,0,
  0,0,0,0,0,0,0,"EFJ_precip", 0,0,0,0,
  0,0,0,0,0,0,0,0,"RED_precip", 0,0,0,
  0,0,0,0,0,0,0,0,0,"RSA_precip", 0,0,
  0,0,0,0,0,0,0,0,0,0,"RSAW_precip",0,
  0,0,0,0,0,0,0,0,0,0,0,"SULF_precip",
  # fires by site -33
  "AB00_Tea",0,0,0,0,0,0,0,0,0,0,0,
  "AB00_Jesusita",0,0,0,0,0,0,0,0,0,0,0,
  0,"AT07_Jesusita",0,0,0,0,0,0,0,0,0,0,
  0,0,"GV01_Gaviota",0,0,0,0,0,0,0,0,0,
  0,0,0,"HO00_Gaviota",0,0,0,0,0,0,0,0,
  0,0,0,"HO00_Sherpa",0,0,0,0,0,0,0,0,
  0,0,0,0,"MC06_Tea",0,0,0,0,0,0,0,
  0,0,0,0,"MC06_Jesusita",0,0,0,0,0,0,0,
  0,0,0,0,0,"RG01_Gaviota",0,0,0,0,0,0,
  0,0,0,0,0,"RG01_Sherpa",0,0,0,0,0,0,
  0,0,0,0,0,0,"RS02_Tea",0,0,0,0,0,
  0,0,0,0,0,0,"RS02_Jesusita",0,0,0,0,0,
  0,0,0,0,0,0,0,"EFJ_Thompson",0,0,0,0,
  0,0,0,0,0,0,0,"EFJ_Conchas",0,0,0,0,
  0,0,0,0,0,0,0,0,"RED_Thompson",0,0,0,
  0,0,0,0,0,0,0,0,0,"RSA_Conchas",0,0,
  0,0,0,0,0,0,0,0,0,0,"RSAW_Thompson",0,
  0,0,0,0,0,0,0,0,0,0,"RSAW_Conchas",0,
  0,0,0,0,0,0,0,0,0,0,0,"SULF_Thompson"
),12,33)


# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observtion model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) 
fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# # fit EM by itself
# fit <- MARSS(y = dat_dep, model = mod_list,
#                 control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) 

# export model fit
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_1state_AB00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_2state_AB00_AT07_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_3state_AB00_AT07_GV01_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_4state_AB00_AT07_GV01_HO00_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_5state_AB00_AT07_GV01_HO00_MC06_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_6state_AB00_AT07_GV01_HO00_MC06_RG01_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_7state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_SP02_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_RED_EM.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_12state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_RED_RSA_RSAW_SULF_BFGS.rds")
#saveRDS(fit, file = "data_working/marss_test_run/fit_120521_12state_cond_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_RED_RSA_RSAW_SULF_mBFGS.rds")
saveRDS(fit, file = "data_working/marss_test_run/fit_120621_12state_cond_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_RED_RSA_RSAW_SULF_mBFGS.rds")


### DIAGNOSES ###
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## check for hidden errors
fit[["errors"]]
# lots of errors when IND sites were included, none when these were removed

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
#est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)

#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_1state_AB00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_2state_AB00_AT07_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_3state_AB00_AT07_GV01_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_4state_AB00_AT07_GV01_HO00_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_5state_AB00_AT07_GV01_HO00_MC06_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_6state_AB00_AT07_GV01_HO00_MC06_RG01_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_7state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_SP02_EM_hessian.rds") ****** NOT SAVED *******
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_8state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_RED_EM_hessian.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_12state_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_RED_RSA_RSAW_SULF_BFGS.rds")
#saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120521_12state_cond_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_RED_RSA_RSAW_SULF_mBFGS.rds")
saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_120621_12state_cond_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_RED_RSA_RSAW_SULF_mBFGS.rds")

CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)


### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00", "AT07", "GV01", "HO00", "MC06", "RG01", "RS02","EFJ","RED","RSA","RSAW","SULF")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>%
  #rename(Parameter = parm) %>%# rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "SP02_precip",
                                             "EFJ_precip", "RED_precip", "RSA_precip", "RSAW_precip",
                                             "SULF_precip",
                                             "AB00_Tea",
                                             "AB00_Jesusita",
                                             "AT07_Jesusita",
                                             "GV01_Gaviota",
                                             "HO00_Gaviota",
                                             "HO00_Sherpa",
                                             "MC06_Tea",
                                             "MC06_Jesusita",
                                             "RG01_Gaviota",
                                             "RG01_Sherpa",
                                             "RS02_Tea",
                                             "RS02_Jesusita",
                                             "EFJ_Thompson",
                                             "EFJ_Conchas",
                                             "RED_Thompson",
                                             "RSA_Conchas",
                                             "RSAW_Thompson",
                                             "RSAW_Conchas",
                                             "SULF_Thompson"
                                             )))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 12/5/2021") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))

CIs_fit_ed2 = CIs_fit_ed[!(CIs_fit_ed$Site=="RSA" & CIs_fit_ed$Parameter=="RSAW_precip"),] 
CIs_fit_ed2 = CIs_fit_ed2[!(CIs_fit_ed2$Site=="RSA" & CIs_fit_ed2$Parameter=="RSAW_Thompson"),] 
CIs_fit_ed2 = CIs_fit_ed2[!(CIs_fit_ed2$Site=="RSA" & CIs_fit_ed2$Parameter=="RSAW_Conchas"),] 
CIs_fit_ed2$region = c(rep("Coastal California",33),rep("Subalpine New Mexico",22))
#CIs_fit_ed2 = CIs_fit_ed2[CIs_fit_ed2$Site!="SP02",]

(RESULTS_ALL <-ggplot(CIs_fit_ed2, aes(Parameter, Est., color=region)) + 
  geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=.7) +
  geom_point(position=position_dodge(width=0.3), size=5) + 
  theme_bw()+
  theme(plot.title = element_text(size = 8)) +
  theme(axis.text = element_text(size = 8)) +
  geom_hline(aes(yintercept=0), linetype="dashed")+
  coord_flip() +
  labs(y = "",
       title = "Sp. Conductivity MARSS modeling results - 12/6/2021") +
  theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
  facet_wrap(vars(region, Site), scales = "free"))

ggsave("figures/MARSS_12states_spc_precip_fire_120621.pdf",RESULTS_ALL)

## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
# don't use residuals() - this will pull an entirely different df
resids <- MARSSresiduals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)
null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"
null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

#           dAIC df
# fit        0.0 50
# null.fit 257.4 36
# RESULT: covar model is better than null

### Plot response vars ###
par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
plot(dat_dep[1,], type="o")
plot(dat_dep[2,], type="o")
plot(dat_dep[3,], type="o")
plot(dat_dep[4,], type="o")
plot(dat_dep[5,], type="o")
plot(dat_dep[6,], type="o")
plot(dat_dep[7,], type="o")
plot(dat_dep[8,], type="o")
plot(dat_dep[8,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  #forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}

# These look ok.
# What you don't want is a consistent lag at 1 or 12.
# Patterns are bad (esp. sinusoidal), random is good.
# Should definitely examine these without seasonal effect to see how necessary this is.
# If there's a pattern here, 

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"),
         pch=16,
         xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  qqline(resids$model.residuals[i,])
  # qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16, 
  #        xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  # qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# state residuals - not looking great
# they are qq plots that should look like a straight line
# shapiro test scores should be closer to 1
# Press back arrow to see all 12 states
# ordered according to data - SB sites look ok
# flat lines likely due to low variation in data sites

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

#### Scenario 1 : all catchments are separate states - 11 catchments ####

# Removing SULF site per 02/04/22 meeting

# Pull out only response var
names(dat_cond_log)
# AB00, AT07, GV01, HO00, MC06, RG01, RS02,
# EFJ, RED, RSA, & RSAW
dat_dep <- t(dat_cond_log[,c(4:10,12,15:17)])
row.names(dat_dep)

# Make covariate inputs
# without short ts sites:
dat_cov <- dat_cond_log[,c(2:3, # seasonal covariates
                           19:25, 27, 30:32, # precip
                           34:45, 48, 49, 47, 52, 50, 51)] # fire (reordering to match c matrix below)
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)

#### make C matrix

# without short ts sites:
CC <- matrix(list( 
  # season 1
  "Season1", "Season1", "Season1", "Season1", 
  "Season1", "Season1", "Season1", "Season1", 
  "Season1", "Season1", "Season1",
  # season 2
  "Season2", "Season2", "Season2", "Season2", 
  "Season2", "Season2", "Season2", "Season2",
  "Season2", "Season2", "Season2",
  # precip by site
  "AB00_precip",0,0,0,0,0,0,0,0,0,0,
  0,"AT07_precip",0,0,0,0,0,0,0,0,0,
  0,0,"GV01_precip",0,0,0,0,0,0,0,0,
  0,0,0,"HO00_precip",0,0,0,0,0,0,0,
  0,0,0,0,"MC06_precip",0,0,0,0,0,0,
  0,0,0,0,0,"RG01_precip",0,0,0,0,0,
  0,0,0,0,0,0,"RS02_precip",0,0,0,0,
  0,0,0,0,0,0,0,"EFJ_precip", 0,0,0,
  0,0,0,0,0,0,0,0,"RED_precip", 0,0,
  0,0,0,0,0,0,0,0,0,"RSA_precip", 0,
  0,0,0,0,0,0,0,0,0,0,"RSAW_precip",
  # fires by site
  "AB00_Tea",0,0,0,0,0,0,0,0,0,0,
  "AB00_Jesusita",0,0,0,0,0,0,0,0,0,0,
  0,"AT07_Jesusita",0,0,0,0,0,0,0,0,0,
  0,0,"GV01_Gaviota",0,0,0,0,0,0,0,0,
  0,0,0,"HO00_Gaviota",0,0,0,0,0,0,0,
  0,0,0,"HO00_Sherpa",0,0,0,0,0,0,0,
  0,0,0,0,"MC06_Tea",0,0,0,0,0,0,
  0,0,0,0,"MC06_Jesusita",0,0,0,0,0,0,
  0,0,0,0,0,"RG01_Gaviota",0,0,0,0,0,
  0,0,0,0,0,"RG01_Sherpa",0,0,0,0,0,
  0,0,0,0,0,0,"RS02_Tea",0,0,0,0,
  0,0,0,0,0,0,"RS02_Jesusita",0,0,0,0,
  0,0,0,0,0,0,0,"EFJ_Thompson",0,0,0,
  0,0,0,0,0,0,0,"EFJ_Conchas",0,0,0,
  0,0,0,0,0,0,0,0,"RED_Thompson",0,0,
  0,0,0,0,0,0,0,0,0,"RSA_Conchas",0,
  0,0,0,0,0,0,0,0,0,0,"RSAW_Thompson",
  0,0,0,0,0,0,0,0,0,0,"RSAW_Conchas"),11,31)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) 

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# # fit EM by itself
# fit <- MARSS(y = dat_dep, model = mod_list,
#                 control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) 

# export model fit
saveRDS(fit, file = "data_working/marss_test_run/fit_062122_11state_cond_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_RED_RSA_RSAW_mBFGS.rds")

### DIAGNOSES ###

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
# better to do parametric/non-parametric bootstrapping once model is decided upon
# Maybe increase to over 100 boots, 100 is standard
# est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)

saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_062122_11state_cond_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_RED_RSA_RSAW_mBFGS.rds")

# formatting confidence intervals into dataframe
CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)

### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00","AT07","GV01","HO00","MC06","RG01","RS02","EFJ","RED","RSA","RSAW")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>%
  #rename(Parameter = parm) %>%# rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "SP02_precip",
                                             "EFJ_precip", "RED_precip", "RSA_precip", "RSAW_precip",
                                             "AB00_Tea",
                                             "AB00_Jesusita",
                                             "AT07_Jesusita",
                                             "GV01_Gaviota",
                                             "HO00_Gaviota",
                                             "HO00_Sherpa",
                                             "MC06_Tea",
                                             "MC06_Jesusita",
                                             "RG01_Gaviota",
                                             "RG01_Sherpa",
                                             "RS02_Tea",
                                             "RS02_Jesusita",
                                             "EFJ_Thompson",
                                             "EFJ_Conchas",
                                             "RED_Thompson",
                                             "RSA_Conchas",
                                             "RSAW_Thompson",
                                             "RSAW_Conchas")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 06/21/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))

# Adding labels for plotting purposes
CIs_fit_ed2 = CIs_fit_ed[!(CIs_fit_ed$Site=="RSA" & CIs_fit_ed$Parameter=="RSAW_precip"),] 
CIs_fit_ed2 = CIs_fit_ed2[!(CIs_fit_ed2$Site=="RSA" & CIs_fit_ed2$Parameter=="RSAW_Thompson"),] 
CIs_fit_ed2 = CIs_fit_ed2[!(CIs_fit_ed2$Site=="RSA" & CIs_fit_ed2$Parameter=="RSAW_Conchas"),] 
CIs_fit_ed2$region = c(rep("Coastal California",33),rep("Subalpine New Mexico",18))

(RESULTS_ALL <-ggplot(CIs_fit_ed2, aes(Parameter, Est., color=region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=.7) +
    geom_point(position=position_dodge(width=0.3), size=5) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 06/21/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(vars(region, Site), scales = "free"))

# ggsave("figures/MARSS_11states_cond_precip_fire_062122.png",
#        width = 40,
#        height = 20,
#        units = "cm")

## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- MARSSresiduals(fit)

# Save residuals dataframe
saveRDS(resids, "data_working/marss_test_run/all_resids_062122_11state_cond_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_RED_RSA_RSAW_mBFGS.rds")

kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
# No C matrix
# Need to be sure the covariates we've added are *actually* explaining
# the variation we are seeing
# Should have a better AIC score in model above

# Future to-do : make new model with only seasonal effects, and compare
# to null model to see how it performs. It may not be explaining much
# more variability than the null model and making the model above overly
# complex. Or could create model with covariates minus seasonal and compare.
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)

# fitting null model
# Note - if we change the structure of the model above, make sure that
# the same code is used to run the null models here below.
null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)
# dAIC- delta AIC
# 0.0 = always the value for the lowest model AIC
#           dAIC df
# fit        0.0 64
# null.fit 292.8 33
# RESULT: covar model is better than null

### Plot response vars ###
par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
plot(dat_dep[1,], type="o")
plot(dat_dep[2,], type="o")
plot(dat_dep[3,], type="o")
plot(dat_dep[4,], type="o")
plot(dat_dep[5,], type="o")
plot(dat_dep[6,], type="o")
plot(dat_dep[7,], type="o")
plot(dat_dep[8,], type="o")
plot(dat_dep[8,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  # forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}

# State residuals don't have any patterns that jump out.

# Error when trying to plot model residuals:
# Error in plot.window(...) : need finite 'ylim' values

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"),
         pch=16,
         xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  qqline(resids$model.residuals[i,])
  # qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16,
  #        xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  # qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# State residuals aren't great, but slightly better at SB.

# Error when trying to plot model residuals:
# Error in shapiro.test(resids$model.residuals[i, ]) : # all 'x' values are identical

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

#### Scenario 1 : all catchments are separate states - 11 catchments - w/o seasonal effect ####

# Removing seasonal effect to compare to 11 catchment model
# WITH seasonal effect above

# Pull out only response var
names(dat_cond_log)
# AB00, AT07, GV01, HO00, MC06, RG01, RS02,
# EFJ, RED, RSA, & RSAW
dat_dep <- t(dat_cond_log[,c(4:10,12,15:17)])
row.names(dat_dep)

# Make covariate inputs
# without short ts sites:
dat_cov <- dat_cond_log[,c(19:25, 27, 30:32, # precip
                           34:45, 48, 49, 47, 52, 50, 51)] # fire (reordering to match c matrix below)
dat_cov <- t(scale(dat_cov))
row.names(dat_cov)

#### make C matrix

# without short ts sites:
CC <- matrix(list( 
  # precip by site
  "AB00_precip",0,0,0,0,0,0,0,0,0,0,
  0,"AT07_precip",0,0,0,0,0,0,0,0,0,
  0,0,"GV01_precip",0,0,0,0,0,0,0,0,
  0,0,0,"HO00_precip",0,0,0,0,0,0,0,
  0,0,0,0,"MC06_precip",0,0,0,0,0,0,
  0,0,0,0,0,"RG01_precip",0,0,0,0,0,
  0,0,0,0,0,0,"RS02_precip",0,0,0,0,
  0,0,0,0,0,0,0,"EFJ_precip", 0,0,0,
  0,0,0,0,0,0,0,0,"RED_precip", 0,0,
  0,0,0,0,0,0,0,0,0,"RSA_precip", 0,
  0,0,0,0,0,0,0,0,0,0,"RSAW_precip",
  # fires by site
  "AB00_Tea",0,0,0,0,0,0,0,0,0,0,
  "AB00_Jesusita",0,0,0,0,0,0,0,0,0,0,
  0,"AT07_Jesusita",0,0,0,0,0,0,0,0,0,
  0,0,"GV01_Gaviota",0,0,0,0,0,0,0,0,
  0,0,0,"HO00_Gaviota",0,0,0,0,0,0,0,
  0,0,0,"HO00_Sherpa",0,0,0,0,0,0,0,
  0,0,0,0,"MC06_Tea",0,0,0,0,0,0,
  0,0,0,0,"MC06_Jesusita",0,0,0,0,0,0,
  0,0,0,0,0,"RG01_Gaviota",0,0,0,0,0,
  0,0,0,0,0,"RG01_Sherpa",0,0,0,0,0,
  0,0,0,0,0,0,"RS02_Tea",0,0,0,0,
  0,0,0,0,0,0,"RS02_Jesusita",0,0,0,0,
  0,0,0,0,0,0,0,"EFJ_Thompson",0,0,0,
  0,0,0,0,0,0,0,"EFJ_Conchas",0,0,0,
  0,0,0,0,0,0,0,0,"RED_Thompson",0,0,
  0,0,0,0,0,0,0,0,0,"RSA_Conchas",0,
  0,0,0,0,0,0,0,0,0,0,"RSAW_Thompson",
  0,0,0,0,0,0,0,0,0,0,"RSAW_Conchas"),11,29)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) 

fit_woseas <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# export model fit
#saveRDS(fit_woseas, file = "data_working/marss_test_run/fit_032222_11state_cond_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_RED_RSA_RSAW_mBFGS_woSeas.rds")

### DIAGNOSES ###

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit_woseas[["errors"]]
# NULL - Yay!

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit_woseas <- MARSSparamCIs(fit_woseas)
# better to do parametric/non-parametric bootstrapping once model is decided upon
# Maybe increase to over 100 boots, 100 is standard
# est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)

#saveRDS(est_fit_woseas, "data_working/marss_test_run/CIs_fit_032222_11state_cond_AB00_AT07_GV01_HO00_MC06_RG01_RS02_EFJ_RED_RSA_RSAW_mBFGS_woseas.rds")

# formatting confidence intervals into dataframe
CIs_fit = cbind(
  est_fit_woseas$par$U,
  est_fit_woseas$par.lowCI$U,
  est_fit_woseas$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)

### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00","AT07","GV01","HO00","MC06","RG01","RS02","EFJ","RED","RSA","RSAW")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>%
  #rename(Parameter = parm) %>%# rename site column
  mutate(Parameter = factor(parm, levels = c("AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "SP02_precip",
                                             "EFJ_precip", "RED_precip", "RSA_precip", "RSAW_precip",
                                             "AB00_Tea",
                                             "AB00_Jesusita",
                                             "AT07_Jesusita",
                                             "GV01_Gaviota",
                                             "HO00_Gaviota",
                                             "HO00_Sherpa",
                                             "MC06_Tea",
                                             "MC06_Jesusita",
                                             "RG01_Gaviota",
                                             "RG01_Sherpa",
                                             "RS02_Tea",
                                             "RS02_Jesusita",
                                             "EFJ_Thompson",
                                             "EFJ_Conchas",
                                             "RED_Thompson",
                                             "RSA_Conchas",
                                             "RSAW_Thompson",
                                             "RSAW_Conchas")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 03/22/2022 - Without Seasonal Term") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))

# Adding labels for plotting purposes
CIs_fit_ed2 = CIs_fit_ed[!(CIs_fit_ed$Site=="RSA" & CIs_fit_ed$Parameter=="RSAW_precip"),] 
CIs_fit_ed2 = CIs_fit_ed2[!(CIs_fit_ed2$Site=="RSA" & CIs_fit_ed2$Parameter=="RSAW_Thompson"),] 
CIs_fit_ed2 = CIs_fit_ed2[!(CIs_fit_ed2$Site=="RSA" & CIs_fit_ed2$Parameter=="RSAW_Conchas"),] 
CIs_fit_ed2$region = c(rep("Coastal California",19),rep("Subalpine New Mexico",10))

(RESULTS_ALL <-ggplot(CIs_fit_ed2, aes(Parameter, Est., color=region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=.7) +
    geom_point(position=position_dodge(width=0.3), size=5) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 02/24/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(vars(region, Site), scales = "free"))

# ggsave("figures/MARSS_11states_cond_woseas_032222.png",
#        width = 40,
#        height = 20,
#        units = "cm")

## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- MARSSresiduals(fit_woseas)
kf=print(fit_woseas, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
# No C matrix
# Need to be sure the covariates we've added are *actually* explaining
# the variation we are seeing
# Should have a better AIC score in model above

mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)

# fitting null model
# Note - if we change the structure of the model above, make sure that
# the same code is used to run the null models here below.
null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit_woseas, null.fit)
# dAIC- delta AIC
# 0.0 = always the value for the lowest model AIC
#                  dAIC df
# fit_woseas        0.0 62
# null.fit        296.3 33
# RESULT: covar model is better than null

# AND COMPARE TO MODEL ABOVE WITH SEASONS-
bbmle::AICtab(fit_woseas, fit)
#            dAIC df
# fit_woseas  0.0 62
# fit         3.5 64
# RESULT: covar model WITHOUT seasons is better than
# one WITH seasons - but only slightly?

### Plot response vars ###
par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
plot(dat_dep[1,], type="o")
plot(dat_dep[2,], type="o")
plot(dat_dep[3,], type="o")
plot(dat_dep[4,], type="o")
plot(dat_dep[5,], type="o")
plot(dat_dep[6,], type="o")
plot(dat_dep[7,], type="o")
plot(dat_dep[8,], type="o")
plot(dat_dep[8,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  #forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}

# Neither state nor model residuals have any patterns that jump out.

# Error when trying to plot model residuals:
# Error in plot.window(...) : need finite 'ylim' values

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"),
         pch=16,
         xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  qqline(resids$model.residuals[i,])
  # qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16,
  #        xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  # qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# State residuals aren't great, but slightly better at SB. Maybe look worse now that seasonal effect is removed?

# Error when trying to plot model residuals:
# Error in shapiro.test(resids$model.residuals[i, ]) : # all 'x' values are identical

# reset plotting window
# par(mfrow=c(1,1),oma = c(0, 0, 0, 0))
dev.off()

#### Scenario 1 : all catchments are separate states - new fire covariate ####

# New fire covariate - exponential decay over 4 years (per Verkaik et al. 2013)

dat_cond_d <- dat_new22 %>%
  select(site, index, Season1, Season2, 
         mean_cond_uScm, cumulative_precip_mm,
         AB00_Tea_d, AB00_Jesusita_d, AT07_Jesusita_d, GV01_Gaviota_d, HO00_Gaviota_d, HO00_Sherpa_d, MC06_Tea_d, MC06_Jesusita_d, RG01_Gaviota_d, RG01_Sherpa_d, RS02_Tea_d, RS02_Jesusita_d, SP02_Gap_d, RED_Thompson_d, EFJ_Thompson_d, EFJ_Conchas_d, RSAW_Thompson_d, RSAW_Conchas_d, RSA_Conchas_d, IND_BB_Conchas_d, SULF_Thompson_d) %>% 
  pivot_wider(names_from = site, values_from = c(mean_cond_uScm, cumulative_precip_mm, AB00_Tea_d, AB00_Jesusita_d, AT07_Jesusita_d, GV01_Gaviota_d, HO00_Gaviota_d, HO00_Sherpa_d, MC06_Tea_d, MC06_Jesusita_d, RG01_Gaviota_d, RG01_Sherpa_d, RS02_Tea_d, RS02_Jesusita_d, SP02_Gap_d, RED_Thompson_d, EFJ_Thompson_d, EFJ_Conchas_d, RSAW_Thompson_d, RSAW_Conchas_d, RSA_Conchas_d, IND_BB_Conchas_d, SULF_Thompson_d)) %>%
  select(index, Season1, Season2, 
         mean_cond_uScm_AB00, mean_cond_uScm_AT07, mean_cond_uScm_GV01, mean_cond_uScm_HO00, mean_cond_uScm_MC06, mean_cond_uScm_RG01, mean_cond_uScm_RS02, mean_cond_uScm_SP02, mean_cond_uScm_EFJ, mean_cond_uScm_IND, mean_cond_uScm_IND_BB, mean_cond_uScm_RED, mean_cond_uScm_RSA, mean_cond_uScm_RSAW, mean_cond_uScm_SULF,
         cumulative_precip_mm_AB00, cumulative_precip_mm_AT07, cumulative_precip_mm_GV01, cumulative_precip_mm_HO00, cumulative_precip_mm_MC06, cumulative_precip_mm_RG01, cumulative_precip_mm_RS02, cumulative_precip_mm_SP02, cumulative_precip_mm_EFJ, cumulative_precip_mm_IND, cumulative_precip_mm_IND_BB, cumulative_precip_mm_RED, cumulative_precip_mm_RSA, cumulative_precip_mm_RSAW, cumulative_precip_mm_SULF,
         AB00_Tea_d_AB00, AB00_Jesusita_d_AB00, AT07_Jesusita_d_AT07, GV01_Gaviota_d_GV01, HO00_Gaviota_d_HO00, HO00_Sherpa_d_HO00, MC06_Tea_d_MC06, MC06_Jesusita_d_MC06, RG01_Gaviota_d_RG01, RG01_Sherpa_d_RG01, RS02_Tea_d_RS02, RS02_Jesusita_d_RS02, SP02_Gap_d_SP02, RED_Thompson_d_RED, EFJ_Thompson_d_EFJ, EFJ_Conchas_d_EFJ, RSAW_Thompson_d_RSAW, RSAW_Conchas_d_RSAW, RSA_Conchas_d_RSA, IND_BB_Conchas_d_IND_BB, SULF_Thompson_d_SULF)

dat_cond_d[is.nan(dat_cond_d)] <- NA

# log and scale transform response var
names(dat_cond_d)
dat_cond_d_log = dat_cond_d
dat_cond_d_log[,4:18] = log10(dat_cond_d_log[,4:18])
dat_cond_d_log[,4:18] = scale(dat_cond_d_log[,4:18])
sum(is.nan(dat_cond_d_log[,4:18]))
sum(is.na(dat_cond_d_log[,4:18]))
range(dat_cond_d_log[,4:18], na.rm = T)

# Pull out only response var
names(dat_cond_d_log)
# AB00, AT07, GV01, HO00, MC06, RG01, RS02,
# EFJ, RED, RSA, & RSAW
dat_dep_d <- t(dat_cond_d_log[,c(4:10,12,15:17)])
row.names(dat_dep_d)

# Make covariate inputs
# without short ts sites:
dat_cov_d <- dat_cond_d_log[,c(2:3, # seasonal covariates
                           19:25, 27, 30:32, # precip
                           34:42, 44, 45, 48, 49, 47, 52, 50, 51)] # fire (reordering to match c matrix below)

dat_cov_d <- t(scale(dat_cov_d))
row.names(dat_cov_d)
sum(is.nan(dat_cov_d)) # RG01_Sherpa causing problems so removed it from above
sum(is.na(dat_cov_d))

#### make C matrix

# without short ts sites:
CC <- matrix(list( 
  # season 1
  "Season1", "Season1", "Season1", "Season1", 
  "Season1", "Season1", "Season1", "Season1", 
  "Season1", "Season1", "Season1",
  # season 2
  "Season2", "Season2", "Season2", "Season2", 
  "Season2", "Season2", "Season2", "Season2",
  "Season2", "Season2", "Season2",
  # precip by site
  "AB00_precip",0,0,0,0,0,0,0,0,0,0,
  0,"AT07_precip",0,0,0,0,0,0,0,0,0,
  0,0,"GV01_precip",0,0,0,0,0,0,0,0,
  0,0,0,"HO00_precip",0,0,0,0,0,0,0,
  0,0,0,0,"MC06_precip",0,0,0,0,0,0,
  0,0,0,0,0,"RG01_precip",0,0,0,0,0,
  0,0,0,0,0,0,"RS02_precip",0,0,0,0,
  0,0,0,0,0,0,0,"EFJ_precip", 0,0,0,
  0,0,0,0,0,0,0,0,"RED_precip", 0,0,
  0,0,0,0,0,0,0,0,0,"RSA_precip", 0,
  0,0,0,0,0,0,0,0,0,0,"RSAW_precip",
  # fires by site
  "AB00_Tea",0,0,0,0,0,0,0,0,0,0,
  "AB00_Jesusita",0,0,0,0,0,0,0,0,0,0,
  0,"AT07_Jesusita",0,0,0,0,0,0,0,0,0,
  0,0,"GV01_Gaviota",0,0,0,0,0,0,0,0,
  0,0,0,"HO00_Gaviota",0,0,0,0,0,0,0,
  0,0,0,"HO00_Sherpa",0,0,0,0,0,0,0,
  0,0,0,0,"MC06_Tea",0,0,0,0,0,0,
  0,0,0,0,"MC06_Jesusita",0,0,0,0,0,0,
  0,0,0,0,0,"RG01_Gaviota",0,0,0,0,0,
  #0,0,0,0,0,"RG01_Sherpa",0,0,0,0,0,
  0,0,0,0,0,0,"RS02_Tea",0,0,0,0,
  0,0,0,0,0,0,"RS02_Jesusita",0,0,0,0,
  0,0,0,0,0,0,0,"EFJ_Thompson",0,0,0,
  0,0,0,0,0,0,0,"EFJ_Conchas",0,0,0,
  0,0,0,0,0,0,0,0,"RED_Thompson",0,0,
  0,0,0,0,0,0,0,0,0,"RSA_Conchas",0,
  0,0,0,0,0,0,0,0,0,0,"RSAW_Thompson",
  0,0,0,0,0,0,0,0,0,0,"RSAW_Conchas"),11,30)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov_d,
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep_d, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) 

fit <- MARSS(y = dat_dep_d, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# # fit EM by itself
# fit <- MARSS(y = dat_dep, model = mod_list,
#                 control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) 

# export model fit
saveRDS(fit, file = "data_working/marss_test_run/fit_062122_11state_cond_decay_mBFGS.rds")

### DIAGNOSES ###

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
# better to do parametric/non-parametric bootstrapping once model is decided upon
# Maybe increase to over 100 boots, 100 is standard
# est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)

saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_062122_11state_cond_decay_mBFGS.rds")

# formatting confidence intervals into dataframe
CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)

### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00","AT07","GV01","HO00","MC06","RG01","RS02","EFJ","RED","RSA","RSAW")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  dplyr::rename(Site = i) %>%
  #rename(Parameter = parm) %>%# rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "SP02_precip",
                                             "EFJ_precip", "RED_precip", "RSA_precip", "RSAW_precip",
                                             "AB00_Tea",
                                             "AB00_Jesusita",
                                             "AT07_Jesusita",
                                             "GV01_Gaviota",
                                             "HO00_Gaviota",
                                             "HO00_Sherpa",
                                             "MC06_Tea",
                                             "MC06_Jesusita",
                                             "RG01_Gaviota",
                                             "RS02_Tea",
                                             "RS02_Jesusita",
                                             "EFJ_Thompson",
                                             "EFJ_Conchas",
                                             "RED_Thompson",
                                             "RSA_Conchas",
                                             "RSAW_Thompson",
                                             "RSAW_Conchas")))

# plot results
(RESULTS_ALL_d <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 06/21/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))

# Adding labels for plotting purposes
CIs_fit_ed2 = CIs_fit_ed[!(CIs_fit_ed$Site=="RSA" & CIs_fit_ed$Parameter=="RSAW_precip"),] 
CIs_fit_ed2 = CIs_fit_ed2[!(CIs_fit_ed2$Site=="RSA" & CIs_fit_ed2$Parameter=="RSAW_Thompson"),] 
CIs_fit_ed2 = CIs_fit_ed2[!(CIs_fit_ed2$Site=="RSA" & CIs_fit_ed2$Parameter=="RSAW_Conchas"),] 
CIs_fit_ed2$region = c(rep("Coastal California",32),rep("Subalpine New Mexico",18))

(RESULTS_ALL_d <-ggplot(CIs_fit_ed2, aes(Parameter, Est., color=region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=.7) +
    geom_point(position=position_dodge(width=0.3), size=5) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 06/21/2022\n11 state - 4-year fire decay term") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(vars(region, Site), scales = "free"))

# ggsave("figures/MARSS_11states_cond_precip_fire_decay_062122.png",
#        width = 40,
#        height = 20,
#        units = "cm")

## Script for diagnoses ###

dat = dat_dep_d
time = c(1:ncol(dat_dep_d))
resids <- MARSSresiduals(fit)

# Save resids dataframe
saveRDS(resids, file = "data_working/marss_test_run/all_resids_062122_11state_cond_decay_mBFGS.rds")

kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
# No C matrix
# Need to be sure the covariates we've added are *actually* explaining
# the variation we are seeing
# Should have a better AIC score in model above

mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)

# fitting null model
# Note - if we change the structure of the model above, make sure that
# the same code is used to run the null models here below.
null.kemfit <- MARSS(y = dat_dep_d, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep_d, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)
# dAIC- delta AIC
# 0.0 = always the value for the lowest model AIC
#           dAIC df
# fit        0.0 63
# null.fit 350.7 33
# RESULT: covar model is better than null

### Plot response vars ###
par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
plot(dat_dep[1,], type="o")
plot(dat_dep[2,], type="o")
plot(dat_dep[3,], type="o")
plot(dat_dep[4,], type="o")
plot(dat_dep[5,], type="o")
plot(dat_dep[6,], type="o")
plot(dat_dep[7,], type="o")
plot(dat_dep[8,], type="o")
plot(dat_dep[8,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  #forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}

# State residuals don't have any patterns that jump out.

# Error when trying to plot model residuals:
# Error in plot.window(...) : need finite 'ylim' values

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"),
         pch=16,
         xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  qqline(resids$model.residuals[i,])
  # qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16,
  #        xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  # qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# State residuals aren't great, but slightly better at SB.

# Error when trying to plot model residuals:
# Error in shapiro.test(resids$model.residuals[i, ]) : # all 'x' values are identical

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

#### Scenario 1 : all catchments are separate states - NO fire covariate ####

dat_cond_nof <- dat_new22 %>%
  select(site, index, Season1, Season2, 
         mean_cond_uScm, cumulative_precip_mm) %>% 
  pivot_wider(names_from = site, values_from = c(mean_cond_uScm, cumulative_precip_mm)) %>%
  select(index, Season1, Season2, 
         mean_cond_uScm_AB00, mean_cond_uScm_AT07, mean_cond_uScm_GV01, mean_cond_uScm_HO00, mean_cond_uScm_MC06, mean_cond_uScm_RG01, mean_cond_uScm_RS02, mean_cond_uScm_SP02, mean_cond_uScm_EFJ, mean_cond_uScm_IND, mean_cond_uScm_IND_BB, mean_cond_uScm_RED, mean_cond_uScm_RSA, mean_cond_uScm_RSAW, mean_cond_uScm_SULF,
         cumulative_precip_mm_AB00, cumulative_precip_mm_AT07, cumulative_precip_mm_GV01, cumulative_precip_mm_HO00, cumulative_precip_mm_MC06, cumulative_precip_mm_RG01, cumulative_precip_mm_RS02, cumulative_precip_mm_SP02, cumulative_precip_mm_EFJ, cumulative_precip_mm_IND, cumulative_precip_mm_IND_BB, cumulative_precip_mm_RED, cumulative_precip_mm_RSA, cumulative_precip_mm_RSAW, cumulative_precip_mm_SULF)

dat_cond_nof[is.nan(dat_cond_nof)] <- NA

# log and scale transform response var
names(dat_cond_nof)
dat_cond_nof_log = dat_cond_nof
dat_cond_nof_log[,4:18] = log10(dat_cond_nof_log[,4:18])
dat_cond_nof_log[,4:18] = scale(dat_cond_nof_log[,4:18])
sum(is.nan(dat_cond_nof_log[,4:18]))
sum(is.na(dat_cond_nof_log[,4:18]))
range(dat_cond_nof_log[,4:18], na.rm = T)

# Pull out only response var
names(dat_cond_nof_log)
# AB00, AT07, GV01, HO00, MC06, RG01, RS02,
# EFJ, RED, RSA, & RSAW
dat_dep_nof <- t(dat_cond_nof_log[,c(4:10,12,15:17)])
row.names(dat_dep_nof)

# Make covariate inputs
# without short ts sites:
dat_cov_nof <- dat_cond_nof_log[,c(2:3, # seasonal covariates
                               19:25, 27, 30:32)] # precip

dat_cov_nof <- t(scale(dat_cov_nof))
row.names(dat_cov_nof)
sum(is.nan(dat_cov_nof)) # RG01_Sherpa causing problems so removed it from above
sum(is.na(dat_cov_nof))

#### make C matrix

# without short ts sites:
CC <- matrix(list( 
  # season 1
  "Season1", "Season1", "Season1", "Season1", 
  "Season1", "Season1", "Season1", "Season1", 
  "Season1", "Season1", "Season1",
  # season 2
  "Season2", "Season2", "Season2", "Season2", 
  "Season2", "Season2", "Season2", "Season2",
  "Season2", "Season2", "Season2",
  # precip by site
  "AB00_precip",0,0,0,0,0,0,0,0,0,0,
  0,"AT07_precip",0,0,0,0,0,0,0,0,0,
  0,0,"GV01_precip",0,0,0,0,0,0,0,0,
  0,0,0,"HO00_precip",0,0,0,0,0,0,0,
  0,0,0,0,"MC06_precip",0,0,0,0,0,0,
  0,0,0,0,0,"RG01_precip",0,0,0,0,0,
  0,0,0,0,0,0,"RS02_precip",0,0,0,0,
  0,0,0,0,0,0,0,"EFJ_precip", 0,0,0,
  0,0,0,0,0,0,0,0,"RED_precip", 0,0,
  0,0,0,0,0,0,0,0,0,"RSA_precip", 0,
  0,0,0,0,0,0,0,0,0,0,"RSAW_precip"),11,13)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov_nof,
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep_nof, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) 

fit <- MARSS(y = dat_dep_nof, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# # fit EM by itself
# fit <- MARSS(y = dat_dep, model = mod_list,
#                 control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) 

# export model fit
saveRDS(fit, file = "data_working/marss_test_run/fit_062122_11state_cond_nofire_mBFGS.rds")

### DIAGNOSES ###

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
# better to do parametric/non-parametric bootstrapping once model is decided upon
# Maybe increase to over 100 boots, 100 is standard
# est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)

saveRDS(est_fit, "data_working/marss_test_run/CIs_fit_062122_11state_cond_nofire_mBFGS.rds")

# formatting confidence intervals into dataframe
CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)

### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00","AT07","GV01","HO00","MC06","RG01","RS02","EFJ","RED","RSA","RSAW")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[1:2,], CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  dplyr::rename(Site = i) %>%
  #rename(Parameter = parm) %>%# rename site column
  mutate(Parameter = factor(parm, levels = c("Season1", "Season2", # relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "SP02_precip",
                                             "EFJ_precip", "RED_precip", "RSA_precip", "RSAW_precip")))

# plot results
(RESULTS_ALL_d <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 06/21/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Site, scales = "free"))

# Adding labels for plotting purposes
CIs_fit_ed2 = CIs_fit_ed[!(CIs_fit_ed$Site=="RSA" & CIs_fit_ed$Parameter=="RSAW_precip"),] 
CIs_fit_ed2$region = c(rep("Coastal California",21),rep("Subalpine New Mexico",12))

(RESULTS_ALL_d <-ggplot(CIs_fit_ed2, aes(Parameter, Est., color=region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=.7) +
    geom_point(position=position_dodge(width=0.3), size=5) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 06/21/2022\n11 state - no fire term") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(vars(region, Site), scales = "free"))

# ggsave("figures/MARSS_11states_cond_precip_nofire_062122.png",
#        width = 40,
#        height = 20,
#        units = "cm")

## Script for diagnoses ###

dat = dat_dep_nof
time = c(1:ncol(dat_dep_nof))
resids <- MARSSresiduals(fit)

# Save resids dataframe
saveRDS(resids, file = "data_working/marss_test_run/all_resids_062122_11state_cond_nofire_mBFGS.rds")

kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
# No C matrix
# Need to be sure the covariates we've added are *actually* explaining
# the variation we are seeing
# Should have a better AIC score in model above
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)

# fitting null model
# Note - if we change the structure of the model above, make sure that
# the same code is used to run the null models here below.
null.kemfit <- MARSS(y = dat_dep_nof, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep_nof, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)
# dAIC- delta AIC
# 0.0 = always the value for the lowest model AIC
#           dAIC df
# fit        0.0 46
# null.fit   257 33
# RESULT: covar model is better than null

### Plot response vars ###
par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
plot(dat_dep[1,], type="o")
plot(dat_dep[2,], type="o")
plot(dat_dep[3,], type="o")
plot(dat_dep[4,], type="o")
plot(dat_dep[5,], type="o")
plot(dat_dep[6,], type="o")
plot(dat_dep[7,], type="o")
plot(dat_dep[8,], type="o")
plot(dat_dep[8,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  #forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}

# State residuals don't have any patterns that jump out.

# Error when trying to plot model residuals:
# Error in plot.window(...) : need finite 'ylim' values

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  # qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"),
  #        pch=16,
  #        xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  # qqline(resids$model.residuals[i,])
  qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16,
         xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# State residuals aren't great, but slightly better at SB.

# Error when trying to plot model residuals:
# Error in shapiro.test(resids$model.residuals[i, ]) : # all 'x' values are identical

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

#### Scenario 1 : all catchments are separate states - 11 catchments w/ interaction and binary fire term ####

# no seasonal variables
# fire variable = 0/1 binary and instantaneous for 1 year
# interaction term = precip * fire

# Pull out only response var
names(dat_cond_log)
# AB00, AT07, GV01, HO00, MC06, RG01, RS02 = SB
# EFJ, RED, RSA, & RSAW = VC
dat_dep <- t(dat_cond_log[,c(4:10,12,15:17)])
row.names(dat_dep)

# Make covariate inputs
# without short ts sites:
dat_cov <- dat_cond_log[,c(19:25, 27, 30:32, # precip
                           34:45, 48, 49, 47, 52, 50, 51)] # fire (reordering to match c matrix below)

# Make new interaction term
dat_cov <- dat_cov %>%
  mutate(interact_AB00 = cumulative_precip_mm_AB00 * (AB00_Tea_AB00 + AB00_Jesusita_AB00),
         interact_AT07 = cumulative_precip_mm_AT07 * AT07_Jesusita_AT07,
         interact_GV01 = cumulative_precip_mm_GV01 * GV01_Gaviota_GV01,
         interact_HO00 = cumulative_precip_mm_HO00 * (HO00_Gaviota_HO00 + HO00_Sherpa_HO00),
         interact_MC06 = cumulative_precip_mm_MC06 * (MC06_Tea_MC06 + MC06_Jesusita_MC06),
         interact_RG01 = cumulative_precip_mm_RG01 * (RG01_Gaviota_RG01 + RG01_Sherpa_RG01),
         interact_RS02 = cumulative_precip_mm_RS02 * (RS02_Tea_RS02 + RS02_Jesusita_RS02),
         interact_EFJ = cumulative_precip_mm_EFJ * (EFJ_Thompson_EFJ + EFJ_Conchas_EFJ),
         interact_RED = cumulative_precip_mm_RED * RED_Thompson_RED,
         interact_RSA = cumulative_precip_mm_RSA * RSA_Conchas_RSA,
         interact_RSAW = cumulative_precip_mm_RSAW * (RSAW_Thompson_RSAW + RSAW_Conchas_RSAW))

dat_cov <- t(scale(dat_cov))
row.names(dat_cov)

#### make C matrix

# without short ts sites:
CC <- matrix(list( 
  # precip by site
  "AB00_precip",0,0,0,0,0,0,0,0,0,0,
  0,"AT07_precip",0,0,0,0,0,0,0,0,0,
  0,0,"GV01_precip",0,0,0,0,0,0,0,0,
  0,0,0,"HO00_precip",0,0,0,0,0,0,0,
  0,0,0,0,"MC06_precip",0,0,0,0,0,0,
  0,0,0,0,0,"RG01_precip",0,0,0,0,0,
  0,0,0,0,0,0,"RS02_precip",0,0,0,0,
  0,0,0,0,0,0,0,"EFJ_precip", 0,0,0,
  0,0,0,0,0,0,0,0,"RED_precip", 0,0,
  0,0,0,0,0,0,0,0,0,"RSA_precip", 0,
  0,0,0,0,0,0,0,0,0,0,"RSAW_precip",
  # fires by site
  "AB00_Tea",0,0,0,0,0,0,0,0,0,0,
  "AB00_Jesusita",0,0,0,0,0,0,0,0,0,0,
  0,"AT07_Jesusita",0,0,0,0,0,0,0,0,0,
  0,0,"GV01_Gaviota",0,0,0,0,0,0,0,0,
  0,0,0,"HO00_Gaviota",0,0,0,0,0,0,0,
  0,0,0,"HO00_Sherpa",0,0,0,0,0,0,0,
  0,0,0,0,"MC06_Tea",0,0,0,0,0,0,
  0,0,0,0,"MC06_Jesusita",0,0,0,0,0,0,
  0,0,0,0,0,"RG01_Gaviota",0,0,0,0,0,
  0,0,0,0,0,"RG01_Sherpa",0,0,0,0,0,
  0,0,0,0,0,0,"RS02_Tea",0,0,0,0,
  0,0,0,0,0,0,"RS02_Jesusita",0,0,0,0,
  0,0,0,0,0,0,0,"EFJ_Thompson",0,0,0,
  0,0,0,0,0,0,0,"EFJ_Conchas",0,0,0,
  0,0,0,0,0,0,0,0,"RED_Thompson",0,0,
  0,0,0,0,0,0,0,0,0,"RSA_Conchas",0,
  0,0,0,0,0,0,0,0,0,0,"RSAW_Thompson",
  0,0,0,0,0,0,0,0,0,0,"RSAW_Conchas",
  # interaction terms by site
  "AB00_interact",0,0,0,0,0,0,0,0,0,0,
  0,"AT07_interact",0,0,0,0,0,0,0,0,0,
  0,0,"GV01_interact",0,0,0,0,0,0,0,0,
  0,0,0,"HO00_interact",0,0,0,0,0,0,0,
  0,0,0,0,"MC06_interact",0,0,0,0,0,0,
  0,0,0,0,0,"RG01_interact",0,0,0,0,0,
  0,0,0,0,0,0,"RS02_interact",0,0,0,0,
  0,0,0,0,0,0,0,"EFJ_interact", 0,0,0,
  0,0,0,0,0,0,0,0,"RED_interact", 0,0,
  0,0,0,0,0,0,0,0,0,"RSA_interact", 0,
  0,0,0,0,0,0,0,0,0,0,"RSAW_interact"),11,40)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov,
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) 

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# # fit EM by itself
# fit <- MARSS(y = dat_dep, model = mod_list,
#                 control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) 

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_071122_11state_cond_winteraction_mBFGS.rds")

# import model fit
fit = readRDS("data_working/marss_test_run/fit_071122_11state_cond_winteraction_mBFGS.rds")

### DIAGNOSES ###

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
# better to do parametric/non-parametric bootstrapping once model is decided upon
# Maybe increase to over 100 boots, 100 is standard
# est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)

saveRDS(est_fit, 
        "data_working/marss_test_run/CIs_fit_071122_11state_cond_winteraction_mBFGS.rds")

# formatting confidence intervals into dataframe
CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)

### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00","AT07","GV01","HO00","MC06","RG01","RS02","EFJ","RED","RSA","RSAW")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  rename(Site = i) %>%
  #rename(Parameter = parm) %>%# rename site column
  mutate(Parameter = factor(parm, levels = c("AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "EFJ_precip",
                                             "RED_precip", "RSA_precip", "RSAW_precip",
                                             "AB00_Tea", "AB00_Jesusita", "AT07_Jesusita",
                                             "GV01_Gaviota", "HO00_Gaviota", "HO00_Sherpa",
                                             "MC06_Tea", "MC06_Jesusita", "RG01_Gaviota",
                                             "RG01_Sherpa", "RS02_Tea", "RS02_Jesusita",
                                             "EFJ_Thompson", "EFJ_Conchas", "RED_Thompson",
                                             "RSA_Conchas", "RSAW_Thompson", "RSAW_Conchas",
                                             "AB00_interact", "AT07_interact", "GV01_interact",
                                             "HO00_interact", "MC06_interact", "RG01_interact",
                                             "RS02_interact", "SP02_interact", "EFJ_interact",
                                             "RED_interact", "RSA_interact", "RSAW_interact")))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 07/11/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(.~Site, scales = "free"))

# Adding labels for plotting purposes
CIs_fit_ed2 = CIs_fit_ed[!(CIs_fit_ed$Site=="RSA" & CIs_fit_ed$Parameter=="RSAW_precip"),] 
CIs_fit_ed2 = CIs_fit_ed2[!(CIs_fit_ed2$Site=="RSA" & CIs_fit_ed2$Parameter=="RSAW_Thompson"),] 
CIs_fit_ed2 = CIs_fit_ed2[!(CIs_fit_ed2$Site=="RSA" & CIs_fit_ed2$Parameter=="RSAW_Conchas"),] 
CIs_fit_ed2 = CIs_fit_ed2[!(CIs_fit_ed2$Site=="RSA" & CIs_fit_ed2$Parameter=="RSAW_interact"),] 
CIs_fit_ed2$region = c(rep("Coastal California",26),rep("Subalpine New Mexico",14))

(RESULTS_ALL <-ggplot(CIs_fit_ed2, aes(Parameter, Est., color=region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=.7) +
    geom_point(position=position_dodge(width=0.3), size=5) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 07/11/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(vars(region, Site), scales = "free"))

# ggsave("figures/MARSS_11states_cond_winteraction_071122.png",
#        width = 40,
#        height = 20,
#        units = "cm")

## Script for diagnoses ###

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- MARSSresiduals(fit)

# Save residuals dataframe
saveRDS(resids, 
        "data_working/marss_test_run/all_resids_071122_11state_cond_winteraction_mBFGS.rds")

kf = print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
# No C matrix
# Need to be sure the covariates we've added are *actually* explaining
# the variation we are seeing
# Should have a better AIC score in model above

mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)

# fitting null model
# Note - if we change the structure of the model above, make sure that
# the same code is used to run the null models here below.
null.kemfit <- MARSS(y = dat_dep, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)

# dAIC- delta AIC
# 0.0 = always the value for the lowest model AIC
#           dAIC df
# fit        0.0 73
# null.fit 325.3 33
# RESULT: covar model is better than null

### Plot response vars ###
# par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
# plot(dat_dep[1,], type="o")
# plot(dat_dep[2,], type="o")
# plot(dat_dep[3,], type="o")
# plot(dat_dep[4,], type="o")
# plot(dat_dep[5,], type="o")
# plot(dat_dep[6,], type="o")
# plot(dat_dep[7,], type="o")
# plot(dat_dep[8,], type="o")
# plot(dat_dep[8,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  # forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}

# State residuals don't have any patterns that jump out.

# Error when trying to plot model residuals:
# Error in plot.window(...) : need finite 'ylim' values

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  # qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"),
  #        pch=16,
  #        xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  # qqline(resids$model.residuals[i,])
  qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16,
         xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# State residuals again aren't great, but slightly better at SB.

# Error when trying to plot model residuals:
# Error in shapiro.test(resids$model.residuals[i, ]) : # all 'x' values are identical

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

### residuals analysis ###
# NOTES:
# Sections 12.3.1 and 12.3.2 in the 2021 edition of the MARSS manual describe two complementary ways (outlier detection & level change detection from standardized smoothed model & state residuals, respectively) to detect change in a time series using residuals. See this example also explained in the lab manual in section 6.5.1, here https://atsa-es.github.io/atsa-labs/sec-uss-basic-diagnostics.html. From looking at the code made available with the lab manual (https://atsa-es.github.io/atsa-labs/chap-univariate-state-space.html#data-and-packages-3), I think you would get these two residual types from a model fit using the code below. 

resids <- residuals(fit, type='tT', standardization="marginal")

# for outlier detection
mresids <- subset(resids, name == "model") 
plot(mresids$t, mresids$.std.resids,
     ylab = "model std smoothationl", xlab = "", main = "data outliers")
abline(h = 0)
abline(h = c(2,-2), lty=2)

# for sudden level changes detection
sresids <- subset(resids, name == "state")
sresids$.rownames = as.factor(sresids$.rownames)
plot(sresids$t, sresids$.std.resids, type="l", col=(sresids$.rownames),
     ylab = "state std smoothation", xlab = "", main = "sudden level changes")
abline(h = 0)
abline(h = c(2,-2), lty=2)

#### Scenario 1 : all catchments are separate states - 11 catchments w/ interaction and decay fire covariate ####

# no seasonal variables
# fire variable = exponential decay over 4 years
# interaction term = precip * fire

dat_cond_d <- dat_new22 %>%
  select(site, index, Season1, Season2, 
         mean_cond_uScm, cumulative_precip_mm,
         AB00_Tea_d, AB00_Jesusita_d, AT07_Jesusita_d, GV01_Gaviota_d, HO00_Gaviota_d, HO00_Sherpa_d, MC06_Tea_d, MC06_Jesusita_d, RG01_Gaviota_d, RG01_Sherpa_d, RS02_Tea_d, RS02_Jesusita_d, SP02_Gap_d, RED_Thompson_d, EFJ_Thompson_d, EFJ_Conchas_d, RSAW_Thompson_d, RSAW_Conchas_d, RSA_Conchas_d, IND_BB_Conchas_d, SULF_Thompson_d) %>% 
  pivot_wider(names_from = site, values_from = c(mean_cond_uScm, cumulative_precip_mm, AB00_Tea_d, AB00_Jesusita_d, AT07_Jesusita_d, GV01_Gaviota_d, HO00_Gaviota_d, HO00_Sherpa_d, MC06_Tea_d, MC06_Jesusita_d, RG01_Gaviota_d, RG01_Sherpa_d, RS02_Tea_d, RS02_Jesusita_d, SP02_Gap_d, RED_Thompson_d, EFJ_Thompson_d, EFJ_Conchas_d, RSAW_Thompson_d, RSAW_Conchas_d, RSA_Conchas_d, IND_BB_Conchas_d, SULF_Thompson_d)) %>%
  select(index, Season1, Season2, 
         mean_cond_uScm_AB00, mean_cond_uScm_AT07, mean_cond_uScm_GV01, mean_cond_uScm_HO00, mean_cond_uScm_MC06, mean_cond_uScm_RG01, mean_cond_uScm_RS02, mean_cond_uScm_SP02, mean_cond_uScm_EFJ, mean_cond_uScm_IND, mean_cond_uScm_IND_BB, mean_cond_uScm_RED, mean_cond_uScm_RSA, mean_cond_uScm_RSAW, mean_cond_uScm_SULF,
         cumulative_precip_mm_AB00, cumulative_precip_mm_AT07, cumulative_precip_mm_GV01, cumulative_precip_mm_HO00, cumulative_precip_mm_MC06, cumulative_precip_mm_RG01, cumulative_precip_mm_RS02, cumulative_precip_mm_SP02, cumulative_precip_mm_EFJ, cumulative_precip_mm_IND, cumulative_precip_mm_IND_BB, cumulative_precip_mm_RED, cumulative_precip_mm_RSA, cumulative_precip_mm_RSAW, cumulative_precip_mm_SULF,
         AB00_Tea_d_AB00, AB00_Jesusita_d_AB00, AT07_Jesusita_d_AT07, GV01_Gaviota_d_GV01, HO00_Gaviota_d_HO00, HO00_Sherpa_d_HO00, MC06_Tea_d_MC06, MC06_Jesusita_d_MC06, RG01_Gaviota_d_RG01, RG01_Sherpa_d_RG01, RS02_Tea_d_RS02, RS02_Jesusita_d_RS02, SP02_Gap_d_SP02, RED_Thompson_d_RED, EFJ_Thompson_d_EFJ, EFJ_Conchas_d_EFJ, RSAW_Thompson_d_RSAW, RSAW_Conchas_d_RSAW, RSA_Conchas_d_RSA, IND_BB_Conchas_d_IND_BB, SULF_Thompson_d_SULF)

dat_cond_d[is.nan(dat_cond_d)] <- NA

# log and scale transform response var
names(dat_cond_d)
dat_cond_d_log = dat_cond_d
dat_cond_d_log[,4:18] = log10(dat_cond_d_log[,4:18])
dat_cond_d_log[,4:18] = scale(dat_cond_d_log[,4:18])
sum(is.nan(dat_cond_d_log[,4:18]))
sum(is.na(dat_cond_d_log[,4:18]))
range(dat_cond_d_log[,4:18], na.rm = T)

# Pull out only response var
names(dat_cond_d_log)
# AB00, AT07, GV01, HO00, MC06, RG01, RS02,
# EFJ, RED, RSA, & RSAW
dat_dep_d <- t(dat_cond_d_log[,c(4:10,12,15:17)])
row.names(dat_dep_d)

# Make covariate inputs
# without short ts sites:
dat_cov_d <- dat_cond_d_log[,c(19:25, 27, 30:32, # precip
                               34:42, 44, 45, 48, 49, 47, 52, 50, 51)] # fire (reordering to match c matrix below)

# Make new interaction term
dat_cov_d <- dat_cov_d %>%
  mutate(interact_AB00 = cumulative_precip_mm_AB00 * (AB00_Tea_d_AB00 + AB00_Jesusita_d_AB00),
         interact_AT07 = cumulative_precip_mm_AT07 * AT07_Jesusita_d_AT07,
         interact_GV01 = cumulative_precip_mm_GV01 * GV01_Gaviota_d_GV01,
         interact_HO00 = cumulative_precip_mm_HO00 * (HO00_Gaviota_d_HO00 + HO00_Sherpa_d_HO00),
         interact_MC06 = cumulative_precip_mm_MC06 * (MC06_Tea_d_MC06 + MC06_Jesusita_d_MC06),
         interact_RG01 = cumulative_precip_mm_RG01 * RG01_Gaviota_d_RG01,
         interact_RS02 = cumulative_precip_mm_RS02 * (RS02_Tea_d_RS02 + RS02_Jesusita_d_RS02),
         interact_EFJ = cumulative_precip_mm_EFJ * (EFJ_Thompson_d_EFJ + EFJ_Conchas_d_EFJ),
         interact_RED = cumulative_precip_mm_RED * RED_Thompson_d_RED,
         interact_RSA = cumulative_precip_mm_RSA * RSA_Conchas_d_RSA,
         interact_RSAW = cumulative_precip_mm_RSAW * (RSAW_Thompson_d_RSAW + RSAW_Conchas_d_RSAW))

dat_cov_d <- t(scale(dat_cov_d))
row.names(dat_cov_d)
sum(is.nan(dat_cov_d)) # RG01_Sherpa causing problems so removed it from above
sum(is.na(dat_cov_d))

#### make C matrix

# without short ts sites:
CC <- matrix(list( 
  # precip by site
  "AB00_precip",0,0,0,0,0,0,0,0,0,0,
  0,"AT07_precip",0,0,0,0,0,0,0,0,0,
  0,0,"GV01_precip",0,0,0,0,0,0,0,0,
  0,0,0,"HO00_precip",0,0,0,0,0,0,0,
  0,0,0,0,"MC06_precip",0,0,0,0,0,0,
  0,0,0,0,0,"RG01_precip",0,0,0,0,0,
  0,0,0,0,0,0,"RS02_precip",0,0,0,0,
  0,0,0,0,0,0,0,"EFJ_precip", 0,0,0,
  0,0,0,0,0,0,0,0,"RED_precip", 0,0,
  0,0,0,0,0,0,0,0,0,"RSA_precip", 0,
  0,0,0,0,0,0,0,0,0,0,"RSAW_precip",
  # fires by site
  "AB00_Tea",0,0,0,0,0,0,0,0,0,0,
  "AB00_Jesusita",0,0,0,0,0,0,0,0,0,0,
  0,"AT07_Jesusita",0,0,0,0,0,0,0,0,0,
  0,0,"GV01_Gaviota",0,0,0,0,0,0,0,0,
  0,0,0,"HO00_Gaviota",0,0,0,0,0,0,0,
  0,0,0,"HO00_Sherpa",0,0,0,0,0,0,0,
  0,0,0,0,"MC06_Tea",0,0,0,0,0,0,
  0,0,0,0,"MC06_Jesusita",0,0,0,0,0,0,
  0,0,0,0,0,"RG01_Gaviota",0,0,0,0,0,
  0,0,0,0,0,0,"RS02_Tea",0,0,0,0,
  0,0,0,0,0,0,"RS02_Jesusita",0,0,0,0,
  0,0,0,0,0,0,0,"EFJ_Thompson",0,0,0,
  0,0,0,0,0,0,0,"EFJ_Conchas",0,0,0,
  0,0,0,0,0,0,0,0,"RED_Thompson",0,0,
  0,0,0,0,0,0,0,0,0,"RSA_Conchas",0,
  0,0,0,0,0,0,0,0,0,0,"RSAW_Thompson",
  0,0,0,0,0,0,0,0,0,0,"RSAW_Conchas",
  # interaction terms by site
  "AB00_interact",0,0,0,0,0,0,0,0,0,0,
  0,"AT07_interact",0,0,0,0,0,0,0,0,0,
  0,0,"GV01_interact",0,0,0,0,0,0,0,0,
  0,0,0,"HO00_interact",0,0,0,0,0,0,0,
  0,0,0,0,"MC06_interact",0,0,0,0,0,0,
  0,0,0,0,0,"RG01_interact",0,0,0,0,0,
  0,0,0,0,0,0,"RS02_interact",0,0,0,0,
  0,0,0,0,0,0,0,"EFJ_interact", 0,0,0,
  0,0,0,0,0,0,0,0,"RED_interact", 0,0,
  0,0,0,0,0,0,0,0,0,"RSA_interact", 0,
  0,0,0,0,0,0,0,0,0,0,"RSAW_interact"),11,39)

# Model setup
mod_list <- list(
  ### inputs to process model ###
  B = "diagonal and unequal",
  U = "zero",
  C = CC, 
  c = dat_cov_d,
  Q = "diagonal and unequal", 
  ### inputs to observation model ###
  Z='identity', 
  A="zero",
  D="zero" ,
  d="zero",
  R = "zero", 
  ### initial conditions ###
  #x0 = matrix("x0"),
  V0="zero" ,
  tinitx=0
)

# Fit model

# fit BFGS with priors
kemfit <- MARSS(y = dat_dep_d, model = mod_list,
                control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) 

fit <- MARSS(y = dat_dep_d, model = mod_list,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit$par)

# # fit EM by itself
# fit <- MARSS(y = dat_dep, model = mod_list,
#                 control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) 

# export model fit
saveRDS(fit, 
        file = "data_working/marss_test_run/fit_071122_11state_cond_decay_winteraction_mBFGS.rds")

### DIAGNOSES ###

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit[["errors"]]
# NULL - Yay!

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit <- MARSSparamCIs(fit)
# better to do parametric/non-parametric bootstrapping once model is decided upon
# Maybe increase to over 100 boots, 100 is standard
# est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)

saveRDS(est_fit, 
        "data_working/marss_test_run/CIs_fit_071122_11state_cond_decay_winteraction_mBFGS.rds")

# formatting confidence intervals into dataframe
CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)

### Plot Results for All Sites ###

# First, create dataset of all outputs
# This works for HO00 alone
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])

# Now to iterate over all sites
my_list <- c("AB00","AT07","GV01","HO00","MC06","RG01","RS02","EFJ","RED","RSA","RSAW")

# Create an empty list for things to be sent to
datalist = list()

for (i in my_list) { # for every site in the list
  df <- rbind(CIs_fit[grepl(i, CIs_fit$parm),]) # create a new dataset
  df$i <- i  # remember which site produced it
  datalist[[i]] <- df # add it to a list
}

CIs_fit_ed <- bind_rows(datalist) %>% # bind all rows together
  dplyr::rename(Site = i) %>%
  #rename(Parameter = parm) %>%# rename site column
  mutate(Parameter = factor(parm, levels = c(# relevel parameters
                                             "AB00_precip", "AT07_precip", "GV01_precip",
                                             "HO00_precip", "MC06_precip", "RG01_precip",
                                             "RS02_precip", "EFJ_precip", 
                                             "RED_precip", "RSA_precip", "RSAW_precip",
                                             "AB00_Tea", "AB00_Jesusita", "AT07_Jesusita",
                                             "GV01_Gaviota", "HO00_Gaviota", "HO00_Sherpa",
                                             "MC06_Tea", "MC06_Jesusita", "RG01_Gaviota",
                                             "RS02_Tea", "RS02_Jesusita", "EFJ_Thompson",
                                             "EFJ_Conchas", "RED_Thompson", "RSA_Conchas",
                                             "RSAW_Thompson", "RSAW_Conchas",
                                             "AB00_interact", "AT07_interact", "GV01_interact",
                                             "HO00_interact", "MC06_interact", "RG01_interact",
                                             "RS02_interact", "SP02_interact", "EFJ_interact",
                                             "RED_interact", "RSA_interact", "RSAW_interact"
                                             )))

# plot results
(RESULTS_ALL_d <- ggplot(CIs_fit_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 07/11/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + 
    facet_wrap(.~Site, scales = "free"))

# Adding labels for plotting purposes
CIs_fit_ed2 = CIs_fit_ed[!(CIs_fit_ed$Site=="RSA" & CIs_fit_ed$Parameter=="RSAW_precip"),] 
CIs_fit_ed2 = CIs_fit_ed2[!(CIs_fit_ed2$Site=="RSA" & CIs_fit_ed2$Parameter=="RSAW_Thompson"),] 
CIs_fit_ed2 = CIs_fit_ed2[!(CIs_fit_ed2$Site=="RSA" & CIs_fit_ed2$Parameter=="RSAW_Conchas"),] 
CIs_fit_ed2 = CIs_fit_ed2[!(CIs_fit_ed2$Site=="RSA" & CIs_fit_ed2$Parameter=="RSAW_interact"),] 
CIs_fit_ed2$region = c(rep("Coastal California",25),rep("Subalpine New Mexico",14))

(RESULTS_ALL_d <-ggplot(CIs_fit_ed2, aes(Parameter, Est., color=region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=.7) +
    geom_point(position=position_dodge(width=0.3), size=5) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 07/11/2022\n11 state - 4-year fire decay term") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(vars(region, Site), scales = "free"))

# ggsave("figures/MARSS_11states_cond_precip_fire_decay_winteract_071122.png",
#        width = 40,
#        height = 20,
#        units = "cm")

## Script for diagnoses ###

dat = dat_dep_d
time = c(1:ncol(dat_dep_d))
resids <- MARSSresiduals(fit)

# Save resids dataframe
saveRDS(resids, file = "data_working/marss_test_run/all_resids_071122_11state_cond_decay_winteraction_mBFGS.rds")

kf = print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
# No C matrix
# Need to be sure the covariates we've added are *actually* explaining
# the variation we are seeing
# Should have a better AIC score in model above

mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)

# fitting null model
# Note - if we change the structure of the model above, make sure that
# the same code is used to run the null models here below.
null.kemfit <- MARSS(y = dat_dep_d, model = mod_list_null,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit <- MARSS(y = dat_dep_d, model = mod_list_null,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit$par)

bbmle::AICtab(fit, null.fit)
# dAIC- delta AIC
# 0.0 = always the value for the lowest model AIC
#           dAIC df
# fit        0.0 72
# null.fit 398.1 33
# RESULT: covar model is better than null

### Plot response vars ###
# par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
# plot(dat_dep[1,], type="o")
# plot(dat_dep[2,], type="o")
# plot(dat_dep[3,], type="o")
# plot(dat_dep[4,], type="o")
# plot(dat_dep[5,], type="o")
# plot(dat_dep[6,], type="o")
# plot(dat_dep[7,], type="o")
# plot(dat_dep[8,], type="o")
# plot(dat_dep[8,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  #forecast::Acf(resids$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  forecast::Acf(resids$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}

# State residuals don't have any patterns that jump out.

# Error when trying to plot model residuals:
# Error in plot.window(...) : need finite 'ylim' values

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  # qqnorm(resids$model.residuals[i,], main=paste(i, "model residuals"),
  #        pch=16,
  #        xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[i,])[1]))
  # qqline(resids$model.residuals[i,])
  qqnorm(resids$state.residuals[i,], main=paste(i, "state residuals"), pch=16,
         xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[i,])[1]))
  qqline(resids$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# State residuals aren't great, but slightly better at SB.
# similar to results above with binary fire term.

# Error when trying to plot model residuals:
# Error in shapiro.test(resids$model.residuals[i, ]) : # all 'x' values are identical

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

#### Scenario 2 : catchments in two ecoregions #### 

# Pull out response var
# If needed, here are mean conductivity values per location
dat_dep2 <- dat_cond_log %>%
  mutate(avg_cond_SB = (mean_cond_uScm_AB00 + mean_cond_uScm_AT07 + mean_cond_uScm_GV01 + mean_cond_uScm_HO00 + mean_cond_uScm_MC06 + mean_cond_uScm_RG01 + mean_cond_uScm_RS02)/7,
         avg_cond_NM = (mean_cond_uScm_EFJ + mean_cond_uScm_RED + mean_cond_uScm_RSA + mean_cond_uScm_RSAW)/4)

names(dat_dep2)
dat_dep2 <- dat_dep2[,c(55:56)] # conductivity
dat_dep2 <- t(scale(dat_dep2))
row.names(dat_dep2)

# AB00, AT07, GV01, HO00, MC06, RG01, RS02,
# EFJ, RED, RSA, & RSAW

# Make covariate inputs
# For now, using average precip and cumulative fire
# so number of watersheds affected by fires is incorporated (need to run by group...)
dat_cond_log2 <- dat_cond_log %>%
  mutate(avg_precip_SB = (cumulative_precip_mm_AB00 + cumulative_precip_mm_AT07 + cumulative_precip_mm_GV01 + cumulative_precip_mm_HO00 + cumulative_precip_mm_MC06 + cumulative_precip_mm_RG01 + cumulative_precip_mm_RS02)/7,
         avg_precip_NM = (cumulative_precip_mm_EFJ +  cumulative_precip_mm_RED + cumulative_precip_mm_RSA + cumulative_precip_mm_RSAW)/4,
         c_fire_SB = (AB00_Tea_AB00 + AB00_Jesusita_AB00 + AT07_Jesusita_AT07 + GV01_Gaviota_GV01 + HO00_Gaviota_HO00 + HO00_Sherpa_HO00 + MC06_Tea_MC06 + MC06_Jesusita_MC06 + RG01_Gaviota_RG01 + RG01_Sherpa_RG01 + RS02_Tea_RS02 + RS02_Jesusita_RS02),
         c_fire_NM = (RED_Thompson_RED + EFJ_Thompson_EFJ + EFJ_Conchas_EFJ + RSAW_Thompson_RSAW + RSAW_Conchas_RSAW + RSA_Conchas_RSA))

names(dat_cond_log2)
dat_cov2 <- dat_cond_log2[,c(2:3, # seasonal
                            55:58)] # precip + fire
dat_cov2 <- t(scale(dat_cov2))
row.names(dat_cov2)

CC2 <- matrix(list(# season 1
  "Season1", "Season1",
  # season 2
  "Season2", "Season2",
  # precip
  "avg_precip_SB", 0, 
  0, "avg_precip_NM",
  # fire
  "c_fire_SB", 0,
  0, "c_fire_NM"), 2, 6)

## Attempt #2:
# Use dat_dep and dat_cov from 11 state model.

CC2 <- "unconstrained"

# Note - Need to check about altering the R component of the matrix.

# Make covariate inputs
# without short ts sites:
dat_s <- dat_cond_log[,c(2:3)] # seasonal covariates
dat_s <- t(scale(dat_s))
row.names(dat_s)

# Model setup
mod_list2 <- list(
  # tinitx = "zero", # setting initial state value to time = 0
  B = "identity",
  U = "zero", # zero: does NOT allow a drift term in process model to be estimated # removing due to lack of anticipated monotonic trend
  C = CC2, # see new matrix above
  c = dat_s, #dat_cov,
  Q = "diagonal and unequal", # diagonal and unequal: allows for and estimates the covariance matrix of process errors
  Z = matrix(c(1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1), nrow = nrow(dat_dep), ncol = 2), # number of estimated state processes = 1
  A = "zero",
  R = "diagonal and equal" # allows for and estimates the covariance matrix of observations errors (may want to provide a number for this from method precision etc if possible) - changed to "zero" on 11/22 to "turn off" observation error
)

# fit BFGS with priors
kemfit2 <- MARSS(y = dat_dep, model = mod_list2,
                control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) 

fit2 <- MARSS(y = dat_dep, model = mod_list2,
             control = list(maxit = 5000), method = "BFGS", inits=kemfit2$par)

# # fit EM by itself
# fit2 <- MARSS(y = dat_dep, model = mod_list2,
#                 control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) 

# export model fit
#saveRDS(fit2, file = "data_working/marss_test_run/fit_022522_2states.rds")

### DIAGNOSES ###

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit2[["errors"]]
# NULL - Yay!

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit2 <- MARSSparamCIs(fit2)
# better to do parametric/non-parametric bootstrapping once model is decided upon
# Maybe increase to over 100 boots, 100 is standard
# est = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 100, silent=F)

#saveRDS(est_fit2, "data_working/marss_test_run/CIs_fit_022522_2state_cond.rds")

# formatting confidence intervals into dataframe
CIs_fit2 = cbind(
  est_fit2$par$U,
  est_fit2$par.lowCI$U,
  est_fit2$par.upCI$U)
CIs_fit2 = as.data.frame(CIs_fit2)
names(CIs_fit2) = c("Est.", "Lower", "Upper")
CIs_fit2$parm = rownames(CIs_fit2)
CIs_fit2[,1:3] = round(CIs_fit2[,1:3], 3)

### Plot Results for All Sites ###

# First, create dataset of all outputs for two states
CIs_SB <- rbind(CIs_fit2[c(1,3),]) %>%
  mutate(Region = "Santa Barbara")
CIs_NM <- rbind(CIs_fit2[c(2,4),]) %>%
  mutate(Region = "Valles Caldera")
  
CIs_fit2_ed <- rbind(CIs_SB, CIs_NM) %>%
  mutate(Parameter = factor(parm))

# plot results
(RESULTS_ALL <- ggplot(CIs_fit2_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 06/06/2022\n2 State (Season Only)") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(.~Region, scales = "free"))

# Adding labels for plotting purposes
CIs_fit2_ed$region <- c(rep("Coastal California",4),rep("Subalpine New Mexico",4))

(RESULTS_ALL2 <-ggplot(CIs_fit2_ed, aes(Parameter, Est., color=region)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=.7) +
    geom_point(position=position_dodge(width=0.3), size=5) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 02/25/2022") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm")) + # need to play with margins to make it all fit
    facet_wrap(vars(region), scales = "free"))

# ggsave("figures/MARSS_2states_cond_precip_fire_060622.png",
#        width = 20,
#        height = 6,
#        units = "cm")

## Script for diagnoses ###

dat2 = dat_dep
time2 = c(1:ncol(dat_dep))
resids2 <- MARSSresiduals(fit2)
kf2 = print(fit2, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
# No C matrix
# Need to be sure the covariates we've added are *actually* explaining
# the variation we are seeing
# Should have a better AIC score in model above

# Future to-do : make new model with only seasonal effects, and compare
# to null model to see how it performs. It may not be explaining much
# more variability than the null model and making the model above overly
# complex. Or could create model with covariates minus seasonal and compare.
# NOTE - double check to be sure this structure is correct for 2 state option.
mod_list_null2 <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)

# fitting null model
# Note - if we change the structure of the model above, make sure that
# the same code is used to run the null models here below.
null.kemfit2 <- MARSS(y = dat_dep, model = mod_list_null2,
                     control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit2 <- MARSS(y = dat_dep, model = mod_list_null2,
                  control = list(maxit = 5000), method = "BFGS", inits=null.kemfit2$par)

bbmle::AICtab(fit2, null.fit2)
# dAIC- delta AIC
# 0.0 = always the value for the lowest model AIC
#           dAIC df
# fit        0.0 11
# null.fit  71.5 33
# RESULT: covar model is better than null

### Plot response vars ###
# par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
# plot(dat_dep[1,], type="o")
# plot(dat_dep[2,], type="o")
# plot(dat_dep[3,], type="o")
# plot(dat_dep[4,], type="o")
# plot(dat_dep[5,], type="o")
# plot(dat_dep[6,], type="o")
# plot(dat_dep[7,], type="o")
# plot(dat_dep[8,], type="o")
# plot(dat_dep[8,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  forecast::Acf(resids2$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  #forecast::Acf(resids2$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}

# NM may have some yearly patterns in state residuals we should check...

# SB data has some autocorrelation a few months out for model residuals,
# but no patterns that seem worrisome overall.

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  qqnorm(resids2$model.residuals[i,], main=paste(i, "model residuals"),
         pch=16,
         xlab=paste("shapiro test: ", shapiro.test(resids2$model.residuals[i,])[1]))
  qqline(resids2$model.residuals[i,])
  # qqnorm(resids2$state.residuals[i,], main=paste(i, "state residuals"), pch=16,
  #        xlab=paste("shapiro test: ", shapiro.test(resids2$state.residuals[i,])[1]))
  # qqline(resids2$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# Similar to 11 state model, state residuals aren't great, but slightly better at SB.
# QQ plots even look a bit better for model residuals than they did in 11 state model.

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

#### Scenario 3 : all catchments in a single state #### 

# Use dat_dep and dat_cov from 11 state model.

CC3 <- "unconstrained"

# Note - Need to check about altering the R component of the matrix.

# Make covariate inputs
# without short ts sites:
dat_s <- dat_cond_log[,c(2:3)] # seasonal covariates
dat_s <- t(scale(dat_s))
row.names(dat_s)

# Model setup
mod_list3 <- list(
  # tinitx = "zero", # setting initial state value to time = 0
  B = "identity",
  U = "zero", # zero: does NOT allow a drift term in process model to be estimated # removing due to lack of anticipated monotonic trend
  C = CC3, # see new matrix above
  c = dat_s, #dat_cov,
  Q = "diagonal and unequal", # diagonal and unequal: allows for and estimates the covariance matrix of process errors
  Z = matrix(1, nrow = nrow(dat_dep), ncol = 1), # number of estimated state processes = 1
  A = "zero",
  R = "diagonal and equal" # allows for and estimates the covariance matrix of observations errors (may want to provide a number for this from method precision etc if possible) - changed to "zero" on 11/22 to "turn off" observation error
)

# Fit model
kemfit3 <- MARSS(y = dat_dep, model = mod_list3,
                 control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) 

fit3 <- MARSS(y = dat_dep, model = mod_list3,
              control = list(maxit = 5000), method = "BFGS", inits = kemfit3$par)

# fit3 <- MARSS(y = dat_dep, model = mod_list3,
#              control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

# export model fit
#saveRDS(fit3, file = "data_working/marss_test_run/fit_022522_1state.rds")

# Compare 1 and 2 state models with initial 11 state model
m11 <- fit$AIC # 2350 - appears to be the winner!
m2 <- fit2$AIC # 2571
m1 <- fit3$AIC # 2691

### DIAGNOSES ###

## check for hidden errors
# some don't appear in output in console
# this should print all of them out, those displayed and those hidden
fit3[["errors"]]
# NULL - Yay!

### Plot coef and coef estimates ###
## estimates
# hessian method is much fast but not ideal for final results
est_fit3 <- MARSSparamCIs(fit3)
# better to do parametric/non-parametric bootstrapping once model is decided upon
# Maybe increase to over 100 boots, 100 is standard
# est_fit3 <- MARSSparamCIs(fit3, method = "parametric", alpha = 0.05, nboot = 100, silent=F)

#saveRDS(est_fit3, "data_working/marss_test_run/CIs_fit_022522_1state_cond.rds")

# formatting confidence intervals into dataframe
CIs_fit3 = cbind(
  est_fit3$par$U,
  est_fit3$par.lowCI$U,
  est_fit3$par.upCI$U)
CIs_fit3 = as.data.frame(CIs_fit3)
names(CIs_fit3) = c("Est.", "Lower", "Upper")
CIs_fit3$parm = rownames(CIs_fit3)
CIs_fit3[,1:3] = round(CIs_fit3[,1:3], 3)

### Plot Results for All Sites ###

CIs_fit3_ed <- CIs_fit3 %>%
  mutate(Parameter = factor(parm))

# plot results
(RESULTS_ALL3 <- ggplot(CIs_fit3_ed, aes(Parameter, Est.)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
    geom_point(position=position_dodge(width=0.3), size=2) + 
    theme_bw()+
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text = element_text(size = 8)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    coord_flip() +
    labs(y = "",
         title = "Sp. Conductivity MARSS modeling results - 06/06/2022\n1 State (Season Only)") +
    theme(plot.margin=unit(c(.2,.2,.05,.05),"cm"))) # need to play with margins to make it all fit

# ggsave("figures/MARSS_1state_cond_precip_fire_060622.png",
#        width = 10,
#        height = 6,
#        units = "cm")

## Script for diagnoses ###

dat3 = dat_dep
time3 = c(1:ncol(dat_dep))
resids3 <- MARSSresiduals(fit3)
kf3 = print(fit2, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
# No C matrix
# Need to be sure the covariates we've added are *actually* explaining
# the variation we are seeing
# Should have a better AIC score in model above

# Future to-do : make new model with only seasonal effects, and compare
# to null model to see how it performs. It may not be explaining much
# more variability than the null model and making the model above overly
# complex. Or could create model with covariates minus seasonal and compare.
# NOTE - double check to be sure this structure is correct for 1 state option.
mod_list_null3 <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity",
  A = "zero",
  R = "zero" 
)

# fitting null model
# Note - if we change the structure of the model above, make sure that
# the same code is used to run the null models here below.
null.kemfit3 <- MARSS(y = dat_dep, model = mod_list_null3,
                      control = list(maxit= 100, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

null.fit3 <- MARSS(y = dat_dep, model = mod_list_null3,
                   control = list(maxit = 5000), method = "BFGS", inits=null.kemfit3$par)

bbmle::AICtab(fit3, null.fit3)
# dAIC- delta AIC
# 0.0 = always the value for the lowest model AIC
#           dAIC df
# fit        0.0 11
# null.fit  48.6  8
# RESULT: covar model is better than null

### Plot response vars ###
# par(mfrow=c(4,2),oma = c(0, 0, 2, 0))
# plot(dat_dep[1,], type="o")
# plot(dat_dep[2,], type="o")
# plot(dat_dep[3,], type="o")
# plot(dat_dep[4,], type="o")
# plot(dat_dep[5,], type="o")
# plot(dat_dep[6,], type="o")
# plot(dat_dep[7,], type="o")
# plot(dat_dep[8,], type="o")
# plot(dat_dep[8,], type="o")

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  #forecast::Acf(resids3$model.residuals[i,], main=paste(i, "model residuals"), na.action=na.pass, lag.max = 24)
  forecast::Acf(resids3$state.residuals[i,], main=paste(i, "state residuals"), na.action=na.pass, lag.max = 24)
  mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
}

# May have some patterns in model residuals we should check...

# State residuals look fine.

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
for(i in c(1:12)){
  # qqnorm(resids3$model.residuals[i,], main=paste(i, "model residuals"),
  #        pch=16,
  #        xlab=paste("shapiro test: ", shapiro.test(resids3$model.residuals[i,])[1]))
  # qqline(resids3$model.residuals[i,])
  qqnorm(resids3$state.residuals[i,], main=paste(i, "state residuals"), pch=16,
         xlab=paste("shapiro test: ", shapiro.test(resids3$state.residuals[i,])[1]))
  qqline(resids3$state.residuals[i,])
  mtext("Are resids normal?", outer = TRUE, cex = 1.5)
}

# QQ plots look even better for model residuals than they did in 2 state model.

# End of script.
