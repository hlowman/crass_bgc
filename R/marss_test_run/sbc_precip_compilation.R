# SB Precipitation Assembly
# October 22, 2021
# Heili Lowman

# The following script will assemble the precipitation datasets available for the SBC LTER stream sites into a single data file for use in the MARSS analysis for the CRASS project.

# Load packages
library(plyr) # needs to be loaded prior to dplyr
library(tidyverse) # contains dplyr
library(lubridate)
library(here)

# Load datasets
# First, precipitation - both at LTER and SBC County sites
# Data downloaded from: https://sbclter.msi.ucsb.edu/data/catalog/ 
# with "Watershed" and "Hydrology" selected
SBprecip <- ldply(
  list.files(here("data_raw", "Precipitation", "SBCLTER_Precip"), 
             pattern = "csv"), 
  function(filename){
  d <- read_csv(here("data_raw", "Precipitation", "SBCLTER_Precip", filename))
  d$file <- filename
  return(d)
})

# Note: -999 is the "NA" record used by the LTER & the county
# Make edits for data assembly purposes
SBprecip_ed <- SBprecip %>%
  mutate(precip_ed = case_when(precipitation_mm < 0 ~ 0,
                               TRUE ~ precipitation_mm)) %>% # make all negative records 0
  mutate(DateTime = ymd_hms(timestamp_local)) %>% # format dates
  mutate(Year = year(DateTime), Month = month(DateTime)) %>%
  mutate(site = factor(file))

# Calculate cumulative monthly precipitation at all sites.
SBprecip_monthly <- SBprecip_ed %>%
  group_by(site, Year, Month) %>%
  summarize(cumulative_precip_mm = sum(precip_ed, na.rm = TRUE)) %>%
  ungroup()

# Add site abbreviations so easier to filter.
SBprecip_monthly <- SBprecip_monthly %>%
  mutate(sitecode = case_when(site == "BaronRanch262_precip_allyears_2019-09-16.csv" ~ "BARA",
                        site == "BotanicGarden321_precip_allyears_2019-09-16.csv" ~ "BOGA",
                        site == "BuelltonFS233_precip_allyears_2019-09-16.csv" ~ "BUEL",
                        site == "Carpinteria208_precip_allyears_2019-09-16.csv" ~ "CARP208",
                        site == "CarpinteriaUSFS383_precip_allyears_2019-09-16.csv" ~ "CARPUSFS",
                        site == "CaterWTP229_precip_allyears_2019-09-16.csv" ~ "CAWTP",     
                        site == "ColdSprings210_precip_allyears_2019-09-16.csv" ~ "COSP",   
                        site == "CP201_precip_allyears_2019-09-11.csv" ~ "CP201",            
                        site == "DosPueblos226_precip_allyears_2019-09-16.csv" ~ "DOPU",    
                        site == "DoultonTunnel231_precip_allyears_2019-09-16.csv" ~ "DOTU", 
                        site == "EdisonTrail252_precip_allyears_2019-09-16.csv" ~ "EDTR",   
                        site == "EL201_precip_allyears_2019-09-11.csv" ~ "EL201",           
                        site == "EL202_precip_allyears_2019-09-11.csv" ~ "EL202",           
                        site == "ElDeseo255_precip_allyears_2019-09-16.csv" ~ "ELDE",       
                        site == "GaviotaSP301_precip_allyears_2019-09-16.csv" ~ "GASP",     
                        site == "GB201_precip_allyears_2019-09-11.csv" ~ "GB201",           
                        site == "GlenAnnieCanyon309_precip_allyears_2019-09-16.csv" ~ "GLAN",
                        site == "GoletaFireStation440_precip_allyears_2019-09-16.csv" ~ "GOFS",  
                        site == "GoletaRdYard211_precip_allyears_2019-09-16.csv" ~ "GORY",  
                        site == "GoletaWaterDistrict334_precip_allyears_2019-09-16.csv" ~ "GOWD",
                        site == "GV202_precip_allyears_2019-09-11.csv" ~ "GV202",           
                        site == "HO201_precip_allyears_2019-09-11.csv" ~ "HO201",           
                        site == "HO202_precip_allyears_2019-09-11.csv" ~ "HO202",           
                        site == "KTYD227_precip_allyears_2019-09-16.csv" ~ "KTYD",          
                        site == "Montecito325_precip_allyears_2019-09-16.csv" ~ "MO325",    
                        site == "Nojoqui236_precip_allyears_2019-09-16.csv" ~ "NO236",      
                        site == "RanchoSJ389_precip_allyears_2019-09-16.csv" ~ "RASJ",      
                        site == "RefugioPass429_precip_allyears_2019-09-16.csv" ~ "REPA",   
                        site == "RG201_precip_allyears_2019-09-11.csv" ~ "RG201",           
                        site == "RG202_precip_allyears_2019-09-11.csv" ~ "RG202",           
                        site == "RG203_precip_allyears_2019-09-11.csv" ~ "RG203",           
                        site == "RG204_precip_allyears_2019-09-11.csv" ~ "RG204",
                        site == "SanMarcosPass212_precip_allyears_2019-09-16.csv" ~ "SMPA",
                        site == "SBCaltrans335_precip_allyears_2019-09-16.csv" ~ "SBCT",  
                        site == "SBEngBldg234_precip_allyears_2019-09-16.csv" ~ "SBEB",
                        site == "StanwoodFS228_precip_allyears_2019-09-16.csv" ~ "STFS",
                        site == "TecoloteCanyon280_precip_allyears_2019-09-16.csv" ~ "TECA",
                        site == "TroutClub242_precip_allyears_2019-09-16.csv" ~ "TRCL",
                        site == "UCSB200_precip_allyears_2019-09-16.csv" ~ "UCSB"))

# And export for MARSS script
saveRDS(SBprecip_monthly, "data_working/SBprecip_edited_110721.rds")

#### Additional data exploration regarding select sites ####

# Filter for sites to be used in the initial MARSS run
desired <- c("HO201", "RG202", "CAWTP")

SBprecip_filtered <- SBprecip_monthly %>%
  filter(sitecode %in% desired)

# Visualize each data record to examine for gaps
(data_coverage <- SBprecip_filtered %>%
  mutate(Day = 1) %>%
  mutate(Date = make_date(Year, Month, Day)) %>%
  ggplot(aes(x = Date)) +
  geom_bar() +
  facet_wrap(.~sitecode)) # ok, so CAWTP is out until we can deal with that gap

# End of script.
