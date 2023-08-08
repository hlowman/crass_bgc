# Valles Caldera National Preserve (vcnp) Stream Chemistry Assembly
# November 8, 2021
# Betsy Summers following Heili Lowman's workflow 

# The following script will assemble the stream chemistry datasets available for the vcnp stream sites into a single data file for use in the MARSS analysis for the CRASS project.
# Note that time stamp of solute collection is roughly monthly so we don't need to take monthly means of solutes like done with the SBC solute dataset. 

# Load packages
library(plyr) # needs to be loaded prior to dplyr
library(tidyverse) # contains dplyr
library(lubridate)
library(here)
library(naniar)

# Load dataset
chem_reg_vcnp <- read_csv("data_raw/VCNP_grab.csv")
str(chem_reg_vcnp)
chem_reg_vcnp$date <-  as.Date(chem_reg_vcnp$Date,'%Y/%m/%d') #convert your variable to Date:

str(chem_full_ed_vcnp)

# Rename and select variables to match with sbc data set. 
# note that units are mg/L and not uM
chem_reg_vcnp_names <- chem_reg_vcnp %>%
  dplyr::rename(site_code = "Stream",
         DTG = "Date + Time", # DateTimeGroup = DTG
         nh4_mgL = "Ammonia as N", # units mg/L,
         nO2_nO3_mgL = "Nitrate + Nitrite as N", # mg/L
         po4_mgL = "Phosphate", # mg/L
         tkn_mgL = "TKN",
         tds = "TDS",
         mean_cond_uScm = "Conductivity"
         )
  
# Select variables of interest 
chem_reg_vcnp_names_short <- chem_reg_vcnp_names %>%
  select(site_code, Date, DTG, nh4_mgL, nO2_nO3_mgL, po4_mgL, tkn_mgL, tds, mean_cond_uScm)

str(chem_reg_vcnp_names_short)

# Make Year and Month variable for data assembly purposes
chem_full_ed_vcnp <- chem_reg_vcnp_names_short %>% # can't remove records because in wide format
  #mutate(Date = as.Date(Date, format="%Y-%m-%d")) %>%
  mutate(DateTime = ymd_hms(DTG)) %>% # format dates
  mutate(Year = year(DateTime), Month = month(DateTime))

str(chem_full_ed_vcnp)

# And export for MARSS script
#saveRDS(chem_full_ed_vcnp, "data_working/VCNPchem_edited_110821.rds")
saveRDS(chem_full_ed_vcnp, "data_working/VCNPchem_edited_120521.rds")

# End of script.
