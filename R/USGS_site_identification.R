## USGS site and parameter identification 
## This code is to identify which sites exist in California and Nevada based on the surfacewater parameter list 
## 8/13/21
## Betsy Summers

#install.packages("dataRetrieval")
library(dataRetrieval)
library(ggplot2)
library(ggsn)
library(sf)
library(dplyr)


# list of parameters can be found here: https://nwis.waterdata.usgs.gov/nwis/pmcodes/
# 00910 = Calcium, water, unfiltered, milligrams per liter as calcium carbonate
# 00653	Phosphate, water, filtered, milligrams per liter as PO4
# 00631	Nitrate plus nitrite, water, filtered, milligrams per liter as nitrogen
# 99135	Total organic carbon, water, in situ, estimated, milligrams per liter

# List parameters of interest
# Note we can include more 
pcode<- c("00910","00653","00631","99135")

# CA sites with specified Parameter(s)
CA_sites <- whatNWISsites(stateCd = "CA", 
                          parameterCd = pcode) # c("00665","00925"))

nrow(CA_sites)

# get Phosphorus data
CA_data <- whatNWISdata(stateCd = "CA",
                        parameterCd = pcode) 

CA_data.st <- CA_data %>% 
  filter(site_tp_cd == "ST")

#filter to get stream water only. Option to add filters based on count and time span
phCA.1 <- CA_sites %>% 
  filter(site_tp_cd == "ST") %>% # sampling sites that are streams only
  filter(dec_lat_va < 38)  # sites that are central to southern california
  # filter(colocated == TRUE) # IF you want sites where all parameters are sampled
  # filter(count_nu > 10) %>% # we need to figure out what is the minimium observations allowed
  # mutate(period = as.Date(end_date) - as.Date(begin_date)) %>%
  # filter(period > 5*365) # condition of at least 5 years of data available
  
  
  # extract data for these filtered sites
  # Not sure if we need to do it in batch mode - are we maxing out row lengths allowed in R?    
  param_CA_data <- readNWISqw(siteNumbers = phCA.1$site_no,
                              parameterCd = pcode)

param_summary_CA <- param_CA_data %>% 
  group_by(site_no) %>% 
  summarize(max = max(result_va, na.rm = TRUE),
            count = n()) %>% 
  ungroup() %>% 
  left_join(attr(param_CA_data, "siteInfo"), 
            by = "site_no")


## View sites on map
usa <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE),
                crs = 4269)

sf_ca <- st_as_sf(param_summary_CA,  # sites with ST data only
                  coords = c("dec_long_va", "dec_lat_va"),
                  crs = 4269)
CA.USGS.sites <- ggplot() +
  geom_sf(data = usa[ usa$ID == "california" ,]) +
  geom_sf(data = sf_ca) + 
  xlab(NULL)+
  ylab(NULL)+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  north(sf_ca, symbol=10, location="bottomleft") +
  scalebar(usa[ usa$ID == "california" ,],
           dist=100, dist_unit="mi", st.size = 3,
           transform=TRUE, model="WGS84")


###### Repeat for other states
## Nevada
# NV sites with specified Parameter(s)
NV_sites <- whatNWISsites(stateCd = "NV", 
                          parameterCd = pcode) 

nrow(NV_sites)

# get Phosphorus data
NV_data <- whatNWISdata(stateCd = "NV",
                        parameterCd = pcode) 

NV_data.st <- NV_data %>% 
  filter(site_tp_cd == "ST")

#filter to get stream water only. Option to add filters based on count and time span
phNV.1 <- NV_sites  %>% 
   filter(site_tp_cd == "ST") #%>% # sampling sites that are streams only
  # filter(dec_lat_va < XX) %>% # sites that fall within a lat boundary
  # filter(colocated == TRUE) # IF you want sites where all parameters are sampled
  # filter(count_nu > 10) %>% # we need to figure out what is the minimium observations allowed
  # mutate(period = as.Date(end_date) - as.Date(begin_date)) %>%
  # filter(period > 5*365) # condition of at least 5 years of data available
  
  
  # extract data for these filtered sites
  # Not sure if we need to do it in batch mode - are we maxing out row lengths allowed in R?    
param_NV_data <- readNWISqw(siteNumbers = phNV.1$site_no,
                              parameterCd = pcode)

param_summary_NV <- param_NV_data %>% 
  group_by(site_no) %>% 
  summarize(max = max(result_va, na.rm = TRUE),
            count = n()) %>% 
  ungroup() %>% 
  left_join(attr(param_NV_data, "siteInfo"), 
            by = "site_no")


## View sites on map
usa <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE),
                crs = 4269)

sf_nv <- st_as_sf(param_summary_NV,  # sites with ST data only
                  coords = c("dec_long_va", "dec_lat_va"),
                  crs = 4269)
NV.USGS.sites <- ggplot() +
  geom_sf(data = usa[ usa$ID == "nevada" ,]) +
  geom_sf(data = sf_nv) + 
  xlab(NULL)+
  ylab(NULL)+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  north(sf_nv, symbol=10, location="bottomleft") +
  scalebar(usa[ usa$ID == "nevada" ,],
           dist=100, dist_unit="mi", st.size = 3,
           transform=TRUE, model="WGS84")


#### Export USGS sites as metadata for lat/long
both_states <- full_join(phCA.1, phNV.1)

# As .csv file
write.csv(both_states, "../data_raw/USGS_sites_metadata.csv") 
