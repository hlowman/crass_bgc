### Merge chem, Q, & fire ###
## Calculate % burn for each catchment*fire
## Identify largest fire
## Identify most recent fire
## Categorize dates to pre/post fire, pre/post largest fire, pre/post most recent fire

## Inputs: 
      # Merged chems & Q: USGS_chem_Q.csv
      # Output of catchments_fires from script 02: catchment area and area burned for each fire

### Fire summary attributes ###
# Is Ig_Date local time or UTC? Assuming UTC here.
catchments_fires$Ig_Date <- as.Date(catchments_fires$Ig_Date, format = "%Y-%m-%d", tz = "UTC")

# Fraction catchment burned in each fire
catchments_fires <- catchments_fires %>% mutate(burn_ext = catchment_area/ws_burn_area_km2) 

# Subset most recent fire each catchment

# Subset largest fire each catchment

# Join most recent and largest

# Join to catchments_summary for cumulative burn extent

### Plot C-Q pre/post burn ###

### Next scripts
## Then apply filters: size, Now add LULC filter & replot C-Q