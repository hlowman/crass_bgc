# README -----------------------------------------------------------------------

# Features the raw SQL queries to extract themed, formatted data from the
# SQLite firearea database (firearea.db). These queries do not access spatial
# data.


# libraries --------------------------------------------------------------------

library(DBI)
library(RSQLite)


# queries ----------------------------------------------------------------------

query_shiny_app_summary <- "
SELECT
  catchments.usgs_site,
  ROUND(catchments.area_km2, 1) AS catchment_area,
  fires_subquery.number_fires,
  fires_subquery.total_area_burned,
  wwtp_subquery.number_wwtp,
  lulc_subquery.dominant_lulc
FROM catchments
LEFT JOIN (
  SELECT
    fires.usgs_site,
    COUNT(fires.usgs_site) AS number_fires,
    ROUND(SUM(fires.ws_burn_area_km2), 1) AS total_area_burned
  FROM fires
  GROUP BY usgs_site
  ) AS fires_subquery ON (fires_subquery.usgs_site = catchments.usgs_site)
LEFT JOIN (
  SELECT
    wwtp.usgs_site,
    COUNT(wwtp.usgs_site) AS number_wwtp
  FROM wwtp
  GROUP BY usgs_site
  ) AS wwtp_subquery ON (wwtp_subquery.usgs_site = catchments.usgs_site)
LEFT JOIN(
  SELECT
    usgs_site,
    MAX(freq),
    type AS dominant_lulc
  FROM lulc
  GROUP BY usgs_site
  ) AS lulc_subquery ON (lulc_subquery.usgs_site = catchments.usgs_site)
;
"


query_catchments_fires <- "
SELECT
  catchments.usgs_site,
  ROUND(catchments.area_km2, 1) AS catchment_area,
  fires.Event_ID,
  fires.Incid_Name,
  fires.Ig_Date,
  fires.total_burn_area_km2,
  fires.ws_burn_area_km2
FROM catchments
LEFT JOIN fires ON (catchments.usgs_site = fires.usgs_site)
;
"


query_catchments_lulc <- "
SELECT
  catchments.usgs_site,
  ROUND(catchments.area_km2, 1) AS catchment_area,
  lulc.freq,
  lulc.type
FROM catchments
LEFT JOIN lulc ON (catchments.usgs_site = lulc.usgs_site)
;
"


query_catchments_chemistry <- "
SELECT
  catchments.usgs_site,
  ROUND(catchments.area_km2, 1) AS catchment_area,
  water_chem.*
FROM catchments
LEFT JOIN water_chem ON (catchments.usgs_site = water_chem.usgs_site)
;
"
