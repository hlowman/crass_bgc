# README -----------------------------------------------------------------------

# Calls raw SQL queries detailed in `firearea_database_queries.R` to extract
# themed, formatted data from the SQLite firearea database (firearea.db). These
# queries do not access spatial data.

# Users should set path_to_database to the path of the directory where the
# firearea.db resides on their machine.


# source SQL queries -----------------------------------------------------------

source("firearea_database_queries.R")


# database utilities -----------------------------------------------------------

## connect to database
path_to_database <- "path_to_directory_where_you_have_downloaded_the_database"
firearea_db      <- DBI::dbConnect(RSQLite::SQLite(), here::here(path_to_database, "firearea.db"), shutdown=TRUE)

## get database details
DBI::dbGetInfo(firearea_db)

## disconnect from database (disconnect when finished with connection)
# DBI::dbDisconnect(firearea_db)


# run queries ------------------------------------------------------------------

catchments_summary <- DBI::dbGetQuery(
  conn      = firearea_db,
  statement = query_shiny_app_summary
)

catchments_fires <- DBI::dbGetQuery(
  conn      = firearea_db,
  statement = query_catchments_fires
)

catchments_lulc <- DBI::dbGetQuery(
  conn      = firearea_db,
  statement = query_catchments_lulc
)

catchments_water_chemistry <- DBI::dbGetQuery(
  conn      = firearea_db,
  statement = query_catchments_chemistry
)

catchments_daily_discharge <- DBI::dbGetQuery(
  conn      = firearea_db,
  statement = query_discharge_daily
)
