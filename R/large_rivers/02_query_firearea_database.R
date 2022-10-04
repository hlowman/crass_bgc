# README -----------------------------------------------------------------------

# Calls raw SQL queries detailed in `firearea_database_queries.R` to extract
# themed, formatted data from the SQLite firearea database (firearea.db). These
# queries do not access spatial data.

# Users should set path_to_database to the path of the directory where the
# firearea.db resides on their machine.

## Libraries
library(here)

# source SQL queries -----------------------------------------------------------

source(here("R", "large_rivers", "01_firearea_database_queries.R"))


# database utilities -----------------------------------------------------------

## connect to database
path_to_database <- "../"
firearea_db      <- DBI::dbConnect(RSQLite::SQLite(), here::here(path_to_database, "firearea.db"), shutdown=TRUE)
DBI::dbExecute(firearea_db, "PRAGMA foreign_keys = ON ;") # enforce foreign keys

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

catchments_wwtp <- DBI::dbGetQuery(
  conn      = firearea_db,
  statement = query_catchments_wwtp
) |>
dplyr::select(-geometry)
