# README -----------------------------------------------------------------------

# Calls raw SQL queries to extract themed, formatted data from the SQLite
# firearea database (firearea.db). These queries do not access spatial data.

# Users should set the database location to the path where the firearea.db
# resides on their machine.


# source SQL queries -----------------------------------------------------------

source("firearea_database_queries.R")


# database utilities -----------------------------------------------------------

## connect to database
firearea_db <- DBI::dbConnect(RSQLite::SQLite(), "path-to-file/firearea.db")

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
