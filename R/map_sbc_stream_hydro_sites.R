# SBC Map Creation
# Heili Lowman
# 9/29/21

# The following script will create the map to display SBC sites.

# Load packages
library(tidyverse)
library(ggplot2)
library(ggmap)
library(ggrepel)
library(sf)
library(USAboundaries)
library(USAboundariesData)
library(ggspatial)
library(here)

# Load data
map_data <- read_csv("data_raw/sbc_sites_stream_hydro.csv")

map_df <- map_data %>%
  mutate(lon = Lon) %>%
  mutate(lat = Lat)

# Create data sf object
map_sf <- st_as_sf(map_df,
                   coords = c("lon", "lat"),
                   remove = F,
                   crs = 4326) # WGS 84 projection

# Base plot to see how things are looking...
plot(map_sf$geometry)

# create CA county dataset for use in map below using USAboundaries
CA_counties <- us_counties(states = "California")

SB_county <- CA_counties %>%
  filter(name == "Santa Barbara")

# create base terrain map tile

# create bounding box
lats <- c(34.368063, 34.569875)
lons <- c(-120.287304, -119.634991)
bb <- make_bbox(lon = lons, lat = lats, f = 0.05)

sb_basemap <- get_stamenmap(bb, 
                      zoom = 12,
                      maptype = 'terrain-background')

ggmap(sb_basemap)
  
fullmap <- ggmap(sb_basemap) + # base google maps tile
  geom_point(data = map_sf, aes(x = lon, y = lat, fill = Data), 
          size = 4,
          shape = 21,
          inherit.aes = FALSE) + # adds points
  scale_fill_manual(values = c("white", "lightgrey", "grey40", "black")) +
  ggspatial::annotation_north_arrow(location = "tr") + # adds compass due north
  ggspatial::annotation_scale() + # adds scale
  geom_text(x = -120, y = 34.40504, label = "Santa Barbara Channel", color = "gray40", size = 4, fontface = "italic") +
  geom_text(x = -119.95, y = 34.5, label = "Santa Ynez Mountains", color = "gray10", size = 4, fontface = "italic") +
  labs(x = "Longitude (WGS84)",
       y = "Latitude",
       fill = "Available Data") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.background = element_rect(fill = "white", size = 0.5, linetype = "solid")) +
  coord_sf(crs = st_crs(4326))

fullmap

# Export map to desktop.

ggsave(fullmap,
       filename = "figures/crass_sbc_hydro_sites.png",
       width = 30,
       height = 15,
       units = "cm"
       )

# End of R script.
