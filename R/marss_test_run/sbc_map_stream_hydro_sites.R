# SBC Map Creation
# Heili Lowman
# 9/29/21

# The following script will create the map to display SBC sites.
# This script has been edited to examine which sites are closest to which precipitation datasets.

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
  mutate(lat = Lat) %>%
  mutate(DataType = factor(Data))

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
  geom_point(data = map_sf, aes(x = lon, y = lat, color = DataType), 
          size = 2,
          inherit.aes = FALSE) + # adds points
  geom_text(data = map_sf, aes(label = sitecode)) +
  ggspatial::annotation_north_arrow(location = "tr") + # adds compass due north
  ggspatial::annotation_scale() + # adds scale
  geom_text(x = -120, y = 34.40504, label = "Santa Barbara Channel", color = "gray40", size = 4, fontface = "italic") +
  geom_text(x = -119.95, y = 34.5, label = "Santa Ynez Mountains", color = "gray10", size = 4, fontface = "italic") +
  labs(x = "Longitude (WGS84)",
       y = "Latitude") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.background = element_rect(fill = "white", size = 0.5, linetype = "solid")) +
  facet_wrap(.~DataType) +
  coord_sf(crs = st_crs(4326))

fullmap

halfmap <- ggmap(sb_basemap) + # base google maps tile
  geom_point(data = map_sf %>% 
               filter(Data == "Precipitation County" |
                        Data == "Stream Chemistry"), 
             aes(x = lon, y = lat, color = DataType), 
             size = 2,
             inherit.aes = FALSE) + # adds points
  geom_text(data = map_sf %>% 
              filter(Data == "Precipitation County" |
                       Data == "Stream Chemistry"),
            aes(label = sitecode)) +
  ggspatial::annotation_north_arrow(location = "tr") + # adds compass due north
  ggspatial::annotation_scale() + # adds scale
  geom_text(x = -120, y = 34.40504, label = "Santa Barbara Channel", color = "gray40", size = 4, fontface = "italic") +
  geom_text(x = -119.95, y = 34.5, label = "Santa Ynez Mountains", color = "gray10", size = 4, fontface = "italic") +
  labs(x = "Longitude (WGS84)",
       y = "Latitude") +
  theme_bw() +
  theme(legend.position = "none",
        legend.background = element_rect(fill = "white", size = 0.5, linetype = "solid")) +
  facet_wrap(.~DataType, nrow = 2) +
  coord_sf(crs = st_crs(4326))
  

halfmap

# Export map to desktop.

# ggsave(halfmap,
#        filename = "figures/crass_sbc_hydro_sites_4.png",
#        width = 30,
#        height = 15,
#        units = "cm"
#        )

# End of R script.
