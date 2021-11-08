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
library(calecopal)

# Load data
map_data <- read_csv("data_raw/sbc_sites_stream_hydro.csv")

# Edit only for sites being used
sfa <- c("AB00", "AT07", "BC02", "DV01", "GV01", "HO00", "MC00", "RG01", "RS02", "SP02", "TO02", "GV202", "HO201", "RG202", "TECA", "GLAN", "SMPA", "GOFS", "GORY", "CAWTP", "SBEB", "ELDE")

map_df1 <- map_data %>%
  filter(sitecode %in% sfa) %>%
  mutate(Type = factor(case_when(Data == "Stream Chemistry" ~ "Stream",
                                 TRUE ~ "Precipitation"))) %>%
  mutate(PairID = factor(case_when(sitecode == "GV202" | sitecode == "GV01" ~ "A",
                                   sitecode == "HO201" | sitecode == "HO00" ~ "B",
                                   sitecode == "RG202" | sitecode == "RG01" ~ "C",
                                   sitecode == "TECA" | sitecode == "TO02" ~ "D",
                                   sitecode == "GLAN" | sitecode == "BC02" ~ "E",
                                   sitecode == "SMPA" | sitecode == "SP02" ~ "F",
                                   sitecode == "GOFS" | sitecode == "DV01" ~ "G",
                                   sitecode == "GORY" | sitecode == "AT07" ~ "H",
                                   sitecode == "CAWTP" | sitecode == "AB00" ~ "I",
                                   sitecode == "SBEB" | sitecode == "MC00" ~ "J",
                                   sitecode == "ELDE" | sitecode == "RS02" ~ "K")))

map_df <- map_data %>%
  mutate(lon = Lon) %>%
  mutate(lat = Lat) %>%
  mutate(DataType = factor(Data))

map_df2 <- map_df1 %>%
  mutate(lon = Lon) %>%
  mutate(lat = Lat)

# Create data sf object
map_sf <- st_as_sf(map_df,
                   coords = c("lon", "lat"),
                   remove = F,
                   crs = 4326) # WGS 84 projection

map_sf2 <- st_as_sf(map_df2,
                   coords = c("lon", "lat"),
                   remove = F,
                   crs = 4326) # WGS 84 projection

# Base plot to see how things are looking...
plot(map_sf2$geometry)

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

kelp_pal <- cal_palette(name = "kelp1", n = 11, type = "continuous")

pairsmap <- ggmap(sb_basemap) + # base google maps tile
  geom_point(data = map_sf2, aes(x = lon, y = lat, 
                                 color = PairID, shape = Type), 
             size = 5,
             inherit.aes = FALSE) + # adds points
  scale_color_manual(values = kelp_pal) +
  #geom_text(data = map_sf, aes(label = sitecode)) +
  ggspatial::annotation_north_arrow(location = "tr") + # adds compass
  ggspatial::annotation_scale() + # adds scale
  geom_text(x = -120, y = 34.40504, label = "Santa Barbara Channel", color = "gray40", size = 4, fontface = "italic") +
  geom_text(x = -119.95, y = 34.5, label = "Santa Ynez Mountains", color = "gray10", size = 4, fontface = "italic") +
  labs(x = "Longitude (WGS84)",
       y = "Latitude") +
  theme_bw() +
  theme(legend.position = "none") +
  coord_sf(crs = st_crs(4326))

pairsmap

# Export map.

# ggsave(pairsmap,
#        filename = "figures/crass_sbc_sites_plus_precip.png",
#        width = 30,
#        height = 15,
#        units = "cm"
#        )

# End of R script.
