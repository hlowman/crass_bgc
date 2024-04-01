#### Marss Stream Distance
##Yvette D. Hastings
##Written: 3/25/2024

##mapping
##load libraries
library(dplyr)
library(sf)
library(ggplot2)
library(ggspatial)
library(raster)
library(riverdist)
library(terra)


# Load Geojson ------------------------------------------------------------
catchment_sbc <- sf::read_sf("Data/sbc_catchments.geojson") ##catchment
catchment_fires_sbc <- sf::read_sf("Data/sbc_fires.geojson") ##fires
catchment_flowline_sbc <- sf::read_sf("Data/sbc_flowlines.geojson") ##flowlines
catchment_sampling_sbc <- sf::read_sf("Data/sbc_pour_points.geojson") ##sampling point

catchment_nm <- sf::read_sf("Data/nm_catchments.geojson")
catchment_fires_nm <- sf::read_sf("Data/nm_fires.geojson") ##fires
catchment_flowline_nm <- sf::read_sf("Data/nm_flowlines.geojson") ##flowlines
catchment_sampling_nm <- sf::read_sf("Data/nm_points.geojson") ##sampling point


# Calculate distance from sampling point to stream point ------------------
##Return sf object with points from streams

##sbc sites
sbc_point_to_stream_dist <- geosphere::dist2Line(p = st_coordinates(catchment_sampling_sbc), 
                             line = st_coordinates(catchment_flowline_sbc)[,1:2]) 


sbc_point_to_stream_dist.sf <- st_as_sf(as.data.frame(sbc_point_to_stream_dist), coords = c('lon', 'lat')) %>%
  st_set_crs(value = 4326) %>%
  mutate(site = c('gav', 'hon', 'bur', 'rat'))

##nm sites
SA_sample <- catchment_sampling_nm[3,]
SA_stream <- catchment_flowline_nm[44,]
Jemez_sample <- catchment_sampling_nm[2,]
Jemez_stream <- catchment_flowline_nm[1,]
RED_sample <-catchment_sampling_nm[1,]
RED_stream <-catchment_flowline_nm[42,]

SA_dist <- geosphere::dist2Line(line = st_coordinates(SA_stream)[,1:2],
                                                p = st_coordinates(SA_sample)) 

Jemez_dist <- geosphere::dist2Line(line = st_coordinates(Jemez_stream)[,1:2],
                                p = st_coordinates(Jemez_sample))

RED_dist <- geosphere::dist2Line(line = st_coordinates(RED_stream)[,1:2],
                                   p = st_coordinates(RED_sample))

nm_point_to_stream_dist <- rbind(SA_dist, Jemez_dist)
nm_point_to_stream_dist <- rbind (nm_point_to_stream_dist, RED_dist)

nm_point_to_stream_dist.sf <- st_as_sf(as.data.frame(nm_point_to_stream_dist), coords = c('lon', 'lat')) %>%
  st_set_crs(value = 4326) %>%
  mutate(site = c('San Antonio (west)', 'Jemez', 'Redondo'))




# Calculate point where stream intersects fire polygon --------------------
##find where streams intersect fire polygons and extract lat/long

##sbc sites
sbc_intersections_lp <- st_intersection(catchment_flowline_sbc, catchment_fires_sbc)

#fire boundary/stream intersect
nhd_17596123_bur_1 <- as.data.frame(st_coordinates(sbc_intersections_lp[1,])) #convert linestring to df
nhd_17596127_bur_2 <- as.data.frame(sbc_intersections_lp[[6]][[2]][[2]]) #convert matrix to df


sbc_polygon_intersect <- data.frame(nhd_id = c('bur Stream Point', 'nhd_17596123_bur_1', 'nhd_17596127_bur_2'),
                                event_id = c('NA', 'CA3447411972520090505', 'CA3447411972520090505'),
                                longitude = c(sbc_point_to_stream_dist[3,2], nhd_17596123_bur_1[35,1], nhd_17596127_bur_2[6,1]),
                                latitude = c(sbc_point_to_stream_dist[3,3], nhd_17596123_bur_1[35,2],nhd_17596127_bur_2[6,2]))

sbc_polygon_intersect_sf <- st_as_sf(sbc_polygon_intersect, coords = c('longitude', 'latitude'), crs = 4326) #set as sf object


##nm sites
nm_intersections_lp <- st_intersection(catchment_flowline_nm, catchment_fires_nm)

#fire boundary/stream intersect
nm_17826286_EFJ <- as.data.frame(st_coordinates(nm_intersections_lp[1,])) #convert linestring to df
nm_17826160_EFJ <- as.data.frame(st_coordinates(nm_intersections_lp[3,])) #convert matrix to df
nm_17827554_EFJ <- as.data.frame(st_coordinates(nm_intersections_lp[4,])) #convert linestring to df
nm_17826190_EFJ <- as.data.frame(st_coordinates(nm_intersections_lp[5,])) #convert linestring to df
nm_17826216_EFJ <- as.data.frame(st_coordinates(nm_intersections_lp[6,])) #convert linestring to df
nm_17826086_EFJ <- as.data.frame(st_coordinates(nm_intersections_lp[7,])) #convert linestring to df
nm_17826188_EFJ <- as.data.frame(st_coordinates(nm_intersections_lp[8,])) #convert linestring to df
nm_17827550_EFJ <- as.data.frame(st_coordinates(nm_intersections_lp[11,])) #convert linestring to df
nm_17827558_EFJ <- as.data.frame(st_coordinates(nm_intersections_lp[12,])) #convert linestring to df
nm_17827556_EFJ <- as.data.frame(st_coordinates(nm_intersections_lp[13,])) #convert linestring to df
nm_17825970_SA <- as.data.frame(st_coordinates(nm_intersections_lp[24,])) #convert linestring to df
nm_17826016_SA <- as.data.frame(st_coordinates(nm_intersections_lp[26,])) #convert linestring to df
nm_17825884_SA <- as.data.frame(st_coordinates(nm_intersections_lp[27,])) #convert linestring to df
nm_17826052_SA <- as.data.frame(st_coordinates(nm_intersections_lp[29,])) #convert linestring to df
nm_17825904_SA <- as.data.frame(st_coordinates(nm_intersections_lp[30,])) #convert linestring to df
nm_17825952_SA <- as.data.frame(st_coordinates(nm_intersections_lp[35,])) #convert linestring to df
nm_17825980_SA <- as.data.frame(st_coordinates(nm_intersections_lp[43,])) #convert linestring to df
nm_17826022_SA <- as.data.frame(st_coordinates(nm_intersections_lp[47,])) #convert linestring to df
nm_17826090_SA <- as.data.frame(st_coordinates(nm_intersections_lp[52,])) #convert linestring to df
nm_17826026_SA <- as.data.frame(st_coordinates(nm_intersections_lp[53,])) #convert linestring to df


nm_polygon_intersect <- data.frame(nhd_id = c('EFJ Stream Point', 'AM stream Point', 'nm_17826286_EFJ', 'nm_17826160_EFJ', 'nm_17827554_EFJ', 'nm_17826190_EFJ',
                                               'nm_17826216_EFJ', 'nm_17826086_EFJ', 'nm_17826188_EFJ', 'nm_17827550_EFJ', 
                                               'nm_17827558_EFJ', 'nm_17827556_EFJ', 'nm_17825970_SA', 'nm_17826016_SA', 'nm_17825884_SA',
                                               'nm_17826052_SA', 'nm_17825904_SA', 'nm_17825952_SA', 'nm_17825980_SA', 'nm_17826022_SA',
                                               'nm_17826090_SA', 'nm_17826026_SA'),
                                    event_id = c('NA', 'NA', 'NM3581210654120110626', 'NM3581210654120110626', 'NM3581210654120110626',
                                                 'NM3581210654120110626', 'NM3581210654120110626', 'NM3581210654120110626','NM3581210654120110626',
                                                 'NM3589210662020130531', 'NM3589210662020130531', 'NM3589210662020130531', 'NM3581210654120110626',
                                                 'NM3581210654120110626', 'NM3581210654120110626', 'NM3581210654120110626', 'NM3581210654120110626',
                                                 'NM3581210654120110626', 'NM3581210654120110626', 'NM3581210654120110626', 'NM3589210662020130531',
                                                 'NM3589210662020130531'),
                                    longitude = c(nm_point_to_stream_dist[2,2], nm_point_to_stream_dist[3,2], nm_17826286_EFJ[15,1], nm_17826160_EFJ[16,1],
                                                  nm_17827554_EFJ[19,1], nm_17826190_EFJ[25,1], nm_17826216_EFJ[12,1], 
                                                  nm_17826086_EFJ[19,1], nm_17826188_EFJ[33,1], nm_17827550_EFJ[6,1], 
                                                  nm_17827558_EFJ[22,1], nm_17827556_EFJ[20,1], nm_17825970_SA[18,1], 
                                                  nm_17826016_SA[14,1], nm_17825884_SA[25,1], nm_17826052_SA[8,1], nm_17825904_SA[44,1],
                                                  nm_17825952_SA[13,1], nm_17825980_SA[15,1], nm_17826022_SA[23,1], nm_17826090_SA[27,1],
                                                  nm_17826026_SA[28,1]),
                                    latitude = c(nm_point_to_stream_dist[2,3], nm_point_to_stream_dist[3,3], nm_17826286_EFJ[15,2], nm_17826160_EFJ[16,2],
                                                 nm_17827554_EFJ[19,2], nm_17826190_EFJ[25,2], nm_17826216_EFJ[12,2], 
                                                 nm_17826086_EFJ[19,2], nm_17826188_EFJ[33,2], nm_17827550_EFJ[6,2], 
                                                 nm_17827558_EFJ[22,2], nm_17827556_EFJ[20,2], nm_17825970_SA[18,2], 
                                                 nm_17826016_SA[14,2], nm_17825884_SA[25,2], nm_17826052_SA[8,2], nm_17825904_SA[44,2],
                                                 nm_17825952_SA[13,2], nm_17825980_SA[15,2], nm_17826022_SA[23,2], nm_17826090_SA[27,2],
                                                 nm_17826026_SA[28,2]))

nm_polygon_intersect_sf <- st_as_sf(nm_polygon_intersect, coords = c('longitude', 'latitude'), crs = 4326) #set as sf object


# Distance along stream lines ---------------------------------------------
##calculate distance along stream lines - https://cran.r-project.org/web/packages/riverdist/vignettes/riverdist_vignette.html
wgs84_UTM <- CRS("+proj=utm +zone=13 +datum=WGS84") ##convert WGS84 to UTM projection

#site sbc - bur is the only site that needs network
catchment_flowline_sbc_bur <- catchment_flowline_sbc[1:9,]
sbc_flow_network <- line2network(sf = catchment_flowline_sbc_bur, reproject=wgs84_UTM)
plot(sbc_flow_network)

topologydots(rivers=sbc_flow_network) #check connectedness; everything looks good

points_xy = as.data.frame(SpatialPoints(cbind(sbc_polygon_intersect[1:3,3],sbc_polygon_intersect[1:3,4]),proj4string=CRS("+proj=longlat")))
points_xy <- as.data.frame(project(as.matrix(points_xy[,c('coords.x1', 'coords.x2')]), 
        "+proj=longlat", "+proj=utm +zone=13 +units=m"))


#converts XY data to stream locations and plots
points_riv <- xy2segvert(x=points_xy$V1, y = points_xy$V2, rivers = sbc_flow_network)
head(points_riv)
# hist(points_riv$snapdist, main = 'Snapping Distance (m)', xlab = 'Distance (m)')

zoomtoseg(seg=c(1,4), rivers = sbc_flow_network) #zoom map
points(sbc_polygon_intersect[1:3,3], sbc_polygon_intersect[1:3,4], pch = 16, col = 'red')
riverpoints(seg = points_riv$seg, vert = points_riv$vert, rivers = sbc_flow_network, pch = 16, col = 'blue')

##calculating distance - need seg and vert from points_riv output
detectroute(start = 1, end = 6, rivers = sbc_flow_network) ##this shows which segments are traversed from start to finish

sample_to_17596123 <- riverdistance(startseg = 1, startvert = 22, endseg = 4, endvert = 34,
                                    rivers = sbc_flow_network, map = TRUE) #returns distance in meters
sample_to_17596127 <- riverdistance(startseg = 1, startvert = 22, endseg = 5, endvert = 25,
                                    rivers = sbc_flow_network, map = TRUE) #returns distance in meters


sbc_distances <- data.frame(site = c('bur', 'bur', 'gav', 'hon', 'rat'),
                        nhdplus_comid = c('17596123', '17596127', NA, NA, NA),
                        event_id = c('CA3447411972520090505', 'CA3447411972520090505', 
                                     'CA3448712019620040605', 'CA3448712019620040605', 'CA3447411972520090505'),
                        sampling_point_to_stream_m = c(sbc_point_to_stream_dist[3,1], sbc_point_to_stream_dist[3,1],
                                                       sbc_point_to_stream_dist[1,1], sbc_point_to_stream_dist[2,1],
                                                       sbc_point_to_stream_dist[4,1]),
                        stream_distance_m = c(sample_to_17596123, sample_to_17596127, 0, 0, 0),
                        distance_sampling_point_to_fire_km = NA)

 
sbc_distances$distance_sampling_point_to_fire_km <- sbc_distances$stream_distance_m/1000



##NM site - Jemez and SA is the only site that needs network
nm_flow_network <- line2network(sf = catchment_flowline_nm, reproject=wgs84_UTM)
plot(nm_flow_network)

topologydots(rivers=nm_flow_network) #check connectedness; everything looks good

points_xy_nm = as.data.frame(SpatialPoints(cbind(nm_polygon_intersect[1:22,3],nm_polygon_intersect[1:22,4]),proj4string=CRS("+proj=longlat")))
points_xy_nm <- as.data.frame(project(as.matrix(points_xy_nm[,c('coords.x1', 'coords.x2')]), 
                                   "+proj=longlat", "+proj=utm +zone=13 +units=m"))


#converts XY data to stream locations and plots
points_riv_nm <- xy2segvert(x=points_xy_nm$V1, y = points_xy_nm$V2, rivers = nm_flow_network)
head(points_riv_nm)
hist(points_riv_nm$snapdist, main = 'Snapping Distance (m)', xlab = 'Distance (m)')

zoomtoseg(seg=c(5,22), rivers = nm_flow_network) #zoom map
points(nm_polygon_intersect[1:22,3], nm_polygon_intersect[1:22,4], pch = 16, col = 'red')
riverpoints(seg = points_riv_nm$seg, vert = points_riv_nm$vert, rivers = nm_flow_network, pch = 16, col = 'blue')

##calculating distance - need seg and vert from points_riv output
detectroute(start = 1, end = 15, rivers = nm_flow_network) ##this shows which segments are traversed from start to finish

sample_to_17826286 <- riverdistance(startseg = 1, startvert = 5, endseg = 5, endvert = 14,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17826160 <- riverdistance(startseg = 1, startvert = 5, endseg = 15, endvert = 15,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17827554 <- riverdistance(startseg = 1, startvert = 5, endseg = 16, endvert = 18,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17826190 <- riverdistance(startseg = 1, startvert = 5, endseg = 17, endvert = 25,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17826216 <- riverdistance(startseg = 1, startvert = 5, endseg = 21, endvert = 12,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17826086 <- riverdistance(startseg = 1, startvert = 5, endseg = 27, endvert = 19,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17826188 <- riverdistance(startseg = 1, startvert = 5, endseg = 30, endvert = 57,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17827550 <- riverdistance(startseg = 1, startvert = 5, endseg = 23, endvert = 6,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17827558 <- riverdistance(startseg = 1, startvert = 5, endseg = 28, endvert = 22,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17827556 <- riverdistance(startseg = 1, startvert = 5, endseg = 29, endvert = 20,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17825970 <- riverdistance(startseg = 43, startvert = 23, endseg = 62, endvert = 18,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17826016 <- riverdistance(startseg = 43, startvert = 23, endseg = 66, endvert = 36,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17825884 <- riverdistance(startseg = 43, startvert = 23, endseg = 68, endvert = 25,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17826052 <- riverdistance(startseg = 43, startvert = 23, endseg = 66, endvert = 1,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17825904 <- riverdistance(startseg = 43, startvert = 23, endseg = 73, endvert = 43,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17825952 <- riverdistance(startseg = 43, startvert = 23, endseg = 79, endvert = 13,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17825980 <- riverdistance(startseg = 43, startvert = 23, endseg = 93, endvert = 15,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17826022 <- riverdistance(startseg = 43, startvert = 23, endseg = 92, endvert = 1,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17826090 <- riverdistance(startseg = 43, startvert = 23, endseg = 58, endvert = 32,
                                    rivers = nm_flow_network, map = TRUE) #returns distance in meters
sample_to_17826026 <- riverdistance(startseg = 43, startvert = 23, endseg = 65, endvert = 28,
                                       rivers = nm_flow_network, map = TRUE) #returns distance in meters

nm_distances <- data.frame(site = c('EFJ', 'EFJ','EFJ', 'EFJ', 'EFJ', 'EFJ', 'EFJ', 'EFJ','EFJ','EFJ','SA','SA',
                                    'SA','SA','SA','SA','SA','SA','SA','SA','Redondo'),
                            nhdplus_comid = c('17826286', '17826160', '17827554', '17826190', '17826216','17826086',
                                              '17826188', '17827550', '17827558', '17827556', '17825970', '17826016',
                                              '17825884', '17826052', '17825904', '17825952', '17825980', '17826022',
                                              '17826090', '17826026', NA),
                            event_id = c('NM3581210654120110626', 'NM3581210654120110626', 'NM3581210654120110626',
                                         'NM3581210654120110626', 'NM3581210654120110626', 'NM3581210654120110626',
                                         'NM3581210654120110626', 'NM3589210662020130531','NM3589210662020130531','NM3589210662020130531',
                                         'NM3581210654120110626', 'NM3581210654120110626','NM3581210654120110626','NM3581210654120110626',
                                         'NM3581210654120110626','NM3581210654120110626','NM3581210654120110626','NM3581210654120110626',
                                         'NM3589210662020130531','NM3589210662020130531', 'NM3589210662020130531'),
                            sampling_point_to_stream_m = c(nm_point_to_stream_dist[2,1], nm_point_to_stream_dist[2,1],nm_point_to_stream_dist[2,1],
                                                           nm_point_to_stream_dist[2,1],nm_point_to_stream_dist[2,1],nm_point_to_stream_dist[2,1],
                                                           nm_point_to_stream_dist[2,1],nm_point_to_stream_dist[2,1],nm_point_to_stream_dist[2,1],
                                                           nm_point_to_stream_dist[2,1],nm_point_to_stream_dist[1,1],nm_point_to_stream_dist[1,1],
                                                           nm_point_to_stream_dist[1,1],nm_point_to_stream_dist[1,1],
                                                           nm_point_to_stream_dist[1,1],nm_point_to_stream_dist[1,1],nm_point_to_stream_dist[1,1],
                                                           nm_point_to_stream_dist[1,1],nm_point_to_stream_dist[1,1],nm_point_to_stream_dist[1,1],
                                                           nm_point_to_stream_dist[3,1]),
                            stream_distance_m = c(sample_to_17826286, sample_to_17826160, sample_to_17827554, sample_to_17826190,
                                                  sample_to_17826216, sample_to_17826086, sample_to_17826188, sample_to_17827550,
                                                  sample_to_17827558, sample_to_17827556, sample_to_17825970, 
                                                  sample_to_17826016, sample_to_17825884, sample_to_17826052,
                                                  sample_to_17825904, sample_to_17825952, sample_to_17825980,
                                                  sample_to_17826022, sample_to_17826090, sample_to_17826026, 0),
                            distance_sampling_point_to_fire_km = NA)


nm_distances$distance_sampling_point_to_fire_km <- nm_distances$stream_distance_m/1000

final_distance_df <- rbind(sbc_distances, nm_distances)

write.csv(final_distance_df, "Marss Sampling Point to Fires.csv")

##maps
##https://www.r-bloggers.com/2020/12/visualizing-geospatial-data-in-r-part-2-making-maps-with-ggplot2/
sbc_lter_gav <- ggplot() +
  geom_sf(data = catchment_sbc, fill = 'antiquewhite1', color = 'black') +
  geom_sf(data = catchment_fires_sbc, color = 'black', aes(fill = 'Fire Boundaries')) +
  geom_sf(data = catchment_flowline_sbc, aes(color = 'Streams'), show.legend = 'line') +
  geom_sf(data = sbc_point_to_stream_dist.sf, aes(color = 'Stream Point'), size = 3, show.legend = 'point') +
  geom_sf(data = catchment_sampling_sbc, aes(color = 'Sampling Point'), size = 3, show.legend = 'point') +
  ggtitle("Gaviota") +
  xlab("Longitude") + ylab("Latitude") +
  theme_classic() +
  scale_x_continuous(limits = c(-120.1, -120.3))+
  scale_y_continuous(limits = c(34.47, 34.58)) +
  scale_fill_manual(name = ' ', values = c('Fire Boundaries' = '#b51963'),
                    guide = guide_legend(override.aes = list(linetype = "blank", shape = NA))) +
  scale_color_manual(name = '', values = c('Streams' = '#0073e6', 'Stream Point' = 'green', 'Sampling Point' = 'black', 'Intersection' = 'black'),
                     guide = guide_legend(override.aes = list(linetype = c("blank", "blank", "solid"),
                                                              shape = c(16,16, NA)))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        legend.margin = margin(0, 0, 0,0, unit = 'cm')) +
  annotation_north_arrow(
    location = "tr",
    width = unit(1, 'cm'),
    height = unit(1, 'cm'),
    pad_x = unit(1, 'cm'),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering,
    which_north = 'true')


gav <- sbc_lter_gav + 
  coord_sf(
  xlim = c(-120.228, -120.23),
  ylim= c(34.485, 34.486),
  expand = FALSE) +
  ggtitle('GAV')+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5))

hon <- sbc_lter_gav + 
  coord_sf(
    xlim = c(-120.1420, -120.1410),
    ylim= c(34.475, 34.476),
    expand = FALSE) +
  ggtitle('HON')+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5))

# cowplot::ggdraw(sbc_lter_gav)+
#   cowplot::draw_plot({gav},
#             x=0, y = 0.65,
#             width = 0.5, height = 0.5) +
#   cowplot::draw_plot({hon},
#                      x=0.5, y = 0.75,
#                      width = 0.4, height = 0.2)



sbc_lter_bur <- ggplot() +
  geom_sf(data = catchment_sbc, fill = 'antiquewhite1', color = 'black') +
  geom_sf(data = catchment_fires_sbc, color = 'black', aes(fill = 'Fire Boundaries')) +
  geom_sf(data = catchment_flowline_sbc, aes(color = 'Streams'), show.legend = 'line') +
  geom_sf(data = sbc_point_to_stream_dist.sf, aes(color = 'Stream Point'), size = 3, show.legend = 'point') +
  geom_sf(data = catchment_sampling_sbc, aes(color = 'Sampling Point'), size = 3, show.legend = 'point') +
  ggtitle("sbc_lter_bur") +
  xlab("Longitude") + ylab("Latitude") +
  theme_classic() +
  scale_x_continuous(limits = c(-119.71, -119.78))+
  scale_y_continuous(limits = c(34.4, 34.5)) +
  scale_fill_manual(name = ' ', values = c('Fire Boundaries' = '#b51963'),
                    guide = guide_legend(override.aes = list(linetype = "blank", shape = NA))) +
  scale_color_manual(name = '', values = c('Streams' = '#0073e6', 'Stream Point' = 'green', 'Sampling Point' = 'black', 'Intersection' = 'black'),
                     guide = guide_legend(override.aes = list(linetype = c("blank", "blank", "solid"),
                                                              shape = c(16,16, NA)))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        legend.margin = margin(0, 0, 0,0, unit = 'cm')) +
  annotation_north_arrow(
    location = "tl",
    width = unit(1, 'cm'),
    height = unit(1, 'cm'),
    pad_x = unit(1, 'cm'),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering,
    which_north = 'true')

bur <- sbc_lter_bur + 
  coord_sf(
    xlim = c(-119.742, -119.738),
    ylim= c(34.406, 34.404),
    expand = FALSE) +
  ggtitle('BUR')+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) 
  

sbc_lter_rat <- ggplot() +
  geom_sf(data = catchment_sbc, fill = 'antiquewhite1', color = 'black') +
  geom_sf(data = catchment_fires_sbc, color = 'black', aes(fill = 'Fire Boundaries')) +
  geom_sf(data = catchment_flowline_sbc, aes(color = 'Streams'), show.legend = 'line') +
  geom_sf(data = sbc_point_to_stream_dist.sf, aes(color = 'Stream Point'), size = 3, show.legend = 'point') +
  geom_sf(data = catchment_sampling_sbc, aes(color = 'Sampling Point'), size = 3, show.legend = 'point') +
  ggtitle("Rattlesnake") +
  xlab("Longitude") + ylab("Latitude") +
  theme_classic() +
  scale_x_continuous(limits = c(-119.7, -119.67))+
  scale_y_continuous(limits = c(34.454, 34.50)) +
  scale_fill_manual(name = ' ', values = c('Fire Boundaries' = '#b51963'),
                    guide = guide_legend(override.aes = list(linetype = "blank", shape = NA))) +
  scale_color_manual(name = '', values = c('Streams' = '#0073e6', 'Stream Point' = 'green', 'Sampling Point' = 'black', 'Intersection' = 'black'),
                     guide = guide_legend(override.aes = list(linetype = c("blank", "blank", "solid"),
                                                              shape = c(16,16, NA)))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        legend.margin = margin(0, 0, 0,0, unit = 'cm')) +
  annotation_north_arrow(
    location = "tl",
    width = unit(1, 'cm'),
    height = unit(1, 'cm'),
    pad_x = unit(1, 'cm'),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering,
    which_north = 'true')

rat <- sbc_lter_rat + 
  coord_sf(
    xlim = c(-119.695, -119.690),
    ylim= c(34.46, 34.456),
    expand = FALSE) +
  ggtitle('Rattlesnake')+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5))
  



nm_lter <- ggplot() +
  geom_sf(data = catchment_nm, fill = 'antiquewhite1', color = 'black') +
  geom_sf(data = catchment_fires_nm, color = 'black', aes(fill = 'Fire Boundaries')) +
  geom_sf(data = catchment_flowline_nm, aes(color = 'Streams'), show.legend = 'line') +
  geom_sf(data = nm_point_to_stream_dist.sf, aes(color = 'Stream Point'), size = 3, show.legend = 'point') +
  geom_sf(data = catchment_sampling_nm, aes(color = 'Sampling Point'), size = 3, show.legend = 'point') +
  ggtitle("New Mexico Sites") +
  xlab("Longitude") + ylab("Latitude") +
  theme_classic() +
  scale_x_continuous(limits = c(-106.4, -106.63))+
  scale_y_continuous(limits = c(35.83, 36.03)) +
  scale_fill_manual(name = ' ', values = c('Fire Boundaries' = '#b51963'),
                    guide = guide_legend(override.aes = list(linetype = "blank", shape = NA))) +
  scale_color_manual(name = '', values = c('Streams' = '#0073e6', 'Stream Point' = 'green', 'Sampling Point' = 'black', 'Intersection' = 'black'),
                     guide = guide_legend(override.aes = list(linetype = c("blank", "blank", "solid"),
                                                              shape = c(16,16, NA)))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        legend.margin = margin(0, 0, 0,0, unit = 'cm')) +
  annotation_north_arrow(
    location = "tl",
    width = unit(1, 'cm'),
    height = unit(1, 'cm'),
    pad_x = unit(1, 'cm'),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering,
    which_north = 'true')
  
  RED <- nm_lter + 
    coord_sf(
      xlim = c(-106.59, -106.605),
      ylim= c(35.86, 35.870),
      expand = FALSE) +
    ggtitle('Redondo')+
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5))
  
  EFJ <- nm_lter + 
    coord_sf(
      xlim = c(-106.4899, -106.491),
      ylim= c(35.848, 35.849),
      expand = FALSE) +
    ggtitle('Jemez')+
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5))
  
  SA <- nm_lter + 
    coord_sf(
      xlim = c(-106.595, -106.60),
      ylim= c(35.972, 35.975),
      expand = FALSE) +
    ggtitle('San Antonio')+
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5))




