# Create a spatial object representing the area of interest, using a self-defined
# range of coordinates (format: longitude, latitude).


# general libraries
library(tidyverse)
# GIS libraries
library(sp)
library(rgdal)

# save the polygon coordinates in a matrix (format: longitude, latitude)
coords <- matrix(c(38.56560025749677, -3.400672109607846,
                   38.140846474271804, -3.567770809227452,
                   38.143517881682136, -3.743724862559294,
                   38.466570524622405, -3.8142801796984336,
                   38.66581180152877, -4.150945007253591,
                   38.96060890710317, -3.965217031121372, 
                   39.05079692594667, -3.7277823254453137
                  ), ncol = 2, byrow = TRUE)

# create the polygon and reformat into SpatialPolygonsDataFrame (for saving as shapefile)
area <- Polygon(coords) %>% 
  list() %>% 
  Polygons(ID = 'ww_kenya') %>% 
  list() %>% 
  SpatialPolygons(proj4string = 
                    CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>% 
  as('SpatialPolygonsDataFrame')

# check the polygon shape
plot(area, axes=TRUE)

# save the polygon as a shapefile & kml
writeOGR(area, 'data', 'study_area', driver = 'ESRI Shapefile', overwrite_layer =TRUE)
writeOGR(area, dsn = 'data/study_area.kml', layer = 'study_area', driver = 'KML', 
         overwrite_layer = TRUE) # can be visualized in Sentinelhub
