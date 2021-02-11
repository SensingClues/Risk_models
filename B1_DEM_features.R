# Retrieve the DEM rasters for the area of interest in Kenya.
# Reproject it to the desired output raster crs and dimensions.
# Extract the DEM values at the locations of interest.


# GIS libraries
library(raster)
library(rgdal)
library(sp)
library(leaflet)
# other
library(tidyverse)
library(ggplot2)


# retrieve locations in area of interest (AOI)
df <- read.csv('output/location_features.csv')


#---------------------------------
# load, reproject and mask DEM
#---------------------------------
# area of interest
aoi <- readOGR(dsn = "data", layer = "study_area")

out <- raster('output/raster_template.grd')

# merge DEM tiles, mask them for the aoi and reproject them to desired resolution/crs (out)
dem_aoi <- raster('data/s04_e038_1arc_v3.tif') %>% 
  merge(raster('data/s05_e038_1arc_v3.tif')) %>% 
  merge(raster('data/s04_e039_1arc_v3.tif')) %>% 
  crop(aoi) %>% 
  mask(aoi) %>% 
  # a 7x7 matrix is the smallest focal window resulting in no NAs in 1000m DEM
  focal(w=matrix(1,7,7), fun=mean, NAonly=TRUE, na.rm=TRUE) %>% 
  projectRaster(to = out, method = 'bilinear')
plot(dem_aoi)

# looking at a smaller area within the aoi:
# plot(crop(dem_aoi, extent(dem_aoi, 350, 550, 1500, 1700)))
# dem_small <- crop(dem_aoi, extent(dem_aoi, 2800, 3200, 2300, 2800))


# if data from sentinelhub: (download failed, the files were empty)
# dem_aoi <- raster(paste0('data/', list.files('data', 'DEM')))


#---------------------------------
# extract DEM values at incident 
# locations (with cell numbers)
#---------------------------------

dem_vals <- raster::extract(dem_aoi, df$cell_id)

features <- cbind(df, elevation_m = dem_vals)


# compare histograms
par(mfrow = c(2,2))
  hist(features$elevation_m[features$incident == 1], breaks = seq(100, 1600, 30), main = 'Incident locations')
  hist(features$elevation_m[features$incident == 0], breaks = seq(100, 1600, 30), main = 'Pseudo-absence locations')
  # compare to DEM values in complete aoi
  hist(values(dem_aoi), breaks = seq(100, 1600, 30), main = 'Complete area of interest')
  hist(features$elevation_m, breaks = seq(100, 1600, 30))
par(mfrow = c(1,1))

# with ggplot instead
# ggplot(features, aes(elevation_m)) +
#   geom_histogram(binwidth = 25, na.rm = TRUE)


#---------------------------------
# other DEM related features at 
# incident locations
#---------------------------------

# slope in degrees
slp_aoi <- terrain(dem_aoi, opt = 'slope', unit = 'degrees')
plot(slp_aoi)

slp_vals <- raster::extract(slp_aoi, df$cell_id)


# aspect (values from 1 to 6, wind directions)
asp_aoi <- terrain(dem_aoi, opt = 'aspect')
plot(asp_aoi)

asp_vals <- raster::extract(asp_aoi, df$cell_id)


features <- cbind(features, slope_deg = slp_vals) %>% 
  cbind(aspect = asp_vals)



# compare histograms for slope
par(mfrow = c(2,2))
hist(features$slope_deg[features$incident == 1], breaks = seq(0, 30, 0.25), main = 'Incident locations')
hist(features$slope_deg[features$incident == 0], breaks = seq(0, 30, 0.25), main = 'Pseudo-absence locations')
# compare to DEM values in complete aoi
hist(values(slp_aoi), breaks = seq(0, 30, 0.25), main = 'Complete area of interest')
par(mfrow = c(1,1))

# compare histograms for aspect
par(mfrow = c(2,2))
hist(features$aspect[features$incident == 1], breaks = seq(0, 6.5, 0.25), main = 'Incident locations')
hist(features$aspect[features$incident == 0], breaks = seq(0, 6.5, 0.25), main = 'Pseudo-absence locations')
# compare to DEM values in complete aoi
hist(values(asp_aoi), breaks = seq(0, 6.5, 0.25), main = 'Complete area of interest')
par(mfrow = c(1,1))



# visualize the DEM raster on a map
map <- leaflet() %>% 
  addTiles() %>% 
  addRasterImage(projectRaster(dem_aoi, crs = crs('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))) %>% 
  addPolygons(data = aoi)
map


# store the new features with the existing location_features file
write.csv(features, 'output/location_features.csv', row.names = FALSE)
