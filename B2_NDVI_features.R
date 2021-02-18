# Retrieve Sentinel-2 band 4 (red) and band 8 (NIR) rasters.
# Derive NDVI rasters from the bands and reproject them to the 
# desired output raster crs and dimensions.
# Extract the NDVI time series of 2020 at the locations of interest
# and get NDVI maximum and mean values per location.


# GIS libraries
library(raster)
library(rgdal)
library(sp)
library(leaflet)
# other
library(tidyverse)
library(ggplot2)
library(Rfast)


df <- read.csv('output/location_features.csv')

out <- raster('output/raster_template.grd')
plot(out)

#---------------------------------
# create NDVI time series RasterStack 
# from Sentinel-2 band 4(red) and band 8(NIR)
#---------------------------------
folder <- 'data/Sentinel_2/'
aoi <- readOGR(dsn = "data", layer = "study_area")

# load Sentinel-2 images (retrieved from EO browser)
b04_stack <- stack(paste0(folder, list.files(folder, 'L2A_B04')))
b08_stack <- stack(paste0(folder, list.files(folder, 'L2A_B08')))
classified_stack <- stack(paste0(folder, list.files(folder, 'classification')), bands = 1)

# check equal nlayers()
nlayers(b04_stack) == nlayers(b08_stack) & nlayers(b04_stack) == nlayers(classified_stack)

ndvi_stack <- stack()
ndvi_names <- NULL
for (i in 1:length(names(b04_stack))) {
  date_b04 <- str_sub(names(b04_stack[[i]]), start = 2, end = 11)
  date_b08 <- str_sub(names(b08_stack[[i]]), start = 2, end = 11)
  if (date_b04 == date_b08) {
    print(paste0(i, ' out of ', length(names(b04_stack))))
    ndvi_names <- c(ndvi_names, date_b04)
    new_ndvi <- (b08_stack[[i]] - b04_stack[[i]])/(b08_stack[[i]] + b04_stack[[i]])
    
    # using sentinel classification to replace cloud pixels with NA value for NDVI
    # value 65535 = cloud high probability
    # value 49345 = cloud medium probability
    cl <- classified_stack[[i]]
    cl[cl == 65535] <- NA
    cl[cl == 49345] <- NA
    new_ndvi <- mask(new_ndvi, cl)
    
    ndvi_stack <- addLayer(ndvi_stack, new_ndvi)
  }
}
names(ndvi_stack) <- ndvi_names

plot(ndvi_stack[[1]])

# remove large objects
remove(b04_stack, b08_stack, classified_stack)


# reproject to res/crs of output raster
ndvi_stack <- projectRaster(ndvi_stack, to = out, method = 'bilinear')
plot(ndvi_stack[[1]])


#---------------------------------
# extract NDVI time series at incident locations
#---------------------------------

ndvi_vals <- raster::extract(ndvi_stack,  df$cell_id)
View(ndvi_vals)

# plot an NDVI time series for a single location
dates <- as.Date(gsub('X', '', colnames(ndvi_vals)), '%Y.%m.%d')

plot(dates, ndvi_vals[1,], type = 'l')
points(dates, ndvi_vals[1,])


# get ndvi max and mean
ndvi_m <- rowMeans(ndvi_vals, na.rm = TRUE)
ndvi_mx <- apply(ndvi_vals, 1, max)

features <- cbind(df, ndvi_mean = ndvi_m) %>% 
  cbind(ndvi_max = ndvi_mx)
summary(features)


# using savitzky-golay filter (concerning cloud cover outliers, however only works if first value is not affected by clouds)
# library(signal)
# sg <- sgolayfilt(ndvi_vals[1,])
# plot(sg, type="l")
# lines(filtfilt(rep(1, 5)/5,1,ndvi_vals[1,]), col = "red") # averaging filter
# points(ndvi_vals[1,], pch = "x")



# run incident_data_Focus.R first to get DATA
# check ndvi at incidents 
map <- leaflet() %>% 
  addTiles() %>% 
  addRasterImage(crop(ndvi_stack[[1]], DATA[1:10])) %>% 
  addAwesomeMarkers(data=DATA[1:10])
map


# store the new features with the existing location_features file
write.csv(features, 'output/location_features.csv', row.names = FALSE)







#---------------------------------
# when using copernicus hub jp2 files
#---------------------------------
# more complicated file structure and have to convert .jp2 format to .tif

aoi <- readOGR(dsn = "data", layer = "study_area") %>% 
  spTransform(crs("+proj=utm +zone=37 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

library(gdalUtils)

gdalUtils::gdal_chooseInstallation('JP2OpenJPEG')

in.files <- list.files(folder,"jp2$")
out.files <- gsub(".jp2", ".tif", in.files)
for(i in 1:length(in.files)){
  raster::writeRaster(raster(rgdal::readGDAL(paste0(folder, in.files[i]))),
                      file.path(str_sub(folder, end = 15), out.files[i]),
                      overwrite = TRUE)
}
raster::writeRaster(raster(rgdal::readGDAL(paste0(folder, 'MSK_CLDPRB_20m.jp2'))),
                    file.path(paste0(folder, 'MSK_CLDPRB_20m.tif')),
                    overwrite = TRUE)

# for(i in 1:length(in.files)){
#   gdalUtils::gdal_translate(in.files[i], out.files[i])
# }


r20200629 <- stack(out.files[1], out.files[2]) %>% 
  crop(aoi) %>% 
  mask(aoi)
names(r20200629) <- c('B04', 'B08')

ndvi_cop <- (r20200629[[2]] - r20200629[[1]])/(r20200629[[2]] + r20200629[[1]])

plot(ndvi_cop)

cloud <- raster(paste0(folder, 'MSK_CLDPRB_20m.tif')) %>% 
  crop(aoi) %>% 
  mask(aoi)
plot(cloud)
