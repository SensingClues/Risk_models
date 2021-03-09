# Retrieve Sentinel-2 imagery as raster objects.
# Derive indices (e.g. NDVI) from the bands and reproject them 
# to the desired output raster crs and dimensions.
# Extract the index time series of 2020 at the locations of interest
# and get index maximum and mean values per location.


# GIS libraries
library(raster)
library(rgdal)
library(sp)
library(leaflet)
# other
library(tidyverse)
library(ggplot2)
library(Rfast)
library(matrixStats)


#------------------------------------------------------------
# setup
#------------------------------------------------------------

# choose a scenario name
sc_name <- 'SC_6month' # other scenarios: 'SC_3month', 'SC_1month', 'SC_alldata'

df <- read.csv(paste0('output/location_features_', sc_name, '.csv'))

aoi <- readOGR(dsn = "data", layer = "study_area")

out <- raster(paste0('output/rasters_', sc_name, '.tif'))


#------------------------------------------------------------
# create indices from Sentinel-2 imagery (from EO browser)
#------------------------------------------------------------
folder <- 'data/Sentinel_2/'

# function to retrieve Sentinel-2 images and create an index (e.g. NDVI)
sentinel2_index <- function(name1, name2, sent2_path, cloudmask = TRUE){
  # IMPORTANT: index is executed: (name1 - name2) / (name1 + name2)
  
  # load Sentinel-2 images (retrieved from EO browser)
  stack1 <- stack(paste0(folder, list.files(folder, name1)))
  stack2 <- stack(paste0(folder, list.files(folder, name2)))
  
  # check equal nlayers()
  if (nlayers(stack1) != nlayers(stack2)){
    warning('Sentinel-2 images are missing, check the designated folder.')
  }
  
  if (cloudmask){
    classified_stack <- stack(paste0(folder, list.files(folder, 'classification')), bands = 1)
    if (nlayers(stack1) != nlayers(classified_stack)){
      warning('Sentinel-2 scene classification images are missing, check the designated folder.')
    }
  }
  
  index_stack <- stack()
  date_names <- NULL
  for (i in 1:length(names(stack1))) {
    date1 <- str_sub(names(stack1[[i]]), start = 2, end = 11)
    date2 <- str_sub(names(stack2[[i]]), start = 2, end = 11)
    if (date1 == date2) {
      print(paste0(i, ' out of ', length(names(stack1))))
      date_names <- c(date_names, date1)
      layer <- (stack1[[i]] - stack2[[i]])/(stack1[[i]] + stack2[[i]])
      
      if (cloudmask){
        # using sentinel classification to replace cloud pixels with NA value for NDVI
        # value 65535 = cloud high probability
        # value 49345 = cloud medium probability
        cl <- classified_stack[[i]]
        cl[cl == 65535] <- NA
        cl[cl == 49345] <- NA
        layer <- mask(layer, cl)
      }
      index_stack <- addLayer(index_stack, layer)
    }
  }
  names(index_stack) <- date_names
  return(index_stack)
}


# vegetation index (NDVI)
ndvi_stack <- sentinel2_index('L2A_B08', 'L2A_B04', sent2_path = folder)
plot(ndvi_stack[[1]])

# moisture index (NDMI)
ndmi_stack <- sentinel2_index('L2A_B8A', 'L2A_B11', sent2_path = folder)
plot(ndmi_stack[[1]])


# reproject to res/crs of output raster
ndvi_stack <- projectRaster(ndvi_stack, to = out, method = 'bilinear')
ndmi_stack <- projectRaster(ndmi_stack, to = out, method = 'bilinear')


#------------------------------------------------------------
# extract a time series of each index at incident locations
#------------------------------------------------------------

# NDVI

ndvi_vals <- raster::extract(ndvi_stack,  df$cell_id)
# View(ndvi_vals)

# plot an ndvi time series for a single location
dates <- as.Date(gsub('X', '', colnames(ndvi_vals)), '%Y.%m.%d')
plot(dates, ndvi_vals[1,], type = 'l')
points(dates, ndvi_vals[1,])

# get ndvi max and mean
ndvi_m <- rowMeans(ndvi_vals, na.rm = TRUE)
ndvi_mx <- rowMaxs(ndvi_vals, na.rm = TRUE) # max because clouds give low NDVI


# NDMI

ndmi_vals <- raster::extract(ndmi_stack,  df$cell_id)

# get ndmi min and mean
ndmi_m <- rowMeans(ndmi_vals, na.rm = TRUE)
ndmi_mn <- rowMins(ndmi_vals, na.rm = TRUE) # min because clouds give high NDMI



features <- cbind(df, ndvi_mean = ndvi_m) %>% 
  cbind(ndvi_max = ndvi_mx) %>% 
  cbind(ndmi_mean = ndmi_m) %>% 
  cbind(ndmi_min = ndmi_mn)
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
write.csv(features, paste0('output/location_features_', sc_name, '.csv'), row.names = FALSE)





#------------------------------------------------------------
# when using copernicus hub jp2 files
#------------------------------------------------------------
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
