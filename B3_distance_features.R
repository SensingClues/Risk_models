# Retrieve the available distance rasters (towns, roads, water bodies) 
# and reproject them to the desired output raster crs and dimensions.
# Extract the distance values at the locations of interest.


# GIS libraries
library(raster)


#------------------------------------------------------------
# setup
#------------------------------------------------------------

# choose a scenario name
sc_name <- 'SC_6month' # other scenarios: 'SC_3month', 'SC_1month', 'SC_alldata'

df <- read.csv(paste0('output/location_features_', sc_name, '.csv'))

out <- raster(paste0('output/rasters_', sc_name, '.tif'))


#------------------------------------------------------------
# extract nearest distance at locations from distance rasters
#------------------------------------------------------------

# distance to towns
# retrieve distance raster and reproject to resolution, crs and extent of output raster
dtowns <- raster('data/Distance_to_Town.grd') %>% 
  projectRaster(to = out, method = 'bilinear')
plot(dtowns)
# extract the distance values at the area of interest location grid cells
dtown_vals <- raster::extract(dtowns, df$cell_id)


# distance to water bodies
dwater <- raster('data/Distance_to_River.grd') %>% 
  projectRaster(to = out, method = 'bilinear') 
plot(dwater)
dwater_vals <- raster::extract(dwater, df$cell_id)


# distance to roads
droads <- raster('data/Distance_to_Road.grd') %>% 
  projectRaster(to = out, method = 'bilinear')
plot(droads)
droad_vals <- raster::extract(droads, df$cell_id)


# add all distance features to locations DataFrame
features <- cbind(df, dtown = dtown_vals) %>% 
  cbind(dwater = dwater_vals) %>% 
  cbind(droad = droad_vals)
summary(features)


# store the new features with the existing location_features file
write.csv(features, paste0('output/location_features_', sc_name, '.csv'), row.names = FALSE)



