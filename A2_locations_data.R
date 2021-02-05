# Retrieve incident locations from the Focus API (Sening Clues).
# Create pseudo-absence locations sampled from the area of interest.
# Save all locations (identified by raster cell number) in a DataFrame with labels:
# incident, pseudo-absent, other


# general libaries
library(tidyverse)
# focus sparql endpoint
library(httr)
library(jsonlite)
library(SPARQL)
# GIS libraries
library(raster)
library(rgdal)
library(sp)
library(leaflet)
# library(mopa) # library for pseudo-absence sampling but has been removed from CRAN, could get an old version


# set a name for the output
name_out <- 'train'

# any output objects are stored in an output/ folder
if (! dir.exists('output')){
  dir.create('output')
}


#---------------------------------
# login to collect data
#---------------------------------
library(getPass)
URL <- "https://focus.sensingclues.org/"
url.login <- paste0(URL,"api/auth/login")
# we set up an authenticated session
rl <- POST(url.login, 
           body = jsonlite::toJSON(list(username = getPass("Enter your Cluey username:"),
                                        password=getPass()), 
                                   auto_unbox = TRUE), 
           encode = "raw",
           content_type_json())

#---------------------------------
# incident data collection
#---------------------------------

incident.type <- '["https://sensingclues.poolparty.biz/SCCSSOntology/97"]'
map.bounds <- '{"south":-5.104080690471747,"west":37.03937894518083,
"north":-2.56592681976295,"east":40.349053193120966}'
date.time.range <- '{"to":"2020-07-31T22:00:00.000Z","from":"2010-01-04T23:00:00.000Z"}'

q <- paste0('{"filters":
{"geoQuery":{"operator":"intersects","mapBounds":', map.bounds, ',"drawings":[]},
"concepts":', incident.type,',"dataSources":[],"dateTimeRange":', date.time.range,'},
"options":{"start":0,"pageLength":20},"start":1,"pageLength":500}')
url.collect <- paste0(URL, 'api/map/all/default/0/features')
incidentDATA.L <- POST(url.collect, body=q, encode="raw", content_type_json())
incidentDATA <- content(incidentDATA.L)

Nobs <- length(incidentDATA$features)
if (Nobs <= 0){
  warning('No incident data has been collected.')
} else {
  print(paste0('There have been ', Nobs, ' incident locations collected in this session.'))
} 
  
DATA <- NULL
for(i in 1:Nobs) {
  coords <- unlist(incidentDATA$features[[i]]$geometry$coordinates)
  # Only read JSON files
  DATA <- rbind(DATA,coords)
}
DATA <- as.data.frame(DATA, row.names = NULL)
names(DATA) <- c("lon","lat")
row.names(DATA) <- NULL
coordinates(DATA) <- ~lon+lat
proj4string(DATA) <- "+proj=longlat +ellps=WGS84 +no_defs"
DATA.T <- spTransform(DATA, crs("+proj=utm +zone=37 +south +datum=WGS84 +units=m +no_defs"))

# quick display of incident locations
map <- leaflet() %>% 
  addTiles() %>% 
  addAwesomeMarkers(data=DATA)
map


#---------------------------------
# create pseudo-absence locations 
# randomly using raster cell numbers
#---------------------------------

# create the desired output raster from the aoi Polygon object
aoi <- readOGR(dsn = "data", layer = "study_area") %>% 
  spTransform(crs("+proj=utm +zone=37 +south +datum=WGS84 +units=m +no_defs"))

out <- raster(aoi, resolution = 1000) %>% # 1 x 1 km spatial resolution
  mask(aoi)


# get all cell numbers in AOI
crs(out)@projargs == crs(aoi)@projargs # check crs
cell_aoi <- cellFromPolygon(out, aoi)[[1]]
length(cell_aoi)


# there are duplicates in the incident data:
length(DATA.T)
length(remove.duplicates(DATA.T))

# get cell numbers of the indicent locations
crs(out)@projargs == crs(DATA.T)@projargs # check crs
cell_pr <- cellFromXY(out, DATA.T) %>% 
  unique() # remove duplicate cells (incidents in same raster cell)
length(cell_pr)


# randomly sample pseudo-absence cell numbers from the non-incident cells
cell_opt <- setdiff(cell_aoi, cell_pr)
length(cell_opt)

set.seed(555)
cell_ab <- sample(cell_opt, length(cell_pr))
length(cell_ab)

sum(cell_pr %in% cell_ab)


# percentage of locations in AOI used in training the model
(length(cell_pr) + length(cell_ab)) / length(cell_aoi) * 100


# add location labels to output raster
out[cell_aoi] <- rep(2, length(cell_aoi))
out[cell_pr] <- rep(1, length(cell_pr))
out[cell_ab] <- rep(0, length(cell_ab))

# save the output raster
writeRaster(out, 'output/raster_template', overwrite = TRUE)

# visualize training incident and pseudo data on map
plot(out, legend = FALSE, col = rev(terrain.colors(3)), main = 'Training data - randomly sampled pseudo-absence')
legend("bottomleft", legend = c("pseudo-absent", "incident", "other"), 
       fill = rev(terrain.colors(3)), title = 'Assigned set')

map <- leaflet() %>% 
  addTiles() %>% 
  addRasterImage(ratify(out))
map


#---------------------------------
# create final DataFrame
#---------------------------------

# create DataFrame to store features in with the cell numbers as identifiers of the location
# and a presence/absence column with 1/0 values and a label column incident/pseudo-absence/other

df <- data.frame(cell_id = cell_aoi, incident = out[cell_aoi]) %>% 
  mutate(cell_label = ifelse(incident == 0, 'pseudo-absent', ifelse(incident == 1, 'incident', 'other')))

# check if labels for cell numbers are correct
head(df %>% filter(incident == 1), 10)$cell_id == sort(cell_pr)[1:10]


# store DataFrame as csv
write.csv(x=df, file= paste0('output/location_features_', name_out, '.csv'), row.names = FALSE)



# Comments:
# xyFromCell() to get center of raster cells from list of cell numbers (incidents)
# possibly dont need xyFromCell() but could just use cell numbers?
# raster::extract() can take: a numeric vector representing cell numbers
# so only need to reproject DEM/NDVI rasters to output raster dimensions/resolution
