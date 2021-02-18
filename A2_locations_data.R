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
# library(mopa) # library for pseudo-absence sampling but has been removed from CRAN
# (could get an old version)


# # set a name for the output
# name_out <- 'tr_Jan10-Jul20_ts_Aug20-Dec20'
# set coordinate reference system
crs_kenya <- crs("+proj=utm +zone=37 +south +datum=WGS84 +units=m +no_defs")

# retrieve the area of interest polygon
aoi <- readOGR(dsn = "data", layer = "study_area") %>% 
  spTransform(crs_kenya)

# any output objects are stored in an output/ folder
if (! dir.exists('output')){
  dir.create('output')
}


#------------------------------------------------------------
# login to collect data
#------------------------------------------------------------
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


#------------------------------------------------------------
# incident data collection
#------------------------------------------------------------

i.type <- '["https://sensingclues.poolparty.biz/SCCSSOntology/97"]'
bounds <- '{"south":-5.104080690471747,"west":37.03937894518083,
            "north":-2.56592681976295,"east":40.349053193120966}'
dt.range <- '{"to":"2020-07-31T22:00:00.000Z","from":"2010-01-04T23:00:00.000Z"}'

get_Focus_data <- function(incident.type, map.bounds, date.time.range){
  q <- paste0('{"filters":{"geoQuery":{"operator":"intersects","mapBounds":', map.bounds, ',"drawings":[]},
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
  DATA <- as.data.frame(DATA)
  row.names(DATA) <- NULL
  names(DATA) <- c("lon","lat")
  coordinates(DATA) <- ~lon+lat
  proj4string(DATA) <- "+proj=longlat +ellps=WGS84 +no_defs"
  
  return(DATA)
}


train <- get_Focus_data(i.type, bounds, dt.range)
test <- get_Focus_data(i.type, bounds, date.time.range = '{"to":"2021-01-31T22:00:00.000Z","from":"2020-08-01T23:00:00.000Z"}')

# transform to projected coordinate system (in meters)
train.T <- spTransform(train, crs_kenya)
test.T <- spTransform(test, crs_kenya)

# remove points outside of aoi
train.T <- train.T[aoi]
test.T <- test.T[aoi]

# quick display of incident locations
map <- leaflet() %>% 
  addTiles() %>% 
  addPolygons(data = spTransform(aoi, crs('+proj=longlat +ellps=WGS84 +no_defs'))) %>% 
  addCircleMarkers(data = train, color = 'blue') %>% 
  addCircleMarkers(data = test, color = 'purple')
map


#------------------------------------------------------------
# create pseudo-absence locations randomly using raster cell numbers
#------------------------------------------------------------

# create the desired output raster from the aoi Polygon object
out <- raster(aoi, resolution = 1000) # 1 x 1 km spatial resolution

# there are duplicates in the incident data:
length(train.T)
length(remove.duplicates(train.T))

get_pseudo_absence_locations <- function(area_sp, area_grid, presence_data_sp, nabsent = NULL, seed = 555){
  
  # check crs of all objects are the same
  # if not, the sp objects will be transformed to have the same crs as the grid object
  if (crs(area_grid)@projargs != crs(area_sp)@projargs | 
      crs(area_grid)@projargs != crs(presence_data_sp)@projargs){
    area_sp <- spTransform(area_sp, crs(area_grid))
    presence_data_sp <- spTransform(presence_data_sp, crs(area_grid))
  }
  
  # get all cell numbers in AOI
  cell_aoi <- cellFromPolygon(area_grid, area_sp)[[1]]
  
  # get cell numbers of the indicent locations
  cell_pr <- cellFromXY(area_grid, presence_data_sp) %>% 
    unique() # remove duplicate cells (incidents in same raster cell)
  
  # randomly sample pseudo-absence cell numbers from the non-incident cells
  cell_opt <- setdiff(cell_aoi, cell_pr)

  set.seed(seed)
  if (is.null(nabsent)){ # get equal amount of absent as present
    cell_ab <- sample(cell_opt, length(cell_pr))
  } else {
    cell_ab <- sample(cell_opt, nabsent)
  }

  return(list(cell_ab, cell_pr, cell_aoi))
}


# assign each location in aoi to absent/present/other for both training and testing
tr_cell_list <- get_pseudo_absence_locations(aoi, out, train.T)
tr_cell_ab <- tr_cell_list[[1]]
tr_cell_pr <- tr_cell_list[[2]]
tr_cell_aoi <- tr_cell_list[[3]]

ts_cell_list <- get_pseudo_absence_locations(aoi, out, test.T, seed = 500)
ts_cell_ab <- ts_cell_list[[1]]
ts_cell_pr <- ts_cell_list[[2]]
ts_cell_aoi <- ts_cell_list[[3]]


# check that absent and present cells do not overlap
sum(tr_cell_pr %in% tr_cell_ab)
sum(ts_cell_pr %in% ts_cell_ab)

# percentage of locations in AOI used in training the model
(length(tr_cell_pr) + length(tr_cell_ab)) / length(tr_cell_aoi) * 100


# add location labels to output raster
out$train <- NA
out$train[tr_cell_aoi] <- rep(2, length(tr_cell_aoi))
out$train[tr_cell_pr] <- rep(1, length(tr_cell_pr))
out$train[tr_cell_ab] <- rep(0, length(tr_cell_ab))

out$test <- NA
out$test[ts_cell_aoi] <- rep(2, length(ts_cell_aoi))
out$test[ts_cell_pr] <- rep(1, length(ts_cell_pr)) # NA in presents because outside aoi!
out$test[ts_cell_ab] <- rep(0, length(ts_cell_ab))


# save the output raster
writeRaster(out, 'output/raster_template', overwrite = TRUE)


# visualize training incident and pseudo data on map
plot(out)

plot(out$train, legend = FALSE, col = c('firebrick3', 'seagreen3', 'lightgrey'), 
     main = 'Training data - randomly sampled pseudo-absence')
legend("bottomleft", legend = c("pseudo-absent", "incident", "other"), 
       fill = c('firebrick3', 'seagreen3', 'lightgrey'), title = 'Assigned set')

map <- leaflet() %>% 
  addTiles() %>% 
  addRasterImage(ratify(out$train))
map


#------------------------------------------------------------
# create final DataFrame
#------------------------------------------------------------

# create DataFrame to store features in with the cell numbers as identifiers of the location
# and a presence/absence column with 1/0 values and a label column incident/pseudo-absence/other

df <- data.frame(cell_id = tr_cell_aoi, train = out$train[tr_cell_aoi], test = out$test[tr_cell_aoi]) %>% 
  mutate(train_label = ifelse(train == 0, 'pseudo-absent', ifelse(train == 1, 'incident', 'other')),
         test_label = ifelse(test == 0, 'pseudo-absent', ifelse(test == 1, 'incident', 'other')))

# check if labels for cell numbers are correct
head(df %>% filter(train == 1), 10)$cell_id == sort(tr_cell_pr)[1:10]


# store DataFrame as csv
write.csv(x = df, file = 'output/location_features.csv', row.names = FALSE)
# write.csv(x = df, file = paste0('output/location_features_', name_out, '.csv'), row.names = FALSE)



# Comments:
# xyFromCell() to get center of raster cells from list of cell numbers (incidents)
# possibly dont need xyFromCell() but could just use cell numbers?
# raster::extract() can take: a numeric vector representing cell numbers
# so only need to reproject DEM/NDVI rasters to output raster dimensions/resolution
