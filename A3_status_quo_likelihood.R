# Create the status quo for predicting incident likelihood:
# a hotspot map based on incident density corrected for patrolling density.
# Code by Koen de Koning.


# focus sparql endpoint
library(httr)
library(jsonlite)
library(SPARQL)
# GIS libraries
library(raster)
library(rgdal)
library(sp)
library("KernSmooth")
library(leaflet)


#------------------------------------------------------------
# setup
#------------------------------------------------------------

# set date range of training set
date_range <- c('2020-02-01', '2020-09-01') # (from, to)

# choose a scenario name
sc_name <- 'SC_6month' # other scenarios: 'SC_3month', 'SC_1month', 'SC_alldata'

# mapping parameters:
HotspotColor = c("#040404B3","#080918B3","#0E0D24B3","#150F2EB3","#1D1135B3","#24123CB3","#2C1242B3","#341348B3","#3C134EB3","#451353B3","#4D1259B3",
                 "#56125DB3","#5F1162B3","#681066B3","#701069B3","#79106DB3","#82106FB3","#8A1172B3","#931373B3","#9B1674B3","#A31A75B3","#AB1E75B3",
                 "#B32375B3","#BA2973B3","#C12F71B3","#C8356FB3","#CF3B6BB3","#D64267B3","#DC4962B3","#E2505BB3","#E85752B3","#ED5F48B3","#F2673AB3",
                 "#F37133B3","#F47B2CB3","#F58426B3","#F58E23B3","#F69622B3","#F79F25B3","#F7A82CB3","#F7B134B3","#F8B93EB3","#F8C149B3","#F8CA54B3",
                 "#F9D25FB3","#F9DB6BB3","#FAE377B3","#FBEC84B3","#FDF490B3","#FFFE9EB3")


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

get_Focus_data <- function(incident.type, map.bounds, date.range){
  q <- paste0('{"filters":{"geoQuery":{"operator":"intersects","mapBounds":', map.bounds, ',"drawings":[]},
              "concepts":', incident.type,',"dataSources":[],
              "dateTimeRange":{"to":"',date.range[2],'T22:00:00.000Z","from":"',date.range[1],'T23:00:00.000Z"}},
              "options":{"start":0,"pageLength":20},"start":1,"pageLength":1000}')
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


DATA <- get_Focus_data(i.type, bounds, date_range)

map <- leaflet() %>% 
  addTiles() %>% 
  addAwesomeMarkers(data = DATA)
map


#------------------------------------------------------------
# tracks data collection
#------------------------------------------------------------

# tracks before 2020 are stored locally, tracks from 2020 onwards are in Focus
if (as.Date(date_range[1]) < as.Date("2020-01-01")){
  print('Retrieving tracks from locally stored files.')
  
  TRACKS <- read.csv('data/processed/tracks_2015_to_2019.csv')
  
} else { # if From Date >= 2020-01-01
  print('Retrieving tracks from Focus.')
  
  url.tracks <- paste0(URL,"api/map/all/track/0/features")
  q <- paste0('{"filters":{"geoQuery":{"operator":"intersects","mapBounds":', bounds, '"drawings":[]},
              "dateTimeRange":{"to":"',date_range[2],'T22:00:00.000Z","from":"',date_range[1],'T23:00:00.000Z"}}
              ,"options":{"start":0,"pageLength":20},"start":1,"pageLength":1000}')
  trackDATA.L <- POST(url.tracks, body=q, encode="raw", content_type_json())
  trackDATA <- content(trackDATA.L)
  Ntracks <- length(trackDATA$features)
  
  TRACKS <- NULL
  
  if(Ntracks!=0) {
    print("Unpacking tracks (this may take a while...)")
    for (i in 1:Ntracks) {
      # get the detailed coordinates of the tracks
      coords <- as.numeric(unlist(trackDATA$features[[i]]$geometry)[-1])
      coords <- t(matrix(unlist(coords),nrow=2,ncol=length(coords)))
      coords <- data.frame(coords)
      names(coords) <- c("long","lat")
      
      coords$Timestamp <- strsplit(trackDATA$features[[i]]$properties$DateTimes, ",")[[1]]
      coords$Minute <- as.numeric(substr(coords$Timestamp, 15, 16))
      # coords$Timestamp <- as.POSIXct(coords$Timestamp, "%Y-%m-%dT%H:%M:%S")
      
      # downsample track records to once every 10 min
      coords$dMINUTE <- NA
      coords$dMINUTE[2:nrow(coords)] <- coords$Minute[2:nrow(coords)]-coords$Minute[2:nrow(coords)-1]
      # keep every 10th minute of the hour, remove double recordings per minute
      Include <- coords$Minute %% 10 == 0 & (coords$dMINUTE != 0 | is.na(coords$dMINUTE))
      coords <- coords[Include,]
      
      # remove unuseful columns
      coords <- dplyr::select(coords, c('long', 'lat', 'Timestamp'))
      
      TRACKS <- rbind(TRACKS,coords)
    }
    print("Successfully loaded the tracks")
    
  } else {
    print("no tracks found")
  }
}

# convert tracks to spatial object
if(!is.null(TRACKS)) {
  coordinates(TRACKS) <- ~lon+lat
  crs(TRACKS) <- crs('+proj=longlat +ellps=WGS84 +no_defs')  
  TRACKS.T <- spTransform(TRACKS, crs("+proj=utm +zone=37 +south +datum=WGS84 +units=m +no_defs"))
} else {
  print("all tracks skipped, no valid tracks loaded")
}


#------------------------------------------------------------
# indicent likelihood map
#------------------------------------------------------------

# generate the raster output
if(!is.null(DATA)) {
  print("Generating observation likelihood map")
  xrange <- bbox(DATA)[3]-bbox(DATA)[1]
  yrange <- bbox(DATA)[4]-bbox(DATA)[2]
  scaled <- sqrt(250000/(xrange*yrange))
  
  StudyArea_Raster <- raster(ncol=round(xrange*scaled), nrow=round(yrange*scaled))
  extent(StudyArea_Raster) <- extent(DATA)+max(c(xrange,yrange))
  proj4string(StudyArea_Raster) <- proj4string(DATA)
  
  binsize <- 20/scaled
  
  # generate the Raster layers
  IncDensity <- bkde2D(coordinates(DATA), 
                       bandwidth=c(binsize,binsize), 
                       gridsize=c(StudyArea_Raster@ncols,StudyArea_Raster@nrows),
                       range.x=list(extent(StudyArea_Raster)[1:2],extent(StudyArea_Raster)[3:4]))
  IncDensity$fhat[IncDensity$fhat<(max(IncDensity$fhat,na.rm=TRUE)/25)] <- NA
  IncDensity.raster <- raster(list(x=IncDensity$x1,y=IncDensity$x2,z=IncDensity$fhat))
  proj4string(IncDensity.raster) <- proj4string(StudyArea_Raster)
  
  if(!is.null(TRACKS)) {
    PatDensity <- bkde2D(coordinates(TRACKS), 
                         bandwidth=c(binsize,binsize), 
                         gridsize=c(StudyArea_Raster@ncols,StudyArea_Raster@nrows),
                         range.x=list(extent(StudyArea_Raster)[1:2],extent(StudyArea_Raster)[3:4]))
    PatDensity$fhat[PatDensity$fhat<(max(PatDensity$fhat,na.rm=TRUE)/25)] <- NA
    PatDensity.raster <- raster(list(x=PatDensity$x1,y=PatDensity$x2,z=PatDensity$fhat))
    proj4string(PatDensity.raster) <- proj4string(StudyArea_Raster)
    
    Hotspots <- IncDensity.raster/PatDensity.raster
    print("Hotspots are corrected for areas visited")
  } else {
    Hotspots <- IncDensity.raster
    print("Hotspots are not corrected for areas visited")
  }
  par(mfrow = c(1,3))
  plot(IncDensity.raster, main = "incident")
  plot(PatDensity.raster, main = 'patrolling')
  plot(Hotspots, main = 'hotspots')
  
  # create a leaflet map
  map <- leaflet() %>% 
    addTiles() %>% 
    addRasterImage(Hotspots,colors = HotspotColor) %>% 
    addCircleMarkers(data=DATA,radius = 4, weight = 2, opacity=1,color="black",fillColor="blue")
  map
}


# save status quo likelihood raster for predictive performance measurements in script C1
writeRaster(Hotspots, filename = paste0('output/status_quo_likelihood_', sc_name,'.tif'), 
            format="GTiff", overwrite = TRUE)
