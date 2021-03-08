# Perfoming necessary preparations for the project:
#
# Create a spatial object representing the area of interest, using a self-defined
# range of coordinates (format: longitude, latitude).
#
# Retrieve locally stored tracks (2010-2020), downsample them to have a 10 min interval
# and store them for later use.


# general libraries
library(tidyverse)
# GIS libraries
library(sp)
library(rgdal)


# the tracks object is stored in a data/processed folder
if (! dir.exists('data/processed')){
  dir.create('data/processed')
}


#------------------------------------------------------------
# Create area or interest (AOI) polygon
#------------------------------------------------------------

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


#------------------------------------------------------------
# tracks data collection: local track folders (2010-2020)
#------------------------------------------------------------

TRACKS <- NULL

# for each file in each folder:
FOLDERS <- list.dirs()
FOLDERS <- FOLDERS[grepl('tracks/ST', FOLDERS, fixed = TRUE)]
for (j in 1:length(FOLDERS)) {
  print(paste0('folder ', j, ' out of ', length(FOLDERS)))
  files <- list.files(FOLDERS[j])[-1]
  files <- substr(files,1,nchar(files)-4)
  for (i in 1:length(files)) {
    FILENAME <- paste(FOLDERS[j],"/",files[i],".LOG",sep="")
    track <- read.table(FILENAME,skip=1)
    track <- as.character(track[substring(track$V1,1,6)=="$GPGGA",])
    TRACK <- NULL
    for(a in 1:length(track)) {
      TRACK <- rbind(TRACK,strsplit(track[a],",")[[1]])
    }
    TRACK <- as.data.frame(TRACK)
    # only use tracks with more than 5 recordings
    if (dim(TRACK)[1] >= 5){
      names(TRACK) <- c("SentenceID","Time","Latitude","N_or_S","Longitude","E_or_W","FixQuality",
                        "NofSatellites","HDOP","Altitude","Altitude_unit","HeightAbvWGS84","unit","DGPSref","Checksum")
      TRACK$Time <- as.numeric(as.character(TRACK$Time))
      TRACK$Latitude <- as.numeric(as.character(TRACK$Latitude))
      TRACK$Longitude <- as.numeric(as.character(TRACK$Longitude))
      TRACK$Altitude <- as.numeric(as.character(TRACK$Altitude))
      
      TRACK$ID <- FILENAME
      TRACK$Year <- as.numeric(paste("20",substr(files[i],2,3),sep=""))
      TRACK$Month <- as.numeric(substr(files[i],4,5))
      TRACK$Day <- as.numeric(substr(files[i],6,7))
      TRACK$Hour <- TRACK$Time %/% 10000
      TRACK$Minute <- TRACK$Time %% 10000
      TRACK$Second <- TRACK$Minute %% 100
      TRACK$Minute <- TRACK$Minute %/% 100
      TRACK$Timestamp <- as.POSIXct(paste0(TRACK$Year,"-",TRACK$Month,"-",TRACK$Day," ",
                                           TRACK$Hour,":",TRACK$Minute,":",TRACK$Second),tz="UTC")
      # ----------------
      # # Fix date when time passes midnight
      # TRACK$dTIME <- NA
      # TRACK$dTIME[2:nrow(TRACK)] <- TRACK$Time[2:nrow(TRACK)]-TRACK$Time[2:nrow(TRACK)-1]
      # newID <- c(TRUE,!TRACK$ID[2:nrow(TRACK)] == TRACK$ID[2:nrow(TRACK)-1])
      # DateChange <- TRACK$dTIME < 0 & !is.na(TRACK$dTIME)&!newID
      # i <- 1
      # while(!i > nrow(TRACK)) {
      #   if (newID[i]) {
      #     check <- FALSE
      #   }
      #   if (DateChange[i]) {
      #     check <- TRUE
      #     print(TRACK$Time[i])
      #   }
      #   if (check) {DateChange[i] <- TRUE}
      #   
      #   i <- i + 1
      #   #print(i)
      # }
      # TRACK$Time[DateChange] <- TRACK$Time[DateChange] + 60*60*24
      # TRACK$dTIME[2:nrow(TRACK)] <- TRACK$Time[2:nrow(TRACK)]-TRACK$Time[2:nrow(TRACK)-1]
      # TRACK$dTIME[newID] <- NA
      # ---------------
      
      # get longitude and latitude in WGS84 format
      testLat <- TRACK$Latitude%/%100+TRACK$Latitude%%100/60
      testLong <- TRACK$Longitude%/%100+TRACK$Longitude%%100/60
      testLat[TRACK$N_or_S=="S"] <- -testLat[TRACK$N_or_S=="S"]
      
      TRACK$long <- testLong
      TRACK$lat <- testLat
      testLong[TRACK$E_or_W=="W"] <- -testLong[TRACK$E_or_W=="W"]
      
      # downsample track records to once every 10 min
      TRACK$dMINUTE <- NA
      TRACK$dMINUTE[2:nrow(TRACK)] <- TRACK$Minute[2:nrow(TRACK)]-TRACK$Minute[2:nrow(TRACK)-1]
      # keep every 10th minute of the hour, remove double recordings per minute
      Include <- TRACK$Minute %% 10 == 0 & (TRACK$dMINUTE != 0 | is.na(TRACK$dMINUTE))
      TRACK <- TRACK[Include,]
      
      # remove unuseful columns
      TRACK <- dplyr::select(TRACK, c('ID', 'long', 'lat', 'Timestamp', 'Year'))
      
      TRACKS <- rbind(TRACKS,TRACK)
    }
    
  }
  
}

write.csv(x = TRACKS, file = 'data/processed/tracks_2015_to_2019.csv', row.names = FALSE)

