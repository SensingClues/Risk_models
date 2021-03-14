library(shiny)
library(shinyjs)#, lib.loc=tempdir())
# focus sparql endpoint
library(httr)
library(jsonlite)
library(SPARQL)
# GIS libraries
library(raster)
library(sp)
library(rgdal)
library(leaflet)
# tidyverse
library(tidyverse)
library(tidymodels)


server <- function(input, output) {
  
  Env <- reactiveValues(
    ObsColor = c("#040404B3","#080918B3","#0E0D24B3","#150F2EB3","#1D1135B3","#24123CB3","#2C1242B3","#341348B3","#3C134EB3","#451353B3","#4D1259B3",
                 "#56125DB3","#5F1162B3","#681066B3","#701069B3","#79106DB3","#82106FB3","#8A1172B3","#931373B3","#9B1674B3","#A31A75B3","#AB1E75B3",
                 "#B32375B3","#BA2973B3","#C12F71B3","#C8356FB3","#CF3B6BB3","#D64267B3","#DC4962B3","#E2505BB3","#E85752B3","#ED5F48B3","#F2673AB3",
                 "#F37133B3","#F47B2CB3","#F58426B3","#F58E23B3","#F69622B3","#F79F25B3","#F7A82CB3","#F7B134B3","#F8B93EB3","#F8C149B3","#F8CA54B3",
                 "#F9D25FB3","#F9DB6BB3","#FAE377B3","#FBEC84B3","#FDF490B3","#FFFE9EB3"),
    
    maxPatFreq = 1.12e-08,
    loggedin=FALSE
  )
  
  observeEvent(input$login,{
    # login
    url.login <- "https://focus.sensingclues.org/api/auth/login"
    # put your cluey/focus credentials here
    username = input$username
    password = input$password
    json_body <- jsonlite::toJSON(list(username = username,password=password), auto_unbox = TRUE)
    # we set up an authenticated session
    rl <- POST(url.login, body = json_body, encode = "raw",content_type_json())
    
    # for debugging
    print("Login attempt")
    
    if(status_code(rl)==200) {
      output$loginfo <- renderText(paste("Logged in as",content(rl)$username))
      hide("username")
      hide("password")
      hide("login")
      enable("MakeMap")
      Env$loggedin = TRUE
      
      # for debugging
      print("Successfully logged in")
    } else {
      output$loginfo <- renderText(content(rl)$message)
      
      # for debugging
      print("Login failed")
    }
    
  })
  
  observeEvent(input$MakeMap,{    
    shinyjs::disable("MakeMap")
  })
  
  re <- eventReactive(input$MakeMap,{
    
    ####################
    # ADD FILTERS HERE #
    ####################
    # location
    if (!is.null(input$MAP_bounds)) {
      bounds <- input$MAP_bounds
    } else {
      bounds <- list(north=-3,east=41,south=-6,west=38)
    }
    Ontology <- c("https://sensingclues.poolparty.biz/SCCSSOntology/97",
                  "https://sensingclues.poolparty.biz/SCCSSOntology/93",
                  "https://sensingclues.poolparty.biz/SCCSSOntology/98",
                  "https://sensingclues.poolparty.biz/SCCSSOntology/99")[c("Charcoaling", "Poaching", "Burning", "Cutting") %in% input$OntologyList]
    print("Filters defined")
    #print(paste("Concept =",ontology.tbl$Concept))
    print("Spatial bounds:")
    print(bounds)
    print(paste("Date range: from",Sys.Date() - 31,"to",Sys.Date()))
    
    ###############
    # END FILTERS #
    ###############
    
    # select the proper source URL
    URL <- "https://focus.sensingclues.org/"
    url.incidents <- paste0(URL,"api/map/all/default/0/features")
    url.tracks <- paste0(URL,"api/map/all/track/0/features")
    
    #########
    # SETUP #
    #########
    
    # set coordinate reference system
    crs_kenya <- crs("+proj=utm +zone=37 +south +datum=WGS84 +units=m +no_defs")
    
    # retrieve the area of interest (AOI) polygon
    aoi <- readOGR(dsn = "data", layer = "study_area") %>% 
      spTransform(crs_kenya)
    
    # create the desired output raster from the aoi Polygon object
    out <- raster(aoi, resolution = 1400) # 1.4 x 1.4 km spatial resolution (from Ea)
    
    # retrieve the calculated environmental feature values per cell in the AOI
    features <- read.csv('data/location_features_AOI.csv')
    
    # function to create pseudo-absence locations
    get_pseudo_absence_locations <- function(area_sp, area_grid, presence_data_sp, nabsent = NULL, 
                                             seed = 555, patrol_sp = NULL, min_revisit = 10){
      
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
      
      if (is.null(patrol_sp)){
        # randomly sample pseudo-absence cell numbers from the non-incident cells
        cell_opt <- setdiff(cell_aoi, cell_pr)
        
      } else {
        
        if (crs(area_grid)@projargs != crs(patrol_sp)@projargs){
          patrol_sp <- spTransform(patrol_sp, crs(area_grid))
        }
        
        cell_pat <- as.data.frame(cellFromXY(area_grid, patrol_sp)) # points outside of area_grid get NA value
        names(cell_pat) <- 'cell_id'
        cell_pat <- cell_pat %>% 
          group_by(cell_id) %>% 
          summarise(count = n()) %>% 
          filter(count > min_revisit) %>% 
          filter(cell_id %in% cell_aoi) # remove cell_ids outside of AOI
        cell_pat <- cell_pat$cell_id
        
        cell_opt <- setdiff(cell_pat, cell_pr)
      }
      
      set.seed(seed)
      if (is.null(nabsent)){ # get equal amount of absent as present
        cell_ab <- sample(cell_opt, length(cell_pr))
      } else {
        cell_ab <- sample(cell_opt, nabsent)
      }
      
      
      return(list(cell_ab, cell_pr, cell_aoi))
    }
    
    
    ######################################
    # get the list with observation data #
    ######################################
    
    q <- paste0('
    {"filters":
    {"geoQuery":
    {"operator":"intersects",
    "mapBounds":{"south":',bounds$south,',"west":',bounds$west,',"north":',bounds$north,',"east":',bounds$east,'},
    "drawings":[]},
    "entities":["Observation"],
    "concepts":["',Ontology,'"],
    "dataSources":["',Env$projects$code[Env$projects$name==input$project],'"],
    "dateTimeRange":{"to":"',Sys.Date(),'T22:00:00.000Z","from":"',Sys.Date()-31,'T23:00:00.000Z"}},
    "options":{"start":0,"pageLength":20},"start":1,"pageLength":500}
    ')
    obsDATA.L <- POST(url.incidents, body=q, encode="raw", content_type_json())
    obsDATA <- content(obsDATA.L)
    #print(obsDATA)
    Env$obsDATA <- obsDATA
    Nobs <- length(obsDATA$features)
    
    DATA <- NULL
    if(Nobs!=0) {
      print(paste0('There have been ', Nobs, ' observations collected in this session'))
      print("Getting coordinates of observations")
      for(i in 1:Nobs) {
        coords <- unlist(obsDATA$features[[i]]$geometry$coordinates)
        # Only read JSON files
        DATA <- rbind(DATA,coords)
      }
      DATA <- as.data.frame(DATA, row.names = NULL)
      names(DATA) <- c("lon","lat")
      row.names(DATA) <- NULL
      coordinates(DATA) <- ~lon+lat
      proj4string(DATA) <- "+proj=longlat +ellps=WGS84 +no_defs"
    } else {
      print("No observations found")
    }
    
    if (Nobs < 25){
      output$observationwarning <- renderText('Warning: The likelihood map is based on few observations.')
    }
    
    
    ###################################
    # get the list with tracking data #
    ###################################
    
    queryNEW <- paste0('
    {"filters":
    {"geoQuery":
    {"operator":"intersects",
    "mapBounds":{"south":',bounds$south,',"west":',bounds$west,',"north":',bounds$north,',"east":',bounds$east,'},
    "drawings":[]},
    "dateTimeRange":{"to":"',Sys.Date(),'T24:00:00.000Z","from":"',Sys.Date()-31,'T00:00:00.000Z"}}}
    ')
    trackDATA.L <- POST(url.tracks, body=queryNEW, encode="raw", content_type_json())
    trackDATA <- content(trackDATA.L)
    Ntracks <- length(trackDATA$features)
    print(paste0('Number of tracks:',Ntracks))
    
    TRACKS <- NULL
    if(Ntracks!=0) {
      print("Unpacking tracks (this may take a while...)")
      for (i in 1:Ntracks) {
        # get the detailed coordinates of the tracks
        coords <- as.numeric(unlist(trackDATA$features[[i]]$geometry)[-1])
        coords <- t(matrix(unlist(coords),nrow=2,ncol=length(coords)))
        coords <- data.frame(coords)
        names(coords) <- c("lon","lat")
        
        #DateTime <- strsplit(trackDATA$features[[i]]$properties$DateTimes, ",")[[1]]
        #Date <- as.Date(substr(DateTime, 1,10))
        #Time <- substr(DateTime, 12,19)
        #DateTime <- paste(Date,Time)
        
        # add timestamp data
        # first try this one:
        timestamps <- strsplit(trackDATA$features[[i]]$properties$DateTimes, ",")[[1]]
        
        # if that doesn't work:
        if(is.null(timestamps)) {
          coords$time <- NA
          coords$timeslot <- NA
        } else {
          coords$time <- timestamps
          # remove duplicates
          coords <- coords[!duplicated(coords$time),]
          # calculate the time difference between subsequent fixes
          coords$time <- `substr<-`(coords$time,11,11," ")
          coords$time <- as.POSIXlt(coords$time)
          coords$dt <- c(as.numeric(coords$time[-nrow(coords)] - coords$time[-1]),0)
          # order the data so that it increases in time
          coords <- coords[order(coords$time),]
          # devide dataset into 5-minute time slots and save only 1 5-minute slot
          coords$timeslot <- cumsum(coords$dt)%/%(60*5)   
          coords <- coords[!duplicated(coords$timeslot),-4]
        }
        
        TRACKS <- rbind(TRACKS,coords)
        
      }
      print("Successfully loaded the tracks")
      # make the data spatial
      if(!is.null(TRACKS)) {
        coordinates(TRACKS) <- ~lon+lat
        proj4string(TRACKS) <- "+proj=longlat +ellps=WGS84 +no_defs"
        # TRACKS.T <- spTransform(TRACKS, crs("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
      } else {
        print("All tracks skipped, no valid tracks loaded")
      }
    } else {
      print("No tracks found")
    }
    
    #############################################################
    # retrieve features, train the model and predict likelihood #
    #############################################################
    
    if (is.null(input$MAP_bounds)&!is.null(DATA)) {
      bounds <- list(north=bbox(DATA)[4],east=bbox(DATA)[3],south=bbox(DATA)[2],west=bbox(DATA)[1])
    }
    
    # m makes the map in a leaflet format
    m <- leaflet() %>% 
      addTiles()
    
    # generate the raster output
    if(!is.null(DATA)) {
      print("Generating observation likelihood map")
      
      # assign each location in AOI to absent/present/other for the training set
      tr_cell_list <- get_pseudo_absence_locations(aoi, out, DATA, patrol_sp = TRACKS, 
                                                   min_revisit = 3)
      tr_cell_ab <- tr_cell_list[[1]]
      tr_cell_pr <- tr_cell_list[[2]]
      tr_cell_aoi <- tr_cell_list[[3]]
      
      # add location labels to output raster
      out$train <- NA
      out$train[tr_cell_aoi] <- rep(2, length(tr_cell_aoi))
      out$train[tr_cell_pr] <- rep(1, length(tr_cell_pr))
      out$train[tr_cell_ab] <- rep(0, length(tr_cell_ab))
      # m <- addRasterImage(m,ratify(out$train))
      
      # create features dataframe, rows representing the grid cells in the AOI
      df <- data.frame(cell_id = tr_cell_aoi, train = out$train[tr_cell_aoi]) %>% 
        mutate(train_label = ifelse(train == 0, 'pseudo-absent', 
                                    ifelse(train == 1, 'incident', 'other')),)
      df <- merge(df, features, by = 'cell_id')
      print(head(df, 3))
      
      # retrieve the training data
      train <- df %>% 
        filter(train != 2) %>% 
        mutate(train = factor(train, levels = c(0,1), labels = c("Unharmed", "Incident")))
      print(table(train$train))
      
      # preprocess the training data (remove NAs etc.)
      recipe <-
        train %>%
        recipe(train ~ .) %>%
        step_rm(train_label) %>% 
        step_naomit(all_predictors()) %>% 
        step_corr(all_numeric(), -all_outcomes(), threshold = 0.9) %>%
        step_zv(all_predictors()) %>%
        step_nzv(all_predictors()) %>%
        prep()
      baked_train <- bake(recipe, train)
      
      # train a Random Forest model
      mod_final <- rand_forest(mode = "classification", 
                               mtry=5,
                               trees=100,
                               min_n=5) %>%
        set_engine("randomForest") %>%
        fit(train ~ . -cell_id, data = baked_train)
      print(mod_final)
      
      
      # preprocess all input locations
      baked_all <- bake(recipe, df %>% dplyr::select(-train))
      
      # predict incident likelihood at each location
      preds_aoi <- predict(mod_final, new_data = baked_all, type = "prob") # get probabilities
      preds_aoi <- cbind(preds_aoi, cell_id =baked_all$cell_id)
      print(head(preds_aoi))
      out$likelihood <- NA
      out$likelihood[preds_aoi$cell_id] <- preds_aoi$.pred_Incident
      
      # incident classification map
      preds_aoi <- mutate(preds_aoi, class = ifelse(.pred_Incident < 0.15, 1, 
                                                    ifelse(.pred_Incident < 0.5, 2, 
                                                           ifelse(.pred_Incident < 0.85, 3, 4))))
      out$prediction <- NA
      out$prediction[preds_aoi$cell_id] <- preds_aoi$class
      
      # color palettes
      pal_lik <- colorNumeric(palette = Env$ObsColor, domain = values(out$likelihood),
        na.color = 'transparent')
      pal_lik_rev <- colorNumeric(palette = Env$ObsColor, domain = values(out$likelihood),
        na.color = 'transparent', reverse = TRUE)
      
      # create map
      m <- m %>% 
        addRasterImage(out$likelihood, group = 'likelihood', colors = pal_lik, opacity = 0.75) %>% 
        addLegend(pal = pal_lik_rev, values = values(out$likelihood),
                  labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)),
                  title = paste0(input$OntologyList, " </br> probability"), group = 'likelihood') %>% 
        addRasterImage(out$prediction, group = 'prediction', opacity = 0.75, 
                       colors = c("#008B45", "#43CD80", "#FF6A6A", "#FF0000")) %>%  # 'springgreen4', 'seagreen3', 'indianred1', 'red'
        addLegend(colors = c("#008B45", "#43CD80", "#FF6A6A", "#FF0000"),
                  labels = c('very unlikely', 'unlikely', 'likely', 'very likely'),
                  title = paste0(input$OntologyList, " </br> incident"), group = 'prediction') %>%
        addLayersControl(
          overlayGroups = c("likelihood", "prediction"),
          options = layersControlOptions(collapsed = FALSE)
        ) %>% 
        hideGroup('prediction') %>% 
        addScaleBar(
          position = c("bottomleft"),
          options = scaleBarOptions()
        )
      # m <- addCircleMarkers(m, data=DATA, radius = 4, weight = 2, opacity=1, color="black", fillColor="blue")
    }
    
    # fix the boundaries that are displayed on the interactive map
    latRng <- range(bounds$north, bounds$south)
    lngRng <- range(bounds$east, bounds$west)
    m[[1]]$limits$lat <- latRng
    m[[1]]$limits$lng <- lngRng
    
    shinyjs::enable("MakeMap")
    
    return(m)
    
  })
  
  output$MAP <- renderLeaflet({
    re()
  })
  
}