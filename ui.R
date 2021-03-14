library(shinyjs)#, lib.loc=tempdir())
library(leaflet)

ui <- fluidPage(
  useShinyjs(),
  sidebarLayout(
    sidebarPanel(
      
      ###############
      # LOGIN PAGE  #
      p("Please login with your Cluey credentials"),
      textOutput("loginfo"),
      tags$head(tags$style("#loginfo{color: red;}")),
      textInput('username',"Username"),
      passwordInput('password',"Password"),
      actionButton('login',"Login"),
      ###############
      
      div(
        selectInput("OntologyList", "Selected concept", c("Charcoaling", "Poaching", "Burning", "Cutting"))),
      
      tags$head(
        tags$style('#MakeMap{background-color:green; color:white; font-size:150%; margin-top:50px}')
      ),
      disabled(actionButton("MakeMap", "(Re)calculate map")),
      
      textOutput("observationwarning")
      
    ),mainPanel(
      #outputs
      h1(strong("Observation likelihood map")),
      leafletOutput("MAP", height = 1000)
    )
  )
)
