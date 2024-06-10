### This is an exercise from example 3, add control widgets ###
###############################################################

library(shiny)

# Define UI ----
ui <- page_sidebar(
  title = "CensusVis",
  sidebar = sidebar(
    "Create demographic maps with information from the 2010 US Census.",
    selectInput(
      "variable",
      label =  "Choose a variable to display",
      choices = list("Percent White" = 1, "Percent Black" = 2, "Percent Hispanic" = 3, "Percent Asian" = 4),
      selected = 1),
    sliderInput(
      "range",
      "Range of interest:",
      min = 0, max = 100, value = c(0, 100))
  )
)
  

# Define server logic ----
server <- function(input, output) {
  
}

# Run the app ----
shinyApp(ui = ui, server = server)


###This is more the solution ###
################################
#ui <- page_sidebar(
#  title = "censusVis",
#  sidebar = sidebar(
#    helpText(
#      "Create demographic maps with information from the 2010 US Census."
#    ),
#    selectInput(
#      "var",
#      label = "Choose a variable to display",
#      choices = 
#        list(
#          "Percent White", 
#          "Percent Black", 
#          "Percent Hispanic", 
#          "Percent Asian"
#        ),
#      selected = "Percent White"
#    ),
#    sliderInput("range",
#                label = "Range of interest:",
#                min = 0, max = 100, value = c(0, 100))
#  )
#)
