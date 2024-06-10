library(shiny)

### This gets the script of example 1 with the code ###
#######################################################
runApp("Example1.R", display.mode = "showcase")

runExample("02_text", display.mode = "showcase") # basic histogram
runExample("03_reactivity", display.mode = "showcase") # interactively changing title
runExample("04_mpg", display.mode = "showcase") # has boxplot, global variables
runExample("05_sliders", display.mode = "showcase") # slider bars
runExample("06_tabsets", display.mode = "showcase") # tabbed panels
runExample("07_widgets", display.mode = "showcase") # help text and submit button
runExample("08_html", display.mode = "showcase") # shiny app built for html
runExample("09_upload", display.mode = "showcase") # file upload wizard
runExample("10_download", display.mode = "showcase") # file download wizard
runExample("11_timer", display.mode = "showcase") # automated timer

### Build a user interface ###
##############################
library(shiny)
library(bslib)
library(bsicons)

# Define UI ----
ui <- page_sidebar(
  title = "title panel",
  sidebar = sidebar("sidebar"),
  "main contents",
  # more about cards: https://rstudio.github.io/bslib/articles/cards/?_gl=1*dzznv7*_ga*ODA3NjkyOTk2LjE3MTM2MzM1OTM.*_ga_8QJS108GF1*MTcxMzk1NjA3MS40LjEuMTcxMzk1NjE4OC4wLjAuMA..*_ga_2C0WZ1JHG0*MTcxMzk1NjA3MS40LjEuMTcxMzk1NjE4OC4wLjAuMA..
  value_box(
    title = "value box",
    value = 100,
    showcase = bsicons::icon("bar-chart"),
    theme = "teal"),
  card(
    card_header(
      class = "bg-dark",
      "card header"),
    "main contents in bordered container")
  
  
)

# Define server logic ----
server <- function(input, output) {
  
}

# Run the app ----
shinyApp(ui = ui, server = server)
