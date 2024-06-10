### This is examle 8, about customizing UI ###
##############################################

library(shiny)
library(bslib)

load("./Data/movies.RData")

# Define UI with tags
ui <- page_fluid( theme = bs_theme(preset = "vapor"),
  tags$h1("First level heading"),
  tags$h2("Second level heading"),
  tags$h3("Third level heading"),
  tags$b("This is bold text"),
  tags$br(), # line break
  tags$a("This app is built with Shiny.", href = "https://shiny.posit.co/"),
  tags$p( # paragraph
    "Lorem ipsum",
    tags$em("dolor"), # schuin
    "sit amet,",
    tags$b("consectetur"), # bold
    "adipiscing elit."
  ),
  tags$code("Monospace text"),
  tags$p(
    "These data were obtained from",
    tags$a("IMBD", href = "http://www.imbd.com/"), "and",
    tags$a("Rotten Tomatoes", href = "https://www.rottentomatoes.com/"), "."
  ),
  tags$p("The data represent", nrow(movies), "randomly sampled movies released between 1972 to 2014 in the United States."),
  titlePanel("An image"),
  tags$img(height = 100, width = 300, src = "./Data/shinyimage.jpg"),
  layout_columns(
    card("Width 3 + height"), # Column 1
    card("Width 4 + heigt"), # Column 2
    card("Width 5 + height"), # Column 3
    col_widths = c(3, 4, -10, 5, -2),
    row_heights = c(1,2)
  ),
  navset_card_underline(
    
    nav_panel("Plot", plotOutput("plot")),
    
    nav_panel("Summary", tableOutput("summary")),
    
    nav_panel("Data", DT::dataTableOutput("data")),
    
    nav_panel(
      "Reference",
      markdown(
        glue::glue(
          "These data were obtained from [IMDB](http://www.imdb.com/) and [Rotten Tomatoes](https://www.rottentomatoes.com/).
  
        The data represent {nrow(movies)} randomly sampled movies released between 1972 to 2014 in the United States.
        "
        )
      )
    )
  )
)

# Define server function that does nothing :)
server <- function(input, output, session) {
  bs_themer()
}

# Create the app object
shinyApp(ui = ui, server = server)