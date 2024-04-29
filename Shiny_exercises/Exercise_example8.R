library(shiny)
library(bslib)

# Define UI with conditionalPanel
ui <- page_sidebar(
  title ="Random number generator",
  sidebar = sidebar(
    width = 325,
    # checkbox to get to the sliderinput
    checkboxInput(
      "select_digits", "Click to select number of digits",
      value = FALSE
    ),
    # slider will show if checkbox has value TRUE
    conditionalPanel(
      condition = "input.select_digits == true",
      sliderInput("digit_count", "How many digits?",
                  min = 1, max = 10, value = 4
      )
    ),
    # actionbutton to generate new number
    actionButton("go", "Generate new random number")
  ),
  
  # this is the output
  "Your random number is",
  h4(textOutput("random_number"))
)

# Define server logic that generates a random number based on user input
server <- function(input, output, session) {
  output$random_number <- renderText({
    input$go
    raw <- runif(1)
    digits <- if (input$select_digits == FALSE) {
      sample(1:10, size = 1)
    } else {
      input$digit_count
    }
    round(raw * 10^digits)
  })
}

shinyApp(ui = ui, server = server)