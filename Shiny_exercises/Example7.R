### This is example 7, building a movie shiny app ###
#####################################################

library(shiny)
library(ggplot2)
library(bslib)
library(DT)
library(dplyr)
library(readr)
library(tools)


# download the data
file <-"https://github.com/rstudio-education/shiny-course/raw/main/movies.RData"
destfile <- "./Data/movies.RData"
download.file(file,destfile)

# load the data
load("./Data/movies.RData")

# Make variables to use in interface
n_total <- nrow(movies)
all_studios <- sort(unique(movies$studio))
min_date <- min(movies$thtr_rel_date)
max_date <- max(movies$thtr_rel_date)
movies <- movies %>% # ratio of critics and audience scores
  mutate(score_ratio = audience_score / critics_score)

# User interface
ui <- page_sidebar(
  title = "Movie browser",
  sidebar = sidebar(
    # To select the y-axis input
    selectInput(
      inputId = "y_axis",
      label = "Y-axis:",
      choices = c("imdb_rating", "imdb_num_votes", "critics_score", "audience_score","runtime"),
      selected = "audience_score"
    ),
    # To select x-axis input
    selectInput(
      inputId = "x_axis",
      label = "X-axis:",
      choices = c("imdb_rating", "imdb_num_votes", "critics_score", "audience_score","runtime"),
      selected = "critics_score"
    ),
    # select the color
    selectInput(
      inputId = "z",
      label = "Color by:",
      choices = c("title_type", "genre", "mpaa_rating", "critics_rating", "audience_rating"),
      selected = "mpaa_rating"
    ),
    # select brightness
    sliderInput(
      inputId = "alpha",
      label = "Alpha:",
      min = 0,
      max = 1,
      value = 0.5
    ),
    # select to show data table
    checkboxInput(
      inputId = "show_data",
      label = "Show the data table",
      value = TRUE
    ),
    # type text with the total amount of rows
    HTML(paste("Enter a value between 1 and", n_total)),
    # ask for numeric input how much rows you want in the table
    numericInput(
      inputId = "n",
      label = "Sample size:",
      value = 5,
      min = 1, max = n_total,
      step = 1
    ),
    # select studio's
    selectInput(
      inputId = "studio",
      label = "Select studio:",
      choices = all_studios,
      selected = "20th Century Fox",
      multiple = TRUE
    ),
    HTML(paste0("Movies released since the following date will be plotted. 
                 Pick a date between ", min_date, " and ", max_date, ".")),
    dateRangeInput(
      inputId = "date",
      label = "Choose data range",
      start = "2013-01-01",
      end = "2014-01-01"
    ),
    checkboxGroupInput(
      inputId = "type",
      label = "Select movie type(s):",
      choices = levels(movies$title_type),
      selected = levels(movies$title_type)
    ),
    
    #Download button
    radioButtons(
      inputId = "filetype",
      label = "select filetype:",
      choices = c("csv","tsv"),
      selected = "csv"
    ),
    checkboxGroupInput(inputId = "selected_var",
                       label = "Select variables:",
                       choices = names(movies),
                       selected = c("title")
    ),
    
    # Select which types of movies to plot
    selectInput(
      inputId = "selected_type",
      label = "Select movie type2.0:",
      choices = levels(movies$title_type),
      selected = "Feature Film"
    ),
    
    # text for plottitle
    textInput(
      inputId = "plot_title",
      label = "Plot title",
      placeholder = "Enter text for plot title"
    ),
    
    # to change the title
    actionButton(
      inputId = "update_plot_title",
      label = "Update plot title"
    ),
    
    # for table output
    numericInput(
      inputId = "n_rows",
      label = "How many rows do you want to see?",
      value = 10
    ),
    # for showing message datatable
    actionButton(
      inputId = "button",
      label = "Show"
    )
  ),
  # output: Show scatterplot
  card(#plotOutput(outputId = "scatterplot"),
      plotOutput(outputId = "scatterplot"), #hover = "plot_hover"), # could also use brush = 
      # plotOutput(outputId = "densityplot", height = 200),
      # dataTableOutput(outputId = "datatable"),
      #tableOutput(outputId = "summarytable"),
      #textOutput(outputId = "correlation")
      #dataTableOutput(outputId = "moviestable")
      #textOutput(outputId = "avg_x"),
      #textOutput(outputId = "avg_y"),
      #htmlOutput(outputId = "avgs"),# shorter version of the two before
      #verbatimTextOutput(outputId = "lmoutput")
      #HTML("select filetype and variables, then hit 'Download data'."),
      #downloadButton("download_data","Download data")
      #uiOutput(outputId = "n"),# print number of obs plotted
      #dataTableOutput(outputId = "moviestable2"),
      tableOutput(outputId = "datatable2")
      )
)

# Server logic
server <- function(input, output, session) {
  # scatterplot output
  output$scatterplot <- renderPlot({
    ggplot(data=movies, aes_string(x = input$x_axis, 
                                   y = input$y_axis,
                                   color = input$z)) +
      geom_point(alpha = input$alpha, size = input$size)
  })
  
  # densityplot output
  output$densityplot <- renderPlot({
    ggplot(data = movies, aes_string(x=input$x_axis)) +
      geom_density()
  })
  
  # table output
  output$datatable <- DT::renderDataTable({
    # hold back output from being calculated if input is missing
    req(input$n)
    # filter on date
    movies_from_selected_date <- movies %>%
      filter(thtr_rel_date >= as.POSIXct(input$date))
    # filter on studios
    movies_from_selected_studios <- movies %>%
      filter(studio %in% input$studio) %>%
      select(title:studio)
    # amount of movies from numeric input
    movies_sample <- movies %>%
      sample_n(input$n) %>%
      select(title:studio)
    # output for the table
    datatable(data = movies_sample,
                  options = list(pageLength = 10),
                  rownames = FALSE)
  })
  
  # make summary table
  output$summarytable <- renderTable({
      movies %>%
        filter(title_type %in% input$type) %>%
        group_by(mpaa_rating) %>%
        summarise(mean_score_ratio = mean(score_ratio), SD = sd(score_ratio), n = n())
    },
    # make the table a bit more beautiful
    striped = TRUE, # altering color between rows
    spacing = "l", # larger row heights
    align = "lccr", # align left, center, center, right
    digits = 4, # number of decimal places to display
    width = "90%", # width of table output
    caption = "Score ratio (audience / critics' scores) summary statistics by MPAA rating."
  )
  
  # create text output stating the correlation between two plots
  output$correlation <- renderText({
    r <- round(cor(movies[,input$x_axis], movies[,input$y_axis], use = "pairwise"), 3)
    paste0(
      "Correlation = ", r,
      ". Note: If the relationship between the two variables is not linear, the correlation coefficient will not be meaningful."
    )
  })
  
  # Print data table
  output$moviestable <- renderDataTable({
    nearPoints(movies, input$plot_hover) %>% # brushedPoints if plot_brush
      select(title, audience_score, critics_score)
  })
  
  #output avg_x
  output$avg_x <- renderText({
    avg_x <- movies %>% pull(input$x_axis) %>% mean() %>% round(2)
    paste("Average", input$x_axis, "=", avg_x)
  })
  
  #output avg_y
  output$avg_y <- renderText({
    avg_y <- movies %>% pull(input$y_axis) %>% mean() %>% round(2)
    paste("Average", input$y_axis, "=", avg_y)
  })
  
  # shorter averages display version
  output$avgs <- renderUI({
    avg_x <- movies %>% pull(input$x_axis) %>% mean() %>% round(2)
    avg_y <- movies %>% pull(input$y_axis) %>% mean() %>% round(2)
    str_x <- paste("Average", input$x_axis, "=", avg_x)
    str_y <- paste("Average", input$y_Axis, "=", avg_y)
    HTML(paste(str_x, str_y, sep = '<br/>'))
  })
  
  # regression output
  output$lmoutput <- renderPrint({
    x <- movies %>% pull(input$x_axis)
    y <- movies %>% pull(input$y_axis)
    print(summary(lm(y ~ x, data=movies)), digits=3, signif.stars = FALSE)
  })
  
  # download file
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("movies.", input$filetype)
    },
    content = function(file) {
      if(input$filetype == "csv"){
        write_csv(movies %>% select(input$selected_var), file)
      }
      if(input$filetype == "tsv"){
        write_tsv(movies %>% select(input$selected_var), file)
      }
    }
  )
  
  # Create a subset of data filtering for chosen title types
  movies_subset <- reactive({
    req(input$selected_type)
    filter(movies, title_type %in% input$selected_type)
  })
  
  # make evenreactive title, change this with previous isolate({pretty_plot_title})
  new_plot_title <- eventReactive(
    eventExpr = input$update_plot_title,
    valueExpr = {
      toTitleCase(input$plot_title)
    }
  )
  # Create scatterplot
  output$scatterplot <- renderPlot({
    ggplot(data = movies_subset(),aes_string(x = input$x_axis, y = input$y_axis)) +
      geom_point() +
      labs(title = isolate({pretty_plot_title()})) # isolate function is to only change title if something from the variables is changed by the user
  })
  
  # server - Print number of movies plotted
  output$n <- renderUI({
    HTML(paste0(
      "The plot displays the relationship between the <br>
              audience and critics' scores of <br>",
      nrow(movies_subset()),
      " <b>", input$selected_type, "</b> movies."
    ))
  })
  
  # Create reactive data frame
  movies_selected <- reactive({
    movies %>% select(input$selected_var)
  })
  
  # Create data table
  output$moviestable2 <- DT::renderDataTable({
    req(input$selected_var)
    datatable(
      data = movies_selected(),
      options = list(pageLength = 10),
      rownames = FALSE
    )
  })
  
  # Convert plot_title toTitleCase
  pretty_plot_title <- reactive({
    toTitleCase(input$plot_title)
  })
  # Create descriptive text
  output$description <- renderText({
    paste0("The plot above titled '", pretty_plot_title(), "' visualizes the relationship between ", input$x, " and ", input$y, ", conditional on ", input$z, ".")
  })
  
  # for datatable 2
  # Pring a message to the console every time button is pressed;
  observeEvent(input$button, {
    cat("Showing", input$n_rows, "rows\n")
  })
  
  # Take a reactive dependency on input$button,
  # but not on any of the stuff inside the function
  df <- eventReactive(input$button, {
    head(movies, input$n_rows)
  })
  output$datatable2 <- renderTable({
    df()
  })
}

# Run the app
shinyApp(ui, server)