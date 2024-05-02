### ShinyApp for single cell RNAseq data ###
############################################


## Packages to load
library(shiny)
library(Seurat)
library(bslib)

options(bitmapType = "cairo") # specific for the HPC to make plots


## Loading in the data
#annot <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/RDSObjects/ada_and_stefNoALK_dim20_res06_ANNOT.rds")
cell <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/RDSObjects/CellsOfInterest_SLB.rds")
ALK <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/RDSObjects/stef_ALK_and_noALK_integrated_SLB.rds")

## variables to use in the interface


####################
## User interface ##
####################
ui <- ui <- page_navbar(
  title = "ShinyApp", # title
  bg = "#2D89C8", # navbar color
  inverse = TRUE, # make the words white on the navbar
  
  # First panel of the WT_DATA
  nav_panel(title = "WT_Data",
            sidebarLayout(
              sidebarPanel(
                selectInput(
                  inputId = "WTplot",
                  label = "Choose a plot:",
                  choices = c("UMAP_Cluster","UMAP_GeneExpression","ViolinPlot"),
                  selected = "UMAP_Cluster"
                ),
                uiOutput("WTgene"),
                # width of the sidebar
                width = 2
              ),
              position = c("left","right"),
              fluid = TRUE,
              # where the output is given
              mainPanel(
                textOutput(outputId = "WT_info"),
                headerPanel("First page content."),
                plotOutput(outputId = "WT_diffplots")
              )
            )
  ),
  ##########################################################
  
  # second panel of the ALK DATA
  nav_panel(title = "ALK_Data",
            sidebarLayout(
              sidebarPanel(
                selectInput(
                  inputId = "ALKplot",
                  label = "Choose a plot",
                  choices = c("UMAP_Cluster","UMAP_GeneExpression","ViolinPlot"),
                  selected = "UMAP_Cluster"),
                uiOutput("ALKgene"),
                width = 2
              ),
              # where the output is given
              mainPanel(
                textOutput(outputId = "ALK_info"),
                headerPanel("Second page content."),
                plotOutput(outputId = "ALK_diffplots")
              )
            )
  )
)

##################
## Server logic ##
##################
server <- function(input, output) {
  
  # giving gene choices only when ViolinPlot is chosen as input
  output$WTgene <- renderUI({
    if (input$WTplot == "ViolinPlot") {
      selectInput(
        inputId = "WTgene",
        label = "Choose gene:",
        choices = c("FOXD3","SOX10","PLP1","ERBB3","ASCL1","PHOX2B","PENK","TH","ELAVL4","ISL1","STMN2","PRPH")
      )
    }
  })

  # Make the plot based on the selected plot type
  output$WT_diffplots <- renderPlot({
    # give the UMAP_Cluster
    if (input$WTplot == "UMAP_Cluster") {
      UMAPPlot(cell, group.by = "CellType")
      
      # give the UMAP_GeneExpression
      
      
      # Give the ViolinPlot, without the dots
    } else if (input$WTplot == "ViolinPlot") {
      VlnPlot(cell, features = input$WTgene, pt.size = 0)
    }
  })

#######################################################
  
  # giving gene choices only when ViolinPlot is chosen as input
  output$ALKgene <- renderUI({
    if (input$ALKplot == "ViolinPlot") {
      selectInput(
        inputId = "ALKgene",
        label = "Choose gene:",
        choices = c("PHOX2B")
      )
    }
  })
  
  output$ALK_diffplots <- renderPlot({
    # give the UMAP_Cluster
    if (input$ALKplot == "UMAP_Cluster") {
      UMAPPlot(ALK, group.by = "Celltype")
      
      # give the UMAP_GeneExpression
      
      
      # Give the ViolinPlot, without the dots
    } else if (input$ALKplot == "ViolinPlot") {
      VlnPlot(ALK, features = input$ALKgene, pt.size = 0)
    }
  })
  
  
}


#################
## Run the app ##
#################
shinyApp(ui, server)