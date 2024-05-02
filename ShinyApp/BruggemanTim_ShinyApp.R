### ShinyApp for single cell RNAseq data ###
############################################

######################
## Packages to load ##
######################
library(shiny)
library(Seurat)
library(bslib)
library(shinydashboard) #install.packages("shinydashboard")

options(bitmapType = "cairo") # specific for the HPC to make plots


#########################
## Loading in the data ##
#########################
#annot <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/RDSObjects/ada_and_stefNoALK_dim20_res06_ANNOT.rds")
cell <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/RDSObjects/CellsOfInterest_SLB.rds")
ALK <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/RDSObjects/stef_ALK_and_noALK_integrated_SLB.rds")


#######################################
## variables to use in the interface ##
#######################################
# features for the genes to chose from
WT_features <- rownames(cell)
ALK_features <- rownames(ALK)
allfeatures <- c(WT_features, ALK_features)


###############
## Define UI ##
###############
ui <- dashboardPage(
  dashboardHeader(title = 'ShinyApp',
                  tags$li(class="dropdown", tags$img(src="PPOL_Handjes.png"))
  ),
  dashboardSidebar(
    sidebarMenu(
      menuItem(text = "First plots",
               selectInput(
                  inputId = "plot",
                  label = "Choose a plot:",
                  choices = c("UMAP_Cluster","UMAP_GeneExpression","ViolinPlot"),
                  selected = "UMAP_Cluster"
                ),
                uiOutput("genes")
      ),
      menuItem(text = "Differential gene analyse",
               selectInput(
                 inputId = "plot2",
                 label = "choose second plot",
                 choices = c("first choice","second choice","third choice")
               )),
      menuItem(text ="Pseudotime analyse",
               textInput(inputId = "text",
                         label = "some random code")
               )
    )
  ),
  dashboardBody(
    navset_card_pill(
      tabsetPanel(
        tabPanel(title = "WT",
                 plotOutput("WT_diffplots")),
        tabPanel(title = "ALK",
                 plotOutput("ALK_diffplots"))
      )
    )
  )
  
)


#########################
## Define server logic ##
#########################
server <- function(input, output) {
  
# General output
  
  # giving gene choices only when ViolinPlot is chosen as input
  output$genes <- renderUI({
    if (input$plot == "ViolinPlot" || input$plot == "UMAP_GeneExpression") {
      selectInput(
        inputId = "gene",
        label = "Choose a gene:",
        choices = allfeatures,
      )
    }
  })
  
################################################################################
  # WT_Data output
  ################
  # Make the plot based on the selected plot type
  output$WT_diffplots <- renderPlot({
    # give the UMAP_cluster plot
    if (input$plot == "UMAP_Cluster") {
      UMAPPlot(cell, group.by = "CellType")
    }
    # give the UMAP_cluster + UMAP_GeneExpression plots
    else if (input$plot == "UMAP_GeneExpression") {
      UMAPPlot(cell, group.by = "CellType") +
        FeaturePlot(cell, features = input$gene, cols = c("lightgrey", "#FF6600","#FF0000"))
    }
    
    # Give the UMAP_cluster + ViolinPlot, without the dots
    else if (input$plot == "ViolinPlot") {
      UMAPPlot(cell, group.by = "CellType") +
        VlnPlot(cell, features = input$gene, pt.size = 0)
    }
  })
  
################################################################################
  # ALK_Data output
  #################
 output$ALK_diffplots <- renderPlot({
   # give the UMAP_cluster plot
   if (input$plot == "UMAP_Cluster") {
     UMAPPlot(ALK, group.by = "Celltype")
   }
   # give the UMAP_cluster + UMAP_GeneExpression plots
   else if (input$plot == "UMAP_GeneExpression") {
     UMAPPlot(ALK, group.by = "Celltype") +
       FeaturePlot(ALK, features = input$gene, cols = c("lightgrey", "#FF6600", "#FF0000"))
   }
   # Give the UMAP_cluster + violinplot, without the dots
   else if (input$plot == "ViolinPlot") {
     UMAPPlot(ALK, group.by = "Celltype") +
       VlnPlot(ALK, features = input$gene, pt.size = 0)
   }
 })
 
 
  
}


#########################
## Run the application ##
#########################
shinyApp(ui = ui, server = server)
