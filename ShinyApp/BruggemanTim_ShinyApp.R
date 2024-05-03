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
  # headerbar
  dashboardHeader(title = 'ShinyApp',
                  # making the logo
                  tags$li(class="dropdown", imageOutput("logo", height = 50, width = 50))
  ),
  # sidebarmenu
  dashboardSidebar(
    # sidebarmenu with identification of the different menu's via menuItem
    sidebarMenu(id = "sidebarid",
                menuItem(text = "first plots", tabName = "first"),
                # Choose between the three plots to visualize
                conditionalPanel("input.sidebarid == 'first'",
                                 selectInput(
                                   inputId = "plot",
                                   label = "choose a plot:",
                                   choices = c("UMAP_Cluster","UMAP_GeneExpression","ViolinPlot"),
                                   selected = "UMAP_Cluster"
                                 ),
                                 uiOutput("genes")
                ),
                # Second menuItem
                menuItem("Differential gene analyse", tabName = "diff_gene"),
                # .... to differentiate gene analys, now as example filled in with selectinput
                conditionalPanel("input.sidebarid == 'diff_gene'",
                                 selectInput(
                                   inputId = "plot2",
                                   label = "choose second plots",
                                   choices = c("first plot", "second plot"," third plot")
                                 )
                ),
                # Third menuItem
                menuItem("Pseudotime analyse", tabName = "pseudo")
    )
  ),
  # Output in the body
  dashboardBody(
    tabItems(
      # output for first menu: the plots visualization
      tabItem(
        tabName = "first",
        navset_pill(
          # Output tab for WT_Data
          tabPanel(title = "WT", 
                   plotOutput("WT_diffplots", height = 500, width = 1250)),
          # Output tab for ALK_Data
          tabPanel(title = "ALK", 
                   plotOutput("ALK_diffplots"))
        )
      ),
      # output for second menu: differential gene analyse
      tabItem(tabName = "diff_gene"),
      # output for third menu: pseudotime analyse
      tabItem(tabName = "pseudo")
    )
  )
)



#########################
## Define server logic ##
#########################
server <- function(input, output, session) {
 # General output
 ################
  # giving gene choices only when ViolinPlot or UMAP_GeneExpression is chosen as input
  output$genes <- renderUI({
    if (input$plot == "ViolinPlot" || input$plot == "UMAP_GeneExpression") {
      selectInput(
        inputId = "gene",
        label = "Choose a gene:",
        choices = allfeatures,
      )
    }
  })
  # Putting the logo in the dashbordheader
  output$logo <- renderImage({
    list(
      src = file.path("PPOL_Handjes.png"),
      contentType = "image/png",
      width = 50,
      height = 50
    )
  }, deleteFile = FALSE)
  
################################################################################
 # WT_Data output
 ################
  # Make the plot based on the selected plot type
  observeEvent(input$plot, {
    if (!is.null(input$plot)) {
      if (input$plot == "UMAP_Cluster") {
        output$WT_diffplots <- renderPlot({
          UMAPPlot(cell, group.by ="CellType")
        })
      } else {
        output$WT_diffplots <- renderPlot({
          if (input$plot == "UMAP_GeneExpression") {
            UMAPPlot(cell, group.by = "CellType") +
              FeaturePlot(cell, features = input$gene, cols = c("lightgrey", "#FF6600", "#FF0000"))
          } else if (input$plot == "ViolinPlot") {
            UMAPPlot(cell, group.by = "CellType") +
              VlnPlot(cell, features = input$gene, pt.size = 0)
          }
        })
      }
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
