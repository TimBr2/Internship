### ShinyApp for single cell RNAseq data ###
############################################

######################
## Packages to load ##
######################
library(shiny)
library(Seurat)
library(bslib)
library(shinydashboard) #install.packages("shinydashboard")
library(DT)
library(dplyr)
library(ggplot2)
library(ggpubr) # for the stat_compare_means
library(slingshot)
library(viridis)
library(destiny) # for the diffusionmap

options(bitmapType = "cairo") # specific for the HPC to make plots


#########################
## Loading in the data ##
#########################
#annot <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/RDSObjects/ada_and_stefNoALK_dim20_res06_ANNOT.rds")
cell <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/RDSObjects/CellsOfInterest_SLB.rds")
ALK <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/RDSObjects/stef_ALK_and_noALK_integrated_SLB.rds")
DEGenes <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/DEGenes.RDS")

#######################################
## variables to use in the interface ##
#######################################
## getting only the RNA data for the assays
DefaultAssay(cell) <- "RNA"
#DefaultAssay(ALK) <- "RNA"

## features for the genes to chose from
WT_features <- rownames(cell)
ALK_features <- rownames(ALK)
allfeatures <- c(WT_features, ALK_features)

## preparation of the pseudotime analysis
# make a matrix
sce1 <- as.SingleCellExperiment(cell)
# Dim 18
dm_S <- DiffusionMap(Embeddings(cell, "pca")[,1:18])
cellLabels <- sce1$CellType
tmp <- data.frame(DC1 = eigenvectors(dm_S)[, 1],
                  DC2 = eigenvectors(dm_S)[, 2],
                  DC3 = eigenvectors(dm_S)[, 3],
                  DC4 = eigenvectors(dm_S)[, 4],
                  Samples = cellLabels)
dimentions_DC <- cbind(dm_S$DC1, dm_S$DC2)
colnames(dimentions_DC) <- c("DM_1", "DM_2")
cell[["DC"]] <- CreateDimReducObject(embeddings = dimentions_DC, key = "DC_", assay = DefaultAssay(cell))
# slingshot function with DC as reduction 
# start from SingleCellExperiment matrix  
sce <- as.SingleCellExperiment(cell)
sce_slingshot <- slingshot(sce, clusterLabels = "CellType", reducedDim = "DC")


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
                menuItem(text = "Gene Expression", tabName = "gene_exp"),
                # Choose between the three plots to visualize
                conditionalPanel("input.sidebarid == 'gene_exp'",
                                 selectInput(
                                   inputId = "plot",
                                   label = "choose a plot:",
                                   choices = c("UMAP_Cluster","UMAP_GeneExpression","ViolinPlot"),
                                   selected = "UMAP_Cluster"
                                 ),
                                 uiOutput("genes"),
                                 tags$br(),
                                 uiOutput("statistic_choice")
                ),
                # Second menuItem
                menuItem("DE Genes", tabName = "diff_gene"),
                # radiobuttions to select the datatable from the celltype
                conditionalPanel("input.sidebarid == 'diff_gene'",
                                 radioButtons(inputId = "selected_var",
                                   label = "Select variables:",
                                   choices = names(DEGenes),
                                   selected = "Sympathoblasts"
                                 )
                ),
                # third menuItem
                menuItem("Signature Score", tabName = "score"),
                # fourth menuItem
                menuItem("Pseudotime analyse", tabName = "pseudo"),
                conditionalPanel("input.sidebarid == 'pseudo'",
                  selectInput(
                    inputId = "pseudogene",
                    label = "Choose a gene:",
                    choices = allfeatures
                  )
                )
    )
  ),
  # Output in the body
  dashboardBody(
    tabItems(
      # output for first menu: the plots visualization
      tabItem(
        tabName = "gene_exp",
        navset_pill(
          # Output tab for WT_Data
          tabPanel(title = "WT", 
                   plotOutput("WT_diffplots", height = 500, width = 1250),# outputting the plots
                   plotOutput("WT_violin2", width = 1250),
                   verbatimTextOutput("statistics")), # outputting statistics for only the violinplot
          # Output tab for ALK_Data
          tabPanel(title = "ALK", 
                   plotOutput("ALK_diffplots"))
        )
      ),
      # output for second menu: differential gene analyse
      tabItem(
        tabName = "diff_gene",
        navset_pill(
          tabPanel(title = "WT",
                   DT::dataTableOutput("DEGenes_table")),
          tabPanel(title = "ALK")
        )
      ),
      # output for third menu: signature score
      tabItem(tabName = "score"),
      # output for fourth menu: pseudotime analyse
      tabItem(tabName = "pseudo",
              navset_pill(
                tabPanel(title = "WT",
                         plotOutput("WT_dimplot"),
                         plotOutput("WT_pseudoplots")),
                tabPanel(title = "ALK")
              ))
    )
  )
)


#########################
## Define server logic ##
#########################
server <- function(input, output, session) {
  # General output
  ################
  # Putting the logo in the dashbordheader
  output$logo <- renderImage({
    list(
      src = file.path("PPOL_Handjes.png"),
      contentType = "image/png",
      width = 50,
      height = 50
    )
  }, deleteFile = FALSE)
  
  # giving gene choices only when ViolinPlot or UMAP_GeneExpression is chosen as input
  output$genes <- renderUI({
    if (input$plot == "ViolinPlot" || input$plot == "UMAP_GeneExpression") {
      selectizeInput(
        inputId = "gene",
        label = "Choose a gene:",
        choices = allfeatures
      )
    }
  })
  
  # giving choices to get a comparison of two celltypes
  output$statistic_choice <- renderUI({
    if (input$plot == "ViolinPlot") {
        tags$div( # put in div to get both selectinputs in here
          tags$div(
            HTML("Options for statistics:"),
            style = "text-align: center;" # align content of div to center
            ),
          selectInput(
            inputId = "comparison1",
            label = "choose first comparison",
            choices = levels(cell$CellType),
            selected = "Sympathoblasts"
          ),
          tags$div(
            HTML("VS"),
            style = "text-align: center;" # Align the content of the div to the center, put in separate div to only get "VS" in center
          ),
          selectInput(
            inputId = "comparison2",
            label = "Choose second comparison",
            choices = levels(cell$CellType),
            selected = "Sympathoblasts"
          )
        )
    }
  })
  
  
  ##############################################################################
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
        FeaturePlot(cell, features = input$gene, cols = c("lightgrey", "#FF6600", "#FF0000"))
    }
    # Give the UMAP_cluster + violinplot, without the dots
    else if (input$plot == "ViolinPlot") {
      UMAPPlot(cell, group.by = "CellType") +
        VlnPlot(cell, features = input$gene, pt.size = 0)
    }
  })
  
  # violin2 output plot for statistics
  output$WT_violin2 <- renderPlot({
    if (input$plot == "ViolinPlot") {
    VlnPlot(object = cell, 
            features = input$gene, pt.size=0, 
            #cols = c("dodgerblue2", "chartreuse3", "indianred3", "darkorange", "mediumorchid")
            )+
      #coord_flip()+ 
      theme_classic()+ 
      theme(axis.text.x = element_text(angle = 45, 
                                       hjust = 1, 
                                       vjust = 1,
                                       size = 10), 
            axis.text.y = element_text(size = 10),  
            plot.title = element_text(size=10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position="none",
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) + 
      # This is for the statistic line
      stat_compare_means(comparisons = list(c(input$comparison1, c(input$comparison2))),label = "p.signif",
                         label.y = 1)  + 
      # place for the statistic calculation output
      stat_compare_means(label.y = 1.5) 
    }
  })

  
  # The WT_Data table from the DEGenes.RDS
  output$DEGenes_table <- renderDataTable({
    req(input$selected_var)
    DEGenes_sample <- DEGenes[[input$selected_var]]
    datatable(
      data = DEGenes_sample,
      options = list(pageLength = 50)
    )
  })
  
  # output dimplot for the pseudotime analysis
  output$WT_dimplot <- renderPlot({
    DimPlot(cell, reduction = "DC")
  })
  
  # output pseudoplots
  output$WT_pseudoplots <- renderPlot({
    # making the pseudotimeplot from SCP to Sympathoblasts
    qplot(sce_slingshot$slingPseudotime_1, as.numeric(cell@assays$RNA@data[input$pseudogene,]), 
          color = sce_slingshot$CellType, ylab = paste("Expression of", input$pseudogene, sep = " "), xlab = "Pseudotime") + 
      theme_bw() + 
      geom_smooth(aes(group = 1), se = FALSE, method = "loess", color = "gray") +
      
    # Making the pseudotimeplot from SCP to ProlifSympathoblasts
    qplot(sce_slingshot$slingPseudotime_2, as.numeric(cell@assays$RNA@data[input$pseudogene,]), 
          color = sce_slingshot$CellType, ylab = paste("Expression of", input$pseudogene, sep = " "), xlab = "Pseudotime") + 
    theme_bw() + 
    geom_smooth(aes(group = 1), se = FALSE, method = "loess", color = "gray") +
    theme(legend.position = "none")
  })
  
  
  ##############################################################################
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