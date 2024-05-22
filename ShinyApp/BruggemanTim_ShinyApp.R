### ShinyApp for single cell RNAseq data ###
############################################

######################
## Packages to load ##
######################
# making a library to store my packages in with path
lib_path <- "/user/gent/476/vsc47620/R/x86_64-pc-linux-gnu-library/4.2"

# install.packages("shiny",lib = lib_path)
library(shiny, lib.loc = lib_path)
# install.packages("Seurat", lib = lib_path)
library(Seurat, lib.loc = lib_path)
# install.packages("bslib", lib = lib_path)
library(bslib, lib.loc = lib_path)
# install.packages("shinydashboard", lib = lib_path)
library(shinydashboard, lib.loc = lib_path)
# install.packages("DT", lib = lib_path)
library(DT, lib.loc = lib_path)
# install.packages("dplyr", lib = lib_path)
library(dplyr, lib.loc = lib_path)
# install.packages("ggplot2", lib = lib_path)
library(ggplot2, lib.loc = lib_path)
# install.packages("ggpubr", lib = lib_path)
library(ggpubr, lib.loc = lib_path) # for the stat_compare_means
# BiocManager::install("slingshot", lib = lib_path)
library(slingshot, lib.loc = lib_path)
# install.packages("viridis", lib = lib_path)
library(viridis, lib.loc = lib_path)
# BiocManager::install("destiny", lib = lib_path)
library(destiny, lib.loc = lib_path) # for the diffusionmap
# install.packages("escape", lib = lib_path)
library(escape, lib.loc = lib_path)
# install.packages("UCell", lib = lib_path)
library(UCell, lib.loc = lib_path) # for the signature score
# install.packages("gridExtra", lib = lib_path)
library(gridExtra, lib.loc = lib_path)

options(bitmapType = "cairo") # specific for the HPC to make plots


#########################
## Loading in the data ##
#########################
#annot <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/RDSObjects/ada_and_stefNoALK_dim20_res06_ANNOT.rds")
cell <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/RDSObjects/CellsOfInterest_SLB.rds")
ALK <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/RDSObjects/stef_ALK_and_noALK_integrated_SLB.rds")
DEGenes <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/DEGenes.RDS")
#ALK_DEGenes <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/DEGenes_ALK.RDS")

#######################################
## variables to use in the interface ##
#######################################
## getting only the RNA data for the assays
DefaultAssay(cell) <- "RNA"
#DefaultAssay(ALK) <- "RNA"

## features for the genes to chose from
WT_features <- rownames(cell)
ALK_features <- rownames(ALK)
#allfeatures <- c(WT_features, ALK_features)

## source the preparation of the pseudotime analysis
source("Pseudotime_Calculations.R", local = TRUE)

###############
## Define UI ##
###############
ui <- dashboardPage(
  # headerbar
  dashboardHeader(title = 'ShinyApp',
                  # Add icon with email link
                  tags$li(class = "dropdown",
                          tags$a(href = "mailto:sarahlee.bekaert@ugent.be",
                                 icon("envelope"),
                                 "Contact")
                  ),
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
                                   choices = c("UMAP_Cluster", "UMAP_GeneExpression", "ViolinPlot"),
                                   selected = "UMAP_Cluster"
                                 ),
                                 uiOutput("genes"),
                                 tags$br(),
                                 uiOutput("statistic_choice")
                ),
                # Second menuItem
                menuItem("DE Genes", tabName = "diff_gene"),
                # Conditional panel for the WT tab under DE Genes
                conditionalPanel("input.sidebarid == 'diff_gene'",
                                 radioButtons(
                                   inputId = "selected_var",
                                   label = "Select variables:",
                                   choices = names(DEGenes),
                                   selected = "Sympathoblasts"
                                 )
                ),
                # third menuItem
                menuItem("Signature Score", tabName = "score"),
                conditionalPanel("input.sidebarid == 'score'",
                                 fileInput("file", "Upload text file",
                                           accept = c(".txt")),
                                 textInput("genes", "Enter gene names (separated by comma)"),
                                 actionButton("calculate", "Calculate Signature Score")
                ),
                # fourth menuItem
                menuItem("Pseudotime analyse", tabName = "pseudo"),
                conditionalPanel("input.sidebarid == 'pseudo'",
                                 selectizeInput(
                                   inputId = "pseudogene",
                                   label = "Choose a gene:",
                                   choices = WT_features
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
          tabsetPanel(
            id = "gene_exp_tabs",
            # Output tab for WT_Data
            tabPanel(value = "gene_exp_WT",
                     title = "WT",
                     tags$br(),
                     uiOutput("WT_plottext"),
                     plotOutput("WT_diffplots", height = 500) # outputting the plots
                     #verbatimTextOutput("statistics")
            ), # outputting statistics for only the violinplot
            # Output tab for ALK_Data
            tabPanel(value = "gene_exp_ALK",
                     title = "ALK",
                     uiOutput("ALK_plottext"),
                     tags$br(),
                     plotOutput("ALK_diffplots", height = 500))
          )
        )
      ),
      # output for second menu: differential gene analyse
      tabItem(
        tabName = "diff_gene",
            tabPanel(title = "WT",
                     textOutput("Datatable_text"),
                     tags$br(),
                     DT::dataTableOutput("DEGenes_table"))
      ),
      # output for third menu: signature score
      tabItem(tabName = "score",
              tabPanel(title = "WT",
                       textOutput("signature_text"),
                       tags$br(),
                       plotOutput("signature_plot", height = 500, width = 800)
              )),
      # output for fourth menu: pseudotime analyse
      tabItem(tabName = "pseudo",
              tabPanel(title = "WT",
                       textOutput("pseudo_text"),
                       tags$br(),
                       plotOutput("WT_dimplot"),
                       plotOutput("WT_pseudoplots"))
      )
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
  
  # Reactive value to hold gene choices based on the selected tab: WT vs ALK
  gene_choices <- reactive({
    if (input$gene_exp_tabs == "gene_exp_WT") {
      WT_features
    } else if (input$gene_exp_tabs == "gene_exp_ALK") {
      ALK_features
    } else {
      NULL
    }
  })
  
  # Render gene selection input based on the plot and tab selection
  output$genes <- renderUI({
    if (input$plot == "ViolinPlot" || input$plot == "UMAP_GeneExpression") {
      selected_gene_choice <- if (input$gene_exp_tabs == "gene_exp_WT") {
        "MRPL20"
      } else if (input$gene_exp_tabs == "gene_exp_ALK") {
        "ALK" #ALK gen nemen
      } else {
        NULL
      }
      
      selectizeInput(
        inputId = "gene",
        label = "Choose a gene:",
        choices = gene_choices(),
        selected = selected_gene_choice
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
  # making a text box above the plots from the first tab for WT
  output$WT_plottext <- renderUI({
    text <- switch(input$plot,
                   "UMAP_Cluster" = "UMAP Cluster Plot for WT Data",
                   "UMAP_GeneExpression" = "UMAP Gene Expression Plot for WT Data",
                   "ViolinPlot" = "Violin Plot for WT Data")
    text
  })
  
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
      # make the violinplot with the selected gene
      WT_ViolinPlot <- VlnPlot(object = cell, 
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
                           label.y = 1)  #+ 
        # place for the statistic calculation output
        #stat_compare_means(label.y = 1.5)
      
      # give the UMAP_cluster + violinplot
      UMAPPlot(cell, group.by = "CellType") +
        WT_ViolinPlot
    }
  })
  
  # making a textbox above the datatable
  output$Datatable_text <- renderText({
    paste("This is a datatable where you can find the genes per WT-celltype. In this table you can also find calculations from the genes.")
  })
  
  # The WT_Data table from the DEGenes.RDS
  output$DEGenes_table <- DT::renderDT({
    req(input$selected_var)
    DEGenes_sample <- DEGenes[[input$selected_var]]
    datatable(
      data = DEGenes_sample,
      options = list(pageLength = 50)
    )
  })
  
  # making text for the signatureplot tab
  output$signature_text <- renderText({
    paste("In this tab a list of genes can be inputted, this list will be calculated and a featureplot will be given")
  })
  
  # making the signature tab output
  observeEvent(input$file, {
    # Check if file is uploaded
    req(input$file$datapath)
    # Check file extension
    file_ext <- tools::file_ext(input$file$name)
    if (file_ext != "txt") {
      showNotification("Please upload a text file with .txt extension", type = "warning")
      return(NULL)
    }
  })
  observeEvent(input$calculate, {
    # Check if genes are provided
    if (input$genes == "" && is.null(input$file)) {
      showNotification("Please enter gene names or upload a text file", type = "warning")
      return(NULL)
    }
    # Extract genes
    if (!is.null(input$file)) {
      df <- read.table(input$file$datapath, header = TRUE)  # Read the entire file
      gene_column <- which(!is.na(df))  # Find non-NA column (assuming gene names are in the first column)
      if (length(gene_column) == 0) {
        showNotification("The uploaded file does not contain gene names", type = "warning")
        return(NULL)
      }
      genes <- trimws(df[[1]])  # Trim leading and trailing whitespace --> otherwise gets " CD24" instead of "CD24"
    } else {
      genes <- strsplit(input$genes, ",")[[1]]
      genes <- trimws(genes)  # Trim leading and trailing whitespace
    }
    
    # Check if genes are provided
    if (length(genes) == 0) {
      showNotification("No genes found", type = "warning")
      return(NULL)
    }
    
    # Define signatures
    signatures <- list(user_genes = genes)
    
    # Show notification for signature score calculation
    showNotification("Signature Score being calculated, please be patient")
    
    # Run signature scoring
    DefaultAssay(cell) <- "RNA"
    u.scores <- enrichIt(obj = cell, gene.sets = signatures, groups = 2000, 
                         cores = 1, method = "UCell")
    
    # Add scores to Seurat object
    u.scores <- as.data.frame(u.scores)
    cell <- AddMetaData(cell, u.scores)
    
    # Feature plot
    output$signature_plot <- renderPlot({
      FeaturePlot(cell, features = "user_genes") ########## to do: asking for header and putting that header on it then............
    })
  })
  
  output$pseudo_text <- renderText({
    paste("This tab gives the pseudotime analysis for the WT")
  })
  
  # output dimplot for the pseudotime analysis
  output$WT_dimplot <- renderPlot({
    DimPlot(cell, reduction = "DC")
  })
  
  # output pseudoplots
  output$WT_pseudoplots <- renderPlot({

  # making the pseudotimeplot from SCP to Sympathoblasts
    ggplot1 <- qplot(sce_slingshot$slingPseudotime_1, as.numeric(cell@assays$RNA@data[input$pseudogene,]), 
                     color = sce_slingshot$CellType, ylab = paste("Expression of", input$pseudogene, sep = " "), xlab = "Pseudotime") + 
      theme_bw() + 
      geom_smooth(aes(group = 1), se = FALSE, method = "loess", color = "gray") +
      theme(legend.position = "none") +
      ggtitle("From SCP to Sympathoblasts")

    # Making the pseudotimeplot from SCP to ProlifSympathoblasts
    ggplot2 <- qplot(sce_slingshot$slingPseudotime_2, as.numeric(cell@assays$RNA@data[input$pseudogene,]), 
                     color = sce_slingshot$CellType, ylab = paste("Expression of", input$pseudogene, sep = " "), xlab = "Pseudotime") + 
      theme_bw() + 
      geom_smooth(aes(group = 1), se = FALSE, method = "loess", color = "gray") +
      theme(legend.position = "none") +
      ggtitle("From SCP to ProlifSympathoblasts") 
    
    # outputting the two pseudotimeplots
    grid.arrange(ggplot1, ggplot2, ncol = 2)
  })
  
  
  ##############################################################################
  # ALK_Data output
  #################
  # making a textbox above the plots from the first tab of ALK
  output$ALK_plottext <- renderUI({
    text <- switch(input$plot,
                   "UMAP_Cluster" = "UMAP Cluster Plot for ALK Data",
                   "UMAP_GeneExpression" = "UMAP Gene Expression Plot for ALK Data",
                   "ViolinPlot" = "Violin Plot for ALK Data")
    text
  })
  
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
      # make the violinplot with the selected gene
      ALK_ViolinPlot <- VlnPlot(object = ALK, 
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
        stat_compare_means(comparisons = list("ALK","SC"),label = "p.signif",
                           label.y = 1)  + 
        # place for the statistic calculation output
        stat_compare_means(label.y = 1.5)
      
      # give the UMAP_cluster + violinplot
      UMAPPlot(ALK, group.by = "Celltype") +
        ALK_ViolinPlot
    }
  })
  
  
}


#########################
## Run the application ##
#########################
shinyApp(ui = ui, server = server)