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
# install.packages("shinyjs", lib = lib_path)
library(shinyjs, lib.loc = lib_path) # to reset input file for signaturetab
# install.packages("openxlsx", lib = lib_path)
library(openxlsx, lib.loc = lib_path)

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
                  # add icon with paper link
                  tags$li(class ="dropdown",
                          tags$a(href = "https://www.sciencedirect.com/science/article/pii/S2589004223021739",
                                 icon("file"),
                                 "Paper"
                                )
                          ),
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
                                   label = "Choose a plot:",
                                   choices = c("UMAP_Cluster", "UMAP_GeneExpression", "ViolinPlot"),
                                   selected = "UMAP_Cluster"
                                 ),
                                 uiOutput("genes"),
                                 tags$br(),
                                 uiOutput("statistic_choice"),
                                 tags$br(),
                                 div(style= "text-align: center;",
                                     uiOutput("download_per_plot"))
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
                                 ),
                                 tags$br(),
                                 div(style = "text-align: center;",
                                     downloadButton("download_degenes", "Download"))
                ),
                # third menuItem
                useShinyjs(),  # Initialize shinyjs
                menuItem("Signature Score", tabName = "score"),
                conditionalPanel("input.sidebarid == 'score'",
                                 fileInput("file", "Upload text file", accept = c(".txt")),
                                 textInput("genes", "Enter gene names (separated by comma and at least 5 genes)"),
                                 textInput("header", "Enter you plot title here:"),
                                 actionButton("calculate", "Calculate Signature Score"),
                                 tags$br(),
                                 div(style = "text-align: center;",
                                     downloadButton("download_signature","Download"))
                ),
                # fourth menuItem
                menuItem("Pseudotime analyse", tabName = "pseudo"),
                conditionalPanel("input.sidebarid == 'pseudo'",
                                 selectizeInput(
                                   inputId = "pseudogene",
                                   label = "Choose a gene:",
                                   choices = WT_features,
                                   selected = "SOX10"
                                 ),
                                 tags$br(),
                                 div(style= "text-align: center;",
                                     downloadButton("download_pseudo_dimplot", "Download the first plotrow")),
                                 tags$br(),
                                 div(style= "text-align: center;",
                                     downloadButton("download_pseudo_ggplots", "Download the second plotrow"))
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
                     uiOutput("WT_plottext"),
                     plotOutput("WT_diffplots", height = 500) # outputting the plots
                     #verbatimTextOutput("statistics")
            ), # outputting statistics for only the violinplot
            # Output tab for ALK_Data
            tabPanel(value = "gene_exp_ALK",
                     title = "ALK",
                     uiOutput("ALK_plottext"),
                     plotOutput("ALK_diffplots", height = 500))
          )
        )
      ),
      # output for second menu: differential gene analyse
      tabItem(
        tabName = "diff_gene",
            tabPanel(title = "WT",
                     tags$h3("Differential Gene Analysis"),
                     uiOutput("Datatable_text"),
                     tags$br(),
                     DT::dataTableOutput("DEGenes_table"))
      ),
      # output for third menu: signature score
      tabItem(tabName = "score",
              tabPanel(title = "WT",
                       tags$h3("Signature Score"),
                       uiOutput("signature_text"),
                       tags$br(),
                       plotOutput("signature_plot")
              )),
      # output for fourth menu: pseudotime analyse
      tabItem(tabName = "pseudo",
              tabPanel(title = "WT",
                       tags$h3("Pseudotime Analysis"),
                       uiOutput("pseudo_text"),
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
        "STMN2"
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
  
  # giving choices to get a comparison of two celltypes for WT
  output$statistic_choice <- renderUI({
    if (input$plot == "ViolinPlot") {
      if (input$gene_exp_tabs == "gene_exp_WT") {
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
    }
  })
  
  # Download buttons for the plots per tab per plot
  output$download_per_plot <- renderUI({
    if (input$gene_exp_tabs == "gene_exp_WT") {
      if (input$plot == "UMAP_Cluster") {
        downloadButton("download_UMAP_WT", "Download UMAP")
      }
      else if (input$plot == "UMAP_GeneExpression") {
        downloadButton("download_GeneExp_WT", "Download GeneExpression")
      }
      else if (input$plot == "ViolinPlot") {
        downloadButton("download_violinplot_WT", "Download ViolinPlot")
      }
    }
    else if (input$gene_exp_tabs == "gene_exp_ALK") {
      if (input$plot == "UMAP_Cluster") {
        downloadButton("download_UMAP_ALK", "Download UMAP")
      }
      else if (input$plot == "UMAP_GeneExpression") {
        downloadButton("download_GeneExp_ALK", "Download GeneExpression")
      }
    }
  })


  ##############################################################################
  # WT_Data output
  ################
  # making a text box above the plots from the first tab for WT
  output$WT_plottext <- renderUI({
    text <- switch(input$plot,
                   "UMAP_Cluster" = HTML(
                     "<h3>UMAP Cluster - Wild Type</h3>",
                     "<strong>Title:</strong> UMAP showing the different populations of interest.<br>",
                     "<strong>What does it show:</strong> A UMAP (Uniform Manifold Approximation and Projection) 
                     helps to visualize the data information and shows patterns in the data. 
                     This way, the dots that form clusters on the UMAP can be interpreted as separate cell types.<br>",
                     "<strong>Legend:</strong> Each dot represents a single cell, here colored by the different populations of interest.
                     The X and Y axes represent the two dimensions in the reduced space created by UMAP, which do not have
                     specific units but capture the structure and relationships within the data. Groupings of dots that are
                     close to each other represent similar data points, indicating potential clusters or subgroups within the dataset."
                   ),
                   "UMAP_GeneExpression" = HTML(
                       "<h3>UMAP Gene Expression - Wild Type</h3>",
                       "<strong>Title:</strong> UMAP showing the expression of a gene in the populations of interest.<br>",
                       "<strong>What does is show:</strong> The expression of a gene of interest in each single cell in the data.<br>",
                       "<strong>Legend:</strong> The gene expression is indicated by the red color. A darker red color means a higher expression in a particular cell type, 
                       while lighter red indicates lower or non-significant expression."
                   ),
                   "ViolinPlot" = HTML(
                     "<h3>Violin Plot - Wild Type</h3>",
                     "<strong>Title:</strong> Violin plot showing the expression of a gene in the different populations of interest.<br>",
                     "<strong>What does it show:</strong> The expression of a gene of interest in each cell type (population).<br>",
                     "<strong>Legend:</strong> The width of each violin represents the density of cells that express a particular gene at different expression levels.
                     Wider sections indicate a higher cell density at these expression levels.<br>",
                     "Statistics is based on the t test. The t-test compares the mean expression levels of a gene between two groups of cells.
                     The statistics are visualised using a line between different cell types. Above this line there is 'ns' or '*'.
                     If the line connecting two cell types is labeled with 'ns', it indicates that there is no significant statistical difference in gene expression 
                     levels between those cell types. On the other hand, if the line is labeled with '*', it indicates a significant statistical difference. 
                     The number of '*' symbols corresponds to the level of significance, with more stars suggesting a higher level of statistical difference."
                   )
    )
    text
  })
  
  # Reactive expression to create the violin plot
  WT_ViolinPlot <- reactive({
    req(input$gene, input$comparison1, input$comparison2)  # Ensure inputs are available
    plot <- VlnPlot(object = cell, 
                    features = input$gene, pt.size=0) +
      theme_classic() + 
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
            axis.title.y = element_blank())
    
    # Calculate the maximum y-value of the violin plot
    max_y <- max(ggplot_build(plot)$data[[1]]$ymax)
    
    # Add stat_compare_means at the top of the plot
    plot + stat_compare_means(comparisons = list(c(input$comparison1, input$comparison2)), 
                              label = "p.signif", 
                              label.y = max_y - 1)  # Adjust the value to position the label
                              # place for the statistic calculation output
                              #stat_compare_means(label.y = 1.5)
    })
  
  
  # outputting the plot based on the selected plot type
  output$WT_diffplots <- renderPlot({
    # give the UMAP_cluster plot
    if (input$plot == "UMAP_Cluster") {
      umap_plot_WT <- UMAPPlot(cell, group.by = "CellType")
      blank_plot <- ggplot() + 
        theme_void() + 
        theme(
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
        )
      grid.arrange(umap_plot_WT, blank_plot, ncol = 2, widths = c(2, 1))
    }
    # give the UMAP_cluster + UMAP_GeneExpression plots
    else if (input$plot == "UMAP_GeneExpression") {
      UMAPPlot(cell, group.by = "CellType") +
        FeaturePlot(cell, features = input$gene, cols = c("lightgrey", "#FF6600", "#FF0000"))
    }
    # Give the UMAP_cluster + violinplot, without the dots
    else if (input$plot == "ViolinPlot") {
      # give the UMAP_cluster + violinplot
      UMAPPlot(cell, group.by = "CellType") +
        WT_ViolinPlot()
    }
  })
  
  # Downloading the different gene expression plots for WT
  output$download_UMAP_WT <- downloadHandler(
    filename = function() {"UMAP_WT.png"},
    content = function(file) {
      png(file, width = 1250, height = 750)
      print(UMAPPlot(cell, group.by = "CellType"))
      dev.off()
    }
  )
  output$download_GeneExp_WT <- downloadHandler(
    filename = function() {paste(input$gene, "GeneExpression_WT.png", sep = "_")},
    content = function(file) {
      png(file, width = 1250, height = 750)
      print(FeaturePlot(cell, features = input$gene, cols = c("lightgrey", "#FF6600", "#FF0000")))
      dev.off()
    }
  )
  output$download_violinplot_WT <- downloadHandler(
    filename = function() {paste(input$gene, "ViolinPlot.png", sep = "_")},
    content = function(file) {
      png(file, width = 1250, height = 750)
      print(WT_ViolinPlot())
      dev.off()
    }
  )
  
  
  # making a textbox above the datatable
  output$Datatable_text <- renderUI({
    HTML("<strong>Title:</strong> Differential gene expression analysis.<br>",
         "<strong>What does is show:</strong> List of differential expressed genes in the selected population vs all other populations.<br>",
         "<strong>Legend:</strong><br>",
         "<ul>
         <li>The first column represents the <strong>gene symbols</strong>.</li>
         <li>The second column contains the <strong>p-values</strong>. This value comes from the statistical test used to determine if there is a significant 
         difference in gene expression between two groups. In this case group one is the selected variable and group two are all the other groups combined. 
         A p-value lower than 0.05 indicates a highly significant result.</li>
         <li>The third column is the <strong>avg_log2FC</strong> or average log2 fold change. A positive value indicates that the gene is upregulated in the first group compared 
         to the second group. A negative value indicates downregulation.</li>
         <li>The fourth column is the <strong>pct.1</strong> and represents the percentage of cells in the first group that express the gene above a certain threshold. 
         It shows how prevalent the gene expression is within the first group.</li>
         <li>The fifth column is the <strong>pct.2</strong> and is similar to the fourth column. It represents the percentage of cells in the second group that express the gene above a certain threshold.</li>
         <li>The last column is the <strong>p_val_adj</strong> or the adjusted p-value. This column contains the p-values for multiple testing. Because many genes are tested simultaneously, 
         adjustments are made to control the false discovery rate(FDR). This reduces the likelihood of type I errors (false positives).</li>
         </ul>"
         )
  })
  
  # The WT_Data table from the DEGenes.RDS
  output$DEGenes_table <- DT::renderDT({
    req(input$selected_var)
    DEGenes_sample <- DEGenes[[input$selected_var]]
    datatable <- datatable(
      data = DEGenes_sample,
      options = list(pageLength = 50))
    datatable
  })
  
  # Download DEGenes table in xlsx
  output$download_degenes <- downloadHandler(
    filename = function() {"DeGenes.xlsx"},
    content = function(file) {
      wb <- createWorkbook()
      for (var in names(DEGenes)) {
        DEGenes_sample <- DEGenes[[var]]
        DEGenes_sample <- data.frame(Gene = rownames(DEGenes_sample), DEGenes_sample)
        addWorksheet(wb, var)
        writeData(wb, var, DEGenes_sample)
      }
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # making text for the signatureplot tab
  output$signature_text <- renderUI({
    HTML("<strong>Title:</strong> UMAP showing the signature score of the gene set in the populations of interest.<br>",
         "<strong>What does it show:</strong> It shows the combined expression of all genes in a gene set, the signature score calculations are based on the UCell method and enrichIt function.<br>",
         "<strong>Legend:</strong> The signature score is indicated by the red color. A darker red color means a higher expression in a particular cell, while lighter red indicates lower or non-significant expression.<br>",
         "<br>",
         "<strong>Please be patient, this plot can take a while to calculate.</strong>"
         )
  })
  
  # Reactive expression to read and process the uploaded file
  genes_from_file <- reactive({
    req(input$file)
    file_ext <- tools::file_ext(input$file$name)
    if (file_ext != "txt") {
      showNotification("Please upload a text file with .txt extension", type = "warning")
      return(NULL)
    }
    df <- read.table(input$file$datapath, header = TRUE)
    gene_column <- which(!is.na(df[[1]]))  # Assuming gene names are in the first column
    if (length(gene_column) == 0) {
      showNotification("The uploaded file does not contain gene names", type = "warning")
      return(NULL)
    }
    trimws(df[[1]])  # Trim leading and trailing whitespace
  })
  
  observeEvent(input$calculate, {
    # Check if genes are provided either by text input or file upload
    genes <- NULL
    if (input$genes != "") {
      genes <- strsplit(input$genes, ",")[[1]]
      genes <- trimws(genes)  # Trim leading and trailing whitespace
    } else if (!is.null(input$file)) {
      genes <- genes_from_file()
    }
    
    # Check if genes are found
    if (is.null(genes) || length(genes) == 0) {
      showNotification("No genes found", type = "warning")
      return(NULL)
    }
    
    # Define signatures
    signatures <- list(user_genes = genes)
    
    # Get header input
    header <- input$header
    
    # Show notification for signature score calculation
    showNotification("Signature Score being calculated, please be patient")
    # Run signature scoring
    DefaultAssay(cell) <- "RNA"
    u.scores <- enrichIt(obj = cell, gene.sets = signatures, groups = 2000, 
                         cores = 1, method = "UCell")
    
    # Add scores to Seurat object
    u.scores <- as.data.frame(u.scores)
    cell <- AddMetaData(cell, u.scores)
    
    # Reset the file input
    reset("file")
    
    # Feature plot
    featureplot <- FeaturePlot(cell, features = "user_genes", cols = c("lightgrey", "#FF6600", "#FF0000"))
    Feature_HeaderPlot<- featureplot + labs(title = header)
    
    # Render the plot with header
    output$signature_plot <- renderPlot({
      UMAPPlot(cell, group.by = "CellType") + Feature_HeaderPlot
    })
    
    # downloadbutton for the signatureplot
    output$download_signature <- downloadHandler(
      filename = function() {
        if (header != "") {
          paste(header, ".png", sep = "")
        } else {
          "signatureplot.png"
        }
      },
      content = function(file) {
        png(file, width = 1250, height = 750)
        # Print the plot with header
        print(Feature_HeaderPlot)
        dev.off()
      }
    )
  })
  
  # making information text for the pseudoanalysis tab
  output$pseudo_text <- renderUI({
      HTML(
        "<strong>Title:</strong> Lineage trajectories<br>",
        "<strong>What does it show:</strong><br>",
        "<ol type='A'>
        <li>Diffusion plot showing the differentiating cells colored and annotated by cell population.</li>
        <li>Diffusion maps of in vitro differentiating cells colored by pseudotime trajectory for the two different differentiation arms.</li>
        <li>Pseudotime analysis plot showing the expression of a gene along the lineage trajectory in the diffusion plot, annotated for the populations of interest, again for the two different differentiation arms.</li>
        </ol>",
        "<strong>Legend:</strong> The diffusion map embeddings were calculated using the destiny R package. Single-cell pseudotime trajectories were constructed using the slingshot R package based on the 
        populations of interest and the DC components from destiny. A node in cluster 12 (SCPs) was selected as the starting point for the trajectory."
      )
  })
  
  # Making the dimplot for the pseudotime analysis
  Dimplot <- DimPlot(cell, reduction = "DC") + NoLegend() + ggtitle("A")
  # output dimplot for the pseudotime analysis
  output$WT_dimplot <- renderPlot({
    grid.arrange(Dimplot, SCP_SAP_PLOT, SCP_PROSAP_PLOT, ncol = 3)
  })
  
  # Reactive expression to generate ggplots
  ggplots <- reactive({
    req(input$pseudogene)  # Ensure input$pseudogene is available
    
    ggplot1 <- qplot(sce_slingshot$slingPseudotime_1, 
                     as.numeric(cell@assays$RNA@data[input$pseudogene,]), 
                     color = sce_slingshot$CellType, 
                     ylab = paste("Expression of", input$pseudogene, sep = " "), 
                     xlab = "Pseudotime") + 
      theme_bw() + 
      geom_smooth(aes(group = 1), se = FALSE, method = "loess", color = "black") +
      theme(legend.position = "none") +
      ggtitle("From SCP to Sympathoblasts") +
      theme(plot.title = element_text(face = "bold"))
    
    ggplot2 <- qplot(sce_slingshot$slingPseudotime_2, 
                     as.numeric(cell@assays$RNA@data[input$pseudogene,]), 
                     color = sce_slingshot$CellType, 
                     ylab = paste("Expression of", input$pseudogene, sep = " "), 
                     xlab = "Pseudotime") + 
      theme_bw() + 
      geom_smooth(aes(group = 1), se = FALSE, method = "loess", color = "black") +
      theme(legend.position = "none") +
      ggtitle("From SCP to ProlifSympathoblasts") +
      theme(plot.title = element_text(face = "bold"))
    
    # Arrange the plots next to each other
    grid.arrange(ggplot1, ggplot2, ncol = 2)
  })
  
  # Output pseudoplots
  output$WT_pseudoplots <- renderPlot({
    ggplots()
  })
  
  # Download the first row of plots (the dimplot)
  output$download_pseudo_dimplot <- downloadHandler(
    filename = function() { "Overall_PseudotimeAnalyse.png" },
    content = function(file) {
      png(file, width = 1750, height = 750)
      grid.arrange(Dimplot, SCP_SAP_PLOT, SCP_PROSAP_PLOT, ncol = 3)
      dev.off()
    }
  )
  
  # Download the second row of plots (the ggplots)
  output$download_pseudo_ggplots <- downloadHandler(
    filename = function() {paste(input$pseudogene, "PseudotimeAnalysis.png", sep = "_") },
    content = function(file) {
      ggplots <- ggplots()
      ggsave(file, plot = ggplots, device = "png", width = 12, height = 6)  # Adjust width and height as needed, here specified in inches, not pixels
    }
  )
  
  ##############################################################################
  # ALK_Data output
  #################
  # making a textbox above the plots from the first tab of ALK
  output$ALK_plottext <- renderUI({
    text <- switch(input$plot,
                   "UMAP_Cluster" = HTML(
                     "<h3>UMAP Cluster - ALK</h3>",
                     "<strong>Title:</strong> UMAP showing the different populations of interest.<br>",
                     "<strong>What does it show:</strong> A UMAP (Uniform Manifold Approximation and Projection) 
                     helps to visualize the data information and shows patterns in the data. 
                     This way, the dots that form clusters on the UMAP can be interpreted as separate cell types.<br>",
                     "<strong>Legend:</strong> Each dot represents a single cell, here colored by the different populations of interest.
                     The X and Y axes represent the two dimensions in the reduced space created by UMAP, which do not have
                     specific units but capture the structure and relationships within the data. Groupings of dots that are
                     close to each other represent similar data points, indicating potential clusters or subgroups within the dataset."
                   ),
                   "UMAP_GeneExpression" = HTML(
                     "<h3>UMAP Gene Expression - ALK</h3>",
                     "<strong>Title:</strong> UMAP showing the expression of a gene in the populations of interest.<br>",
                     "<strong>What does is show:</strong> The expression of a gene of interest in each single cell in the data.<br>",
                     "<strong>Legend:</strong> The gene expression is indicated by the red color. A darker red color means a higher expression in a particular cell type, 
                       while lighter red indicates lower or non-significant expression."
                   ),
                   "ViolinPlot" = HTML(
                     "<h3>Violin Plot - ALK</h3>",
                     "Violin plots were not created for the ALK data, for more information please use the contact button at the top and we are happy to discuss this."
                   )
    )
    text
  })
  
  output$ALK_diffplots <- renderPlot({
    # give the UMAP_cluster plot
    if (input$plot == "UMAP_Cluster") {
      umap_plot_ALK <- UMAPPlot(ALK, group.by = "Celltype")
      blank_plot <- ggplot() + 
        theme_void() + 
        theme(
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
        )
      grid.arrange(umap_plot_ALK, blank_plot, ncol = 2, widths = c(2, 1))
    }
    # give the UMAP_cluster + UMAP_GeneExpression plots
    else if (input$plot == "UMAP_GeneExpression") {
      UMAPPlot(ALK, group.by = "Celltype") +
        FeaturePlot(ALK, features = input$gene, cols = c("lightgrey", "#FF6600", "#FF0000"))
    }
    # No Violinplot for the ALK data
  })
  
  # Downloading the different gene expression plots for ALK
  output$download_UMAP_ALK <- downloadHandler(
    filename = function() {"UMAP_ALK.png"},
    content = function(file) {
      png(file, width = 1250, height = 750)
      print(UMAPPlot(ALK, group.by = "Celltype"))
      dev.off()
    }
  )
  output$download_GeneExp_ALK <- downloadHandler(
    filename = function() {paste(input$gene, "GeneExpression_ALK.png", sep = "_")},
    content = function(file) {
      png(file, width = 1250, height = 750)
      print(FeaturePlot(ALK, features = input$gene, cols = c("lightgrey", "#FF6600", "#FF0000")))
      dev.off()
    }
  )
  
  
  
}


#########################
## Run the application ##
#########################
shinyApp(ui = ui, server = server)