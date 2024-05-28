### Pseudotime_Calculations.R to calculate the pseudotime and use as source for the ShinyApp ###
################################################################################################

# Loading in the data
cell <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/RDSObjects/CellsOfInterest_SLB.rds")


## preparation of the pseudotime analysis
# make a matrix
set.seed(1) # otherwise the plot is turned left to right
sce0 <- as.SingleCellExperiment(cell)
# Dim 18
dm_S <- DiffusionMap(Embeddings(cell, "pca")[,1:18])
cellLabels <- sce0$CellType
tmp <- data.frame(DC1 = eigenvectors(dm_S)[, 1],
                  DC2 = eigenvectors(dm_S)[, 2],
                  DC3 = eigenvectors(dm_S)[, 3],
                  DC4 = eigenvectors(dm_S)[, 4],
                  Samples = cellLabels)
dimentions_DC <- cbind(dm_S$DC1, dm_S$DC2)
colnames(dimentions_DC) <- c("DM_1", "DM_2")
cell[["DC"]] <- CreateDimReducObject(embeddings = dimentions_DC, key = "DC_", assay = DefaultAssay(cell))


## /2: Pseudotime from SCP to Sympathoblasts
############################################
SCP_SAP <- subset(cell, subset = seurat_clusters_old == c("12", "22", "11", "9"))

SCP_SAP <- as.SingleCellExperiment(SCP_SAP)
# slingshot function with DC as reduction 
SCP_SAP_slingshot <- slingshot(SCP_SAP, clusterLabels = "seurat_clusters_old", reducedDim = "DC")

# Make this a ggplot! 
df_sce_SCP_SAP <- as.data.frame(reducedDims(SCP_SAP_slingshot)$DC)

SCP_SAP_PLOT <- ggplot(df_sce_SCP_SAP, aes(x = DC_1, y = DC_2, color = SCP_SAP_slingshot$slingPseudotime_1)) + 
  geom_point() + 
  theme_classic() + 
  scale_color_viridis() + 
  ggtitle("Pseudotime analysis SCP to SAPs") + 
  geom_smooth(se = FALSE, method = "loess", color = "black") +
  labs(color = "Pseudotime")

## /3: Pseudotime from SCP to ProlifSympathoblasts
##################################################
SCP_pro_SAP <- subset(cell, subset = seurat_clusters_old == c("12", "22", "11", "18"))

SCP_pro_SAP <- as.SingleCellExperiment(SCP_pro_SAP)
# slingshot function with DC as reduction 
SCP_pro_SAP_slingshot <- slingshot(SCP_pro_SAP, clusterLabels = "seurat_clusters_old", reducedDim = "DC")

# Make this a ggplot! 
df_sce_SCP_pro_SAP <- as.data.frame(reducedDims(SCP_pro_SAP_slingshot)$DC)
# Reverse the pseudotime values
SCP_pro_SAP_slingshot$slingPseudotime_1_reverse <- max(SCP_pro_SAP_slingshot$slingPseudotime_1) - SCP_pro_SAP_slingshot$slingPseudotime_1

SCP_PROSAP_PLOT <- ggplot(df_sce_SCP_pro_SAP, aes(x = DC_1, y = DC_2, color = SCP_pro_SAP_slingshot$slingPseudotime_1)) + 
  geom_point() + 
  theme_classic() + 
  scale_color_viridis() + 
  ggtitle("Pseudotime analysis SCP to Proliferating SAPs") + 
  geom_smooth(se = FALSE, method = "loess", color = "black") +
  theme(legend.position = "none")


## /4: Pseudotime trajectorie per gene
######################################
# slingshot function with DC as reduction 
# start from SingleCellExperiment matrix  
sce <- as.SingleCellExperiment(cell)
sce_slingshot <- slingshot(sce, clusterLabels = "CellType", reducedDim = "DC")