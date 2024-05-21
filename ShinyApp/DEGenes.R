### Diff gene analyse, script to make an RDSObject from with the output ###
###########################################################################

## load in packages##
#####################
library(Seurat)


## Loading in the data ##
#########################
#annot <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/RDSObjects/ada_and_stefNoALK_dim20_res06_ANNOT.rds")
cell <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/RDSObjects/CellsOfInterest_SLB.rds")
ALK <- readRDS("/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/RDSObjects/stef_ALK_and_noALK_integrated_SLB.rds")


## WT diff gene analyse ##
##########################
Idents(cell) <- "CellType"

DimPlot(cell, reduction = "umap", label = TRUE)

# Define your comparison groups
comparison_groups <- c("Sympathoblasts", "CPC", "SCP", "ProlifSympathoblasts", "Bridging Cells")

# Initialize an empty list to store the results
all_results <- list()

# Loop through each comparison group and perform FindMarkers
for (group in comparison_groups) {
  result <- FindMarkers(cell, ident.1 = group)
  all_results[[group]] <- result
}
saveRDS(all_results, file = "DEGenes.RDS")

## ALK diff gene analyse ##
###########################
Idents(ALK) <- c("ALK","SC")
DimPlot(ALK, reduction = "umap", label = TRUE)
UMAPPlot(ALK, group.by = "Celltype")

comparison_groups <- c("ALK", "SC")
all_ALK_results <- list()

# Loop through each comparison group and perform FindMarkers
for (group in comparison_groups) {
  result <- FindMarkers(ALK, ident.1 = group)
  all_ALK_results[[group]] <- result
}
saveRDS(all_ALK_results, file = "DEGenes_ALK.RDS")
