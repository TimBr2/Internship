### SCTransform in Seurat ###
#############################
library(Seurat)
library(ggplot2)
library(sctransform)

### Load in the data ###
########################
setwd("C:/Users/timbr/OneDrive/Documenten/Howest/Internship/seurat_exercises")
pbmc_data <- Read10X(data.dir = "pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc_data)

### Normalization with SCTransform ###
######################################
# Single cmd to replace Normalizedata(), scaledata() and findvariablefeatrues()
# store mitochondrial percentage in object meta data
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
# to solve errors do first:
install.packages('BiocManager')
BiocManager::install('glmGamPoi')
# now this will work
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)


### Perform dimensionality reduction by PCA and UMAP embedding ###
##################################################################
# These are now standard steps in the Seurat workflow for visualization and clustering
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)

pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE)
DimPlot(pbmc, label = TRUE)

# individually annotate clusters based on canonical markers
# These are now standard steps in the Seurat workflow for visualization and clustering
# Visualize canonical marker genes as violin plots.
VlnPlot(pbmc, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7", "ISG15", "CD3D"),
        pt.size = 0.2, ncol = 4)

# Visualize canonical marker genes on the sctransform embedding.
FeaturePlot(pbmc, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7"), pt.size = 0.2,
            ncol = 3)

FeaturePlot(pbmc, features = c("CD3D", "ISG15", "TCL1A", "FCER2", "XCL1", "FCGR3A"), pt.size = 0.2,
            ncol = 3)
