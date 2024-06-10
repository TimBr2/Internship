### Seurat tutorial ###
#######################
# https://satijalab.org/seurat/articles/pbmc3k_tutorial
setwd("C:/Users/timbr/OneDrive/Documenten/Howest/Internship/seurat_exercises")
library(dplyr)
library(Seurat)
library(patchwork)

### Setup Seurat object ###
###########################
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# How does the data look?
# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

dense.size <- object.size(as.matrix(pbmc.data))
dense.size

# sparse because most values in scRNA-seq matrix are 0
sparse.size <- object.size(pbmc.data)
sparse.size

### Pre-processing ###
######################
## QC and selecting cells for further analysis##
################################################
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# we calculate mitochondrial QC metrics with the PercentageFeatureSet function and MT- as a set of mitochondrial genes, low-quality/dying cells have mostly extensive mitochondrial contamination
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

# Where are QC metrics stored?
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# we filter unique feature counts over 2500 and less than 200, also filter cells >5% mitochondrial counts
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc


### Normalize the data ###
##########################
# globalscaling normalization with LogNormalize and default scaling factor 10 000
# This normalizes the feature expansion measurement for each cell by the total expression, multiplies by a scale factor (10,000 by default), and log-transforms the result.
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
# this is the same as above with the default method, could also use SCTransform
pbmc <- NormalizeData(pbmc)

###  Identification of highly variable features ###
###################################################
# Identify the 2000 most highly variable genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


### Scaling the data ###
########################
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# How to remove unwanted sources of variation
# remove the percent.mt column (unwanted sources of variation)
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


### Perform linear dimensional reduction ###
############################################
# output list of genes with most positive and negative loadings
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca") + NoLegend()
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


### Determine dimensionality of the dataset ###
###############################################
# to rank principle components based on percentage of variance, look at the elbow, here 9-10, so true signal captured at first 10PCs
ElbowPlot(pbmc)


#### Cluster the cells ###
##########################
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


### Run non-linear dimensional reduction (UMAP/tSNE) ###
########################################################
# to visualise and explore datasets
pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
#saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")


### Find differntially expressed features (cluster biomarkers)###
#################################################################
# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# find markers for every cluster compared to all remaining cells, report only the positive
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster0.markers, n = 5)
# visualizing marker expression, shows expression probability distributions across clusters
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))
# generate expression heatmap for given cells and featurs, top 20 markers for each cluster
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


### Assigning cell type identity to clusters ###
################################################
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
