# continue to Seurat from Scanpy (based on output from P02_Scanpy_processing.ipynb)
library(Seurat)
library(dplyr)
library(SeuratDisk)

# NOTE! predefine the "processed_data" path
USER_PATH = ".../SpatioTemporal_malaria_liver_stage_atlas/processed_data"
  
######################### 1. Load data ##########################

# load Data

Convert(file.path(USER_PATH, "CS_data_processed_prescaling.h5ad"), dest = "h5seurat", overwrite = TRUE)
CS_data <- LoadH5Seurat(file.path(USER_PATH, "CS_data_processed_prescaling.h5seurat"))

PBA_genes=read.csv(file.path(gsub("processed_data","input_data",USER_PATH),"protists_ensembl_48_PBANKA_geneID-name.csv"))
PBA_genes$external_gene_name = gsub(PBA_genes$external_gene_name,pattern = "_",replacement = "-")

CS_data[["PBA"]]=CreateAssayObject(CS_data@assays$RNA@counts[rownames(CS_data) %in% PBA_genes$external_gene_name,])
CS_data[["MUS"]]=CreateAssayObject(CS_data@assays$RNA@counts[!(rownames(CS_data) %in% PBA_genes$external_gene_name),])

CS_data

######################### 2. Normalization & UMAP  ##########################
  
#normalized data by LogNormalize method
CS_data <- NormalizeData(CS_data)

# feature selection: keep the top 2000 highly variable features
CS_data <- FindVariableFeatures(CS_data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(CS_data), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(CS_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)

# Data scaling and regression
all.genes<-rownames(CS_data)
CS_data <- ScaleData(CS_data, features = all.genes)

# PCA
CS_data <- RunPCA(CS_data, npcs = 50, verbose = FALSE)

ElbowPlot(CS_data,ndims = 30)

# select no. of PCs for UMAP reduction
dim = 10 

# UMAP and Clustering
CS_data <- FindNeighbors(CS_data, reduction = "pca", dims = 1:dim)
CS_data <- FindClusters(CS_data, resolution = 0.4)

CS_data <- RunUMAP(CS_data, reduction = "pca", dims = 1:dim)

# flip UMAP coordinates
CS_data[["umap"]]@cell.embeddings[,2]=-CS_data[["umap"]]@cell.embeddings[,2]

# plot UMAP - Time
DimPlot(CS_data, reduction = "umap", group.by = c("coarse_time"),pt.size = 2)

######################### 3. INFECTED Cells subset ##########################
PBA_data = subset(CS_data,subset = infected=="TRUE")

DefaultAssay(object = PBA_data) <- "PBA"

#normalized data by LogNormalize method
PBA_data <- NormalizeData(PBA_data)

# feature selection: keep the top 2000 highly variable features
PBA_data <- FindVariableFeatures(PBA_data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(PBA_data), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(PBA_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)


# Data scaling and regression
all.genes<-rownames(PBA_data)
PBA_data <- ScaleData(PBA_data, features = all.genes)

# PCA
PBA_data <- RunPCA(PBA_data, npcs = 50, verbose = FALSE)


ElbowPlot(PBA_data,ndims = 30)
dim = 10 

# UMAP and Clustering
PBA_data <- FindNeighbors(PBA_data, reduction = "pca", dims = 1:dim)
PBA_data <- FindClusters(PBA_data, resolution = 0.4)

PBA_data <- RunUMAP(PBA_data, reduction = "pca", dims = 1:dim, metric="euclidean")

# flip UMAP coordinates
PBA_data[["umap"]]@cell.embeddings[,1]=-PBA_data[["umap"]]@cell.embeddings[,1]

# plot UMAP - Time
DimPlot(PBA_data, reduction = "umap", group.by = c("coarse_time"),pt.size = 2)