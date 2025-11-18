setwd("C:\\Users\\liv_u\\Desktop\\GitHub\\Omics_Projects\\scRNA_analysis")
# scRNA Data Analysis in R
# Following [Arpudhamary V.'s](https://github.com/Marydoss-25/scRNA-Data-analysis-) tutorial for beginners in scRNA analysis, this notebook will learn to use the Seurat package and will perform:,
    # Seurat object conversion,
    # QC,
    # Filtering,
    # Normalisation,
    # ID highly variable genes,
    # Data scaling,
    # Dimensionality reduxn,
    # Clustering,
    # Non-linear dimensionality reduction,
    # Cluster visualisation"

# install.packages(c("Seurat", "reticulate"))
# reticulate::py_install("umap-learn")

# load libraries 
library(Seurat) # for sc analysis
library(tidyverse)

# DATASET ---- 
## Credits ----
# Source: 20k Mixture of NSCLC DTCs from 7 donors, 3' v3.1 (without intronic reads), Universal 3' dataset analyzed using Cell Ranger 6.1.2, 10x Genomics, (2021, August 09).
# Link: https://www.10xgenomics.com/datasets/20-k-mixture-of-nsclc-dt-cs-from-7-donors-3-v-3-1-3-1-standard-6-1-0

# Gene Expression - Feature / cell matrix HDF5 (raw)
# Pooled multiplexed sample:   
#   Estimated number of cells: 16,443
#   Cells assigned to a sample: 12,231
# This dataset is licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0) license.

## Load data ----
  # Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  # -> Read count matrix from 10X CellRanger hdf5 file. This can be used to read both scATAC-seq and scRNA-seq matrices.
nsclc_raw <- Read10X_h5(filename="./Data/20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")
# select Gene Expression table
nsclc_counts<- nsclc_raw$`Gene Expression`
# view structure 
str(nsclc_counts) 
# rm(nsclc_raw)

# create seurat obj
seurat_nsclc <- CreateSeuratObject(counts = nsclc_counts, project="NSCLC", min.cells = 3, min.features = 200)

# Code ----
## QC -------------------------
# Current sc extraction methods take min 6h for tissue dissociation, sc dopamine capture, and mRNA extraction.
# Extraction procedures can result in stressed or dead cells == low quality cells 
# --> remove low quality cells during QC

# QC metrics - nFeature_RNA / nCount_RNA - [A UMI (Unique Modular Identifier) - a short, random DNA sequence added to each RNA molecule during scRNA sequencing] / percent.mt
# rows = cell UMI
# cols = transcripts (nCount), genes (nFeature)
# 1. Low quality cells/empty droplets (no cells) often have very few genes == low Feature RNA & low Count RNA
# 2. Cell doublets or multiples have high values of Feature RNA & Count RNA
# 3. low quality/dying cells often have high % of mitochondrial genes (percent.mt)
# 

# view QC metrics
View(seurat_nsclc@meta.data) # each row represents 1 cell --> shows the number of transcripts/UMIs and genes in one cell
range(seurat_nsclc$nCount_RNA) # range of transcript values
range(seurat_nsclc$nFeature_RNA) # range of genes expression values

# percent mitochondrial genes
seurat_nsclc[["percent.MT"]] <- PercentageFeatureSet(seurat_nsclc, pattern="^MT-")
View(seurat_nsclc@meta.data)

# visualise QC metrics
# individual feature violin plots, ncol=3
VlnPlot(seurat_nsclc, features=c("nCount_RNA", "nFeature_RNA", "percent.MT"), ncol=3)
ggsave("./Data/ViolinPlots.png", width = 8, height = 6)
# correlation bw features w lin reg trendline 
FeatureScatter(seurat_nsclc, feature1 = "nCount_RNA", feature2="nFeature_RNA")+geom_smooth(method="lm")
ggsave("./Data/ScatterPlots.png", width = 8, height = 6)
# based on the Feature Scatter plot, we can further filter artefacts
# A) if cells accumulate in lower right corner of the plot [Count=150k, Fetaure=10k], indicates exp captured few genes which are sequenced repeatedly, leading to high transcript counts
# B) if cells accumulate in the top left corner [Count=20k, Feature=30k], indicates exp captured a high number of genes but they were not deeply sequenced 
# ideally, would want transcript-gene ratio to follow the linear trendline.

# In our data, some cells have few genes with high transcripts (scenario A), which we want to filter out



## Filtering Cells (Additional QC Metrics) ------------------
# 1. Higher presence of ribosomal genes = lower quality cell
# 2. doublet = 2+ cells captured in the same droplet during capture --> mixed gene expression profile --> want to filter out doublets in sample

# filter for cells with: less than 5% MT genes & 200 to 2500 genes (<200 is low quality cell, >2500 is likely low depth?) 
    # NOTE: the threshold you sets depends on your data and tissue; heart tissue has higher MT content while lymphocytes have lower MT content.

seurat_nsclc <- subset(seurat_nsclc, subset=nFeature_RNA>200 & nFeature_RNA<2500 & percent.MT<5)


## Normalisation for Comparability -------------------
# for expression values to be comparable across cells, must normalise to account for differences in sequencing depth & library size.
# Some cells have more total reads per UMI than others; if not normalised, highly sequenced cells would dominate the results.

# Default normalisation method is log normalisation of (gene expr in each cell / total gene expression) * scaling factor
seurat_nsclc<- NormalizeData(seurat_nsclc, normalization.method = "LogNormalize", scale.factor = 10000)



## Filtering Highly Variable Features -----------------
# Highly variable genes exhibit cell to cell variation, highlighting the biological signal in sc datasets.

# ID highly variable genes 
seurat_nsclc<- FindVariableFeatures(seurat_nsclc, selection.method = "vst", nfeatures=2000)
# plot top 10 variable genes
top10 <- head(VariableFeatures(seurat_nsclc), 10)
plot1<- VariableFeaturePlot(seurat_nsclc)
LabelPoints(plot=plot1, points = top10, repel=TRUE)
ggsave("./Data/VariableFeaturesPlot.png", width = 8, height = 6)
  # variable gene count = 2000
  # non-variable gene count = 27553

## Scaling ---------------------------------
# To remove unwanted sources of variation due to biology, like being in different stages of the cell cycle or technical noise from batch effects, we scale.
# scaling --> highly expressed genes dont dominate ds analysis (PCA, clustering, etc)

# scale across all genes
allgenes <- rownames(seurat_nsclc)
head(allgenes)
seurat_nsclc <- ScaleData(seurat_nsclc, features=allgenes)
str(seurat_nsclc)

# under RNA ASSAY, 3 slots:
# Count = raw sparse matrix
# data = log normalised counts
# scale.data = scaled data


## Linear Dimensionality Reduction - PCA --------------
seurat_nsclc <- RunPCA(seurat_nsclc, features = VariableFeatures(seurat_nsclc))
# visualise
print(seurat_nsclc[["pca"]], dim=1:5, nfeatures=5)
DimHeatmap(seurat_nsclc, dims=1:5, cells=500, balanced=TRUE)
# heatmap:
#   X== cells, y== PC components, yellow== high expression, purple == low expression 
ggsave("./Data/PC_Heatmap.png", width = 8, height = 6)

# determine dimensionality of data 
# use elbow plot to ID significant PCs which capture majority of biological signals 
ElbowPlot(seurat_nsclc)

## Clustering ----------------------------
# Find neighbours 
# -> Find clusters; resolution of clusters: lower number == fewer clusters; resolution ranges 0 to 1 or more to see the best clusters
# -> Choose the best resolution using DimPlot

seurat_nsclc <- FindNeighbors(seurat_nsclc, dims=1:15)
seurat_nsclc <- FindClusters(seurat_nsclc, resolution = c(0.1, 0.3, 0.5, 0.7, 0.9, 1))
# visualise the differnet resolutions 
DimPlot(seurat_nsclc, group.by = "RNA_snn_res.0.1", label-TRUE)
DimPlot(seurat_nsclc, group.by = "RNA_snn_res.0.3", label=TRUE)
DimPlot(seurat_nsclc, group.by = "RNA_snn_res.0.5", label=TRUE)
# Here, we see that resolution of 0.1 has good clustering and the subsequent increases in res do not change the clusters much, so we can stick to 0.1

# set ID of clusters
Idents(seurat_nsclc) # levels indicate number of clusters 
Idents(seurat_nsclc) <- "RNA_snn_res.0.1" # set identity w resolution 
Idents(seurat_nsclc)

## Non Linear Clustering - TSNE//UMAP --------------
# After clustering, group cells of similar types together at low dim
seurat_nsclc <- RunUMAP(seurat_nsclc, dims=1:15)
DimPlot(seurat_nsclc, reduction = "umap", label=TRUE)
ggsave("./Data/UMAP.png", width = 8, height = 6)