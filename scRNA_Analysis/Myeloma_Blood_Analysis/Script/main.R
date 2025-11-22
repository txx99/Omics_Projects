# data was"normalized to an equal sequencing depth (approximately 32k read pairs per cell) using the cellranger aggr pipeline"
setwd("C:\\Users\\liv_u\\Desktop\\GitHub\\Omics_Projects\\scRNA_Analysis\\Myeloma_Blood_Analysis")

library(Seurat)
library(tidyverse)

# integrated transcriptome count data post depth equalization (depth equalization != gene count normalisation)
df = Read10X_h5("./Data/20k_WholeBlood_Fixation_Leukocyte_Isolation_2donor_aggr_count_filtered_feature_bc_matrix.h5")

# CreateSeuratObject
# (counts = count_df, min.cells= min_cell_count_for_feature_inclusion, min.features= min features count for cell inclusion, project)
sdf = CreateSeuratObject(counts = df, min.cells = 3, min.features = 200, project="Aggr Leukocytes")
rm(df)
# Preprocessing --------------
View(sdf@meta.data)
range(sdf$nCount_RNA)
range(sdf$nFeature_RNA)

# PercentageFeatureSet(df, pattern="^MT-" (for human mitochon genes))
sdf[['percent.MT']] <- PercentageFeatureSet(sdf, pattern="^MT-")
range(sdf$percent.MT) # max is <3%

VlnPlot(sdf, features=c("nCount_RNA", "nFeature_RNA", "percent.MT"), ncol=3)
ggsave("./Data/Plots/ViolinPlots_3_col.png")
FeatureScatter(sdf, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+geom_smooth(method = "lm")
ggsave("./Data/Plots/nCount_nFeature_Scatterplot.png")

# see some genes high transcript counts --> filter cells based on Violin tails
sdf <- subset(sdf, subset =  nCount_RNA < 20000 & nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.MT<10) #percent.MT filtering not useful here due to small range

sdf <- NormalizeData(sdf, normalization.method = 'LogNormalize', scale.factor = 10000)

# find most varied features 
sdf <- FindVariableFeatures(sdf, selection.method = "vst", nfeatures=2000)
top10 <- head(VariableFeatures(sdf), 10)
varFeatPlot<- VariableFeaturePlot(sdf)
LabelPoints(plot=varFeatPlot, points = top10, repel=TRUE)
ggsave("./Data/Plots/VariableFeaturePlot.png")

# scale data
allgenes <- rownames(sdf) 
head(allgenes)
sdf <- ScaleData(sdf, features=allgenes)
str(sdf)
# Formal class Assay5 layer
# NormalizeData layer
# FindVariableFeatures layer
# ScaleData layer

# linear dim redux - PCA
sdf <- RunPCA(sdf, features = VariableFeatures(sdf))
# dimensionality of data
ElbowPlot(sdf)

# Clustering -----------------
# 15 - 25 clusters for leukocytes

# > FindNeighbors()
# -> RunUMAP()
sdf <- FindNeighbors(sdf)
sdf<- RunUMAP(sdf, dims = 1:25)
DimPlot(sdf, reduction = "umap", label = TRUE)

# --> FindClusters(), determine best resolution 
sdf <- FindClusters(sdf, resolution = c(0.5, 0.7, 0.9, 1.1, 1.2, 1.3))

# DimPlot(sdf, group.by = "RNA_snn_res.0.5", label=TRUE)
# DimPlot(sdf, group.by = "RNA_snn_res.0.7", label=TRUE)
# DimPlot(sdf, group.by = "RNA_snn_res.0.9", label=TRUE)
DimPlot(sdf, group.by = "RNA_snn_res.1.1", label=TRUE)
# selecting .1.1 as the most reasonable clustering resolution
ggsave("./Data/Plots/UMAP_clustered.png")
# DimPlot(sdf, group.by = "RNA_snn_res.1.2", label=TRUE)
# DimPlot(sdf, group.by = "RNA_snn_res.1.3", label=TRUE)

# ID cells to their cluster
Idents(sdf) <- "RNA_snn_res.1.1" 
head(Idents(sdf))

saveRDS(sdf, file= "./Data/sdf_clustered.rds")




# # Cell type identification/cluster (FindAllMarkers()) -----------------
# # for faster Wil Rank Sum Test @FindAllMarkers()
# install.packages('devtools')
# devtools::install_github('immunogenomics/presto')
# # find the markers that are highly expressed in each cluster -> compare each cluster against all other clusters
#   # FindMarkers() -> compare one cluster against another cluster
#   # FindAllMarkers() -> compare each cluster against all other clusters
#   # FindConservedMarkers() -> find markers conserved across a condition even in different clusters  
#     # looking at clusters of cell type and w/i cluster variations of condition type; 
# 
# markers<-FindAllMarkers(sdf) #cluster-based
# 
# FeaturePlot()
# DotPlot()

# Manual Annotation: markers vs cell signatures ----

# T Cells
# CD4 T cells: CCR7, IL7R, LST1 low
# CD8 T cells: CD8A, CD8B
# Naive T: CCR7, LTB, IL7R
# Cytotoxic T: GZMB, NKG7, PRF1
# Regulatory T: FOXP3, IL2RA, CTLA4
# γδ T: TRDC, TRGC1/2
# 
# NK Cells: NKG7, GZMB, GNLY, KLRD1 (CD94), KLRF1, CD56, NCAM1
# 
# Monocytes
# Classical: LYZ, S100A8, S100A9, LST1, CTSS
# Non-classical: FCGR3A (CD16), MS4A7
# Inflammatory: IL1B, CCL3, CCL4
# 
# Dendritic Cells
# cDC1: CLEC9A, CADM1, XCR1
# cDC2: FCER1A, CD1C
# pDC: GZMB, IRF7, LILRA4, TCF4
# 
# B Cells
# Naive B: MS4A1, CD79A/B
# Memory B: CD27, TNFRSF13B
# Plasmablasts: MZB1, XBP1, PRDM1
# 
# MAIT cells: KLRB1, SLC4A10, IL7R
# 
# Neutrophils (if captured): S100A8, S100A9, CSF3R, LTF


# Subcluster Major Lineages ----
# for deeper profiling, collect all T cells / B cells and rerun PCA+UMAP --> find naive, memory, exhausted, proliferating subsets



# OPT ----
# Trajectory or pseudotime analysis 
# integrate Seurat with Monocle3 or RNA velocity tools to infer differentiation trajectories within specific cell lineages (like B cell maturation or T cell polarization).