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
ggsave("./Data/ViolinPlots_3_col.png")
FeatureScatter(sdf, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+geom_smooth(method = "lm")
ggsave("./Data/nCount_nFeature_Scatterplot.png")

# see some genes high transcript counts --> filter cells based on Violin tails
sdf <- subset(sdf, subset =  nCount_RNA < 20000 & nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.MT<10) #percent.MT filtering not useful here due to small range

sdf <- NormalizeData(sdf, normalization.method = 'LogNormalize', scale.factor = 10000)

# find most varied features 
sdf <- FindVariableFeatures(sdf, selection.method = "vst", nfeatures=2000)
top10 <- head(VariableFeatures(sdf), 10)
varFeatPlot<- VariableFeaturePlot(sdf)
LabelPoints(plot=varFeatPlot, points = top10, repel=TRUE)
ggsave("./Data/VariableFeaturePlot.png")

# scale data
allgenes <- rownames(sdf) 
head(allgenes)
sdf <- ScaleData(sdf, features=allgenes)
str(sdf)

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

saveRDS(sdf, file= "sdf_clustered.rds")


