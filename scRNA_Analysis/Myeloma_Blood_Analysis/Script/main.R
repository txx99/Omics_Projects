# data was"normalized to an equal sequencing depth (approximately 32k read pairs per cell) using the cellranger aggr pipeline"
# integrated transcriptome count data post depth equalization (depth equalization != gene count normalisation)

setwd("C:\\Users\\liv_u\\Desktop\\GitHub\\Omics_Projects\\scRNA_Analysis\\Myeloma_Blood_Analysis")

library(Seurat)
library(tidyverse)
# for faster Wilcox Rank Sum Test @FindAllMarkers():
# install.packages('devtools')
devtools::install_github('immunogenomics/presto')

# load data
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
# # DimPlot(sdf, group.by = "RNA_snn_res.0.7", label=TRUE)
# # DimPlot(sdf, group.by = "RNA_snn_res.0.9", label=TRUE)
DimPlot(sdf, group.by = "RNA_snn_res.0.7", label=TRUE)
# selecting .0.7 as the most reasonable clustering resolution
ggsave("./Data/Plots/UMAP_clustered_07.png")
# # DimPlot(sdf, group.by = "RNA_snn_res.1.2", label=TRUE)
# # DimPlot(sdf, group.by = "RNA_snn_res.1.3", label=TRUE)

# ID cells to their cluster
Idents(sdf) <- "RNA_snn_res.0.7" 
head(Idents(sdf))

saveRDS(sdf, file= "./Data/seuratObject_clustered.rds")




# Manual cluster annotation -----------------

# find the markers that are highly expressed in each cluster -> compare each cluster against all other clusters
  # FindMarkers() -> compare one cluster against another cluster
  # FindAllMarkers() -> compare each cluster against all other clusters
  # FindConservedMarkers() -> find markers conserved across a condition even in different clusters
    # looking at clusters of cell type and w/i cluster variations of condition type;


# FindAllMarkers(data, # ensure default assay is RNA w DefaultAssay(data)
#                logfc.threshold = 0.25, # 0.25 standard; min log2fc avg expression for considering a gene differentially expressed in one cluster over all the other clusters 
#                min.pct = 0.1, # 0.1 standard; min % of cells in the cluster that must present the gene 
#                only.pos= FALSE, # only want markers that are upregulated (not downregulated)
#                test.use = 'DESeq2', # if you want to specify package/test to use
#                slot = 'counts' # auto adjusts if not specified; appropriate input layer for the test being used (DESeq2 uses raw counts)
#                )

# Output Interpretation:
# avg_log2FC = log2FC seen in current cluster
# pct.1 = in current cluster, % of cells presenting target gene 
# pct.2 = in all other clusters, avg % of cells presenting target gene

# ideally would be using FindConservedMarkers, but both patients' data is aggregated without condition labels
# when annotating, may need to adjust resolution or FindMarkers thresholds
markers<-FindAllMarkers(sdf,
                        logfc.threshold = 0.25, # min log2fc avg expression for considering a gene differentially expressed in one cluster over all the other clusters 
                        min.pct = 0.1,
                        only.pos = TRUE) # min % of cells in cluster that must present the gene

for (i in 0:18){
  genes <- head(markers[markers$cluster==i & markers$avg_log2FC>0, "gene"], 4)
  assign(paste0("group_", i), genes)
}
# FeaturePlots of top genes of each cluster  
print(FeaturePlot(sdf, features = group_18, min.cutoff = "q10"))



# Manual Annotation: markers vs cell signatures --------

# FeatuerPlots of specific cell type's markers
geneset<- c('GZMB', 'TCF4', 'IRF7', 'CLEC4C', 'IL3RA', 'LILRA4')
print(FeaturePlot(sdf, features = geneset, min.cutoff = "q10"))
ggsave('./Data/Plots/Marker_FeaturePlots/pDC_markers.png')

sdf<- RenameIdents(sdf,)

sdf<- RenameIdents(sdf, `0`= 'Neutrophils',
                   `1`='Interferon Neutrophils', 
                   `2`='Inflammatory Neutrophils',
                   `4`= 'High-stress Myeloids',
  
                   `13`= 'Megakaryocyte',
                   `17`= 'Immature neutrophils',
                   
                   `11`= 'Naive B cell',
                   `16` = 'Basophil',
                   `9`= 'Tumor-activated classical monocyte',
                   `8`= 'DC-like classical monocyte',
                   `15`= 'TAM-like non-classical monocytes',
                   `18`= 'Plasmacytoid DC',
                   `3`= 'CD8+ T (naive)',
                   `7`= 'CD4+ T (naive)',
                   `5`= 'Regulatory T cells',
                   `14` = 'MAIT',
                   `12`= 'γδ T cells',
                   `6`= "Cytotoxic T",
                   `10`= 'NK')

sdf@active.ident
DimPlot(sdf, reduction = "umap", label = TRUE, label.size = 3, cols = 'alphabet')
ggsave('./Data/Plots/UMAP_namedClusters_color2.png', height=8, width=10)

saveRDS(sdf, file= "./Data/seuratObject_namedClusters.rds")



# Subcluster Major Lineages ----
# for deeper profiling, collect all T cells / B cells and rerun PCA+UMAP --> find naive, memory, exhausted, proliferating subsets



# OPT ----
# Trajectory or pseudotime analysis 
# integrate Seurat with Monocle3 or RNA velocity tools to infer differentiation trajectories within specific cell lineages (like B cell maturation or T cell polarization).