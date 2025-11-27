# Human peripheral blood mononuclear cells (PBMCs) of a healthy male.

setwd("path/to/folder")

library('Seurat')
library('tidyverse')

# load data
df <- Read10X_h5('./Data/sc5p_v2_hs_PBMC_10k_filtered_feature_bc_matrix.h5', use.names = TRUE)
df<-df$`Gene Expression`
# rename cols to remove '_totalC'
for (i in 1:length(rownames(df))){
  rownames(df)[i]<- sub('_TotalC','', rownames(df)[i])
}
rm(i)

# CreateSeuratObject(counts = count_df, min.cells= min_cell_count_for_feature_inclusion, min.features= min features count for cell inclusion, project= "Project name, shows as orig.ident col")
sdf<- CreateSeuratObject(df, min.cells = 5, min.features = 200, project = "healthy PBMC")
rm(df)
# Preprocess & perform cell clustering + annotation ------------------
# --> remove low quality cells

# percent mitochondrial genes
sdf[["percent.MT"]] <- PercentageFeatureSet(sdf, pattern="^MT-") 

# visualise QC metrics
VlnPlot(sdf, feature=c('nCount_RNA', 'nFeature_RNA', 'percent.MT'), ncol = 3)
ggsave("./Data/Plots/ViolinPlots.png")
FeatureScatter(sdf, feature1 = "nCount_RNA", feature2="nFeature_RNA")+geom_smooth(method="lm")
ggsave("./Data/Plots/ScatterPlot.png")


# Filtering Cells (Additional QC Metrics

# subset(df, subset= nCount_RNA <> condition1 & nFeature_RNA <> condition 2 & percent.MT <> condition3)
# cutoffs based on Violin tails + normal range for this cell type
range(sdf$nCount_RNA)
range(sdf$nFeature_RNA)
range(sdf$percent.MT)
sdf <- subset(sdf, subset=nFeature_RNA>200 & nFeature_RNA<3200 & nCount_RNA < 20000 & percent.MT<10)


# Normalisation for Comparability
sdf<- NormalizeData(sdf, normalization.method = "LogNormalize", scale.factor = 10000)



## Filtering Highly Variable Features
# Highly variable genes exhibit cell to cell variation, highlighting the biological signal in sc datasets.

# ID highly variable genes 
sdf<- FindVariableFeatures(sdf, selection.method = "vst", nfeatures=2000)
# plot top 10 variable genes
top10 <- head(VariableFeatures(sdf), 10)
plot1<- VariableFeaturePlot(sdf)
LabelPoints(plot=plot1, points = top10, repel=TRUE)
ggsave("./Data/Plots/VariableFeaturesPlot.png", width = 8, height = 6)
rm(plot1)

## Scaling

allgenes <- rownames(sdf)
head(allgenes)
sdf <- ScaleData(sdf, features=allgenes)
str(sdf)
rm(allgenes)


## Linear Dimensionality Reduction - PCA ----------------------------
sdf <- RunPCA(sdf, features = VariableFeatures(sdf))

print(sdf[["pca"]], dim=1:5, nfeatures=5)
# heatmap:
#   X== cells, y== PC components, yellow== high expression, purple == low expression 
DimHeatmap(sdf, dims=1:5, cells=500, balanced=TRUE)
ggsave("./Data/Plots/PC_VariableFeatures_Heatmap.png", width = 8, height = 6)

# determine dimensionality of data 
ElbowPlot(sdf)


## Clustering -------------------------------------------------------

sdf <- FindNeighbors(sdf, dims=1:15)

## TSNE//UMAP -- BEFORE resolution, group cells of similar types together at low dim
sdf <- RunUMAP(sdf, dims=1:25)
DimPlot(sdf, reduction = "umap", label=TRUE)

sdf <- FindClusters(sdf, resolution = c(0.1, 0.3, 0.5, 0.7, 0.9, 1))

# visualise different resolutions + select the best as Idents 
DimPlot(sdf, group.by = "RNA_snn_res.0.1", label=TRUE)
DimPlot(sdf, group.by = "RNA_snn_res.0.3", label=TRUE)
DimPlot(sdf, group.by = "RNA_snn_res.0.5", label=TRUE)
DimPlot(sdf, group.by = "RNA_snn_res.0.7", label=TRUE)



DimPlot(sdf, group.by = "RNA_snn_res.0.3", label=TRUE)
ggsave("./Data/Plots/UMAP_Clustered.png", width = 8, height = 6)

# ID cells to their clusters
Idents(sdf) <- "RNA_snn_res.0.3" # set identity w resolution 
Idents(sdf)






# Cluster annotation -------------------------------------------------

# FindMarkers() -> compare one cluster against another cluster
# FindAllMarkers() -> compare each cluster against all other clusters
# FindConservedMarkers() -> find markers conserved across a condition even in different clusters
# looking at clusters of cell type and w/i cluster variations of condition type;


markers <- FindAllMarkers(sdf, # ensure default assay is RNA w DefaultAssay(data)
                          logfc.threshold = 0.25, 
                          min.pct = 0.1, 
                          only.pos= FALSE, 
)


# create list of top upregulated DE genes per cluster
for (i in 0:(length(levels(Idents(sdf)))-1)){
  genes <- head(markers[markers$cluster==i, "gene"], 4)
  assign(paste0("cluster_", i), genes)
  print(FeaturePlot(sdf, features = genes, min.cutoff = "q10"))
  
}


# FeaturePlot of top upregulated genes of each cluster
print(FeaturePlot(sdf, features = cluster_0, min.cutoff = "q10"))

# FeatuerPlot of specific cell type marker genes

geneset<- c('TUBB1', 'PPBP', 'PF4', 'CAVIN2')
print(FeaturePlot(sdf, features = geneset, min.cutoff = "q10"))
ggsave('./Data/Plots/marker_featurePlots/megakaryocyte_markers.png')



# rename clusters
sdf<- RenameIdents(sdf, 
                   `0`='antigen-presenting classical monocytes',
                   `1`= 'naive/central memory T' ,
                   `2`= 'memory CD4+ T',
                   `3`= 'naive CD8+ T'  ,
                   `4`= 'B cells',
                   `5`= 'inflammatory classical monocytes',
                   `6`= 'NK', # + gamma delta T mixed in here somewhat
                   `7` = 'effector memory T' ,
                   `8`= 'cytotoxic T/ MAIT ',
                   `9`= 'non-classical monocyte', 
                   `10`= 'megakaryocyte',
                   `11`= 'pDC')


DimPlot(sdf, reduction = "umap", label = TRUE, label.size = 3.5, pt.size = 0.7)

ggsave('./Data/Plots/UMAP_namedClusters.png', height=8, width=10)

# save data
saveRDS(sdf, file= "./Data/seuratObject_namedClusters.rds")

