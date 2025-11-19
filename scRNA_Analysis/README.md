# scRNA Data Analysis in R
Following [Arpudhamary V.'s](https://github.com/Marydoss-25/scRNA-Data-analysis-) tutorial for beginners in scRNA analysis, this notebook will learn to use the Seurat package for single-cell analysis and will perform: \
    - Seurat object conversion, \
    - QC, \
    - Filtering, \
    - Normalisation, \
    - ID highly variable genes, \
    - Data scaling, \
    - Dimensionality reduction, \
    - Clustering (DimPlot), \
    - Non-linear dimensionality reduction (PCA), and \
    - Cluster visualisation (UMAP) 
## Credits
Code: This notebook follows the Seurat tutorial by [Arpudhamary V.'s](https://github.com/Marydoss-25/scRNA-Data-analysis-). \
Data Source: *20k Mixture of NSCLC DTCs from 7 donors, 3' v3.1 (without intronic reads)*, Universal 3' dataset analyzed using Cell Ranger 6.1.2, 10x Genomics, (2021, August 09). \
Data Link: https://www.10xgenomics.com/datasets/20-k-mixture-of-nsclc-dt-cs-from-7-donors-3-v-3-1-3-1-standard-6-1-0 \
Specific dataset: Gene Expression - Feature / cell matrix HDF5 (raw) \
Pooled multiplexed sample:   
  Estimated number of cells: 16,443 \
  Cells assigned to a sample: 12,231

This dataset is licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0) license.
