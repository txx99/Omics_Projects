# scRNA Seq: Myeloma Blood Analysis
Using the Seurat package, scRNA sequencing data of blood samples from 1 healthy and 1 myeloma patient is preprocessed and visualised by UMAP clustering. 
The 19 clusters are then manually annotated by cell types based on differential expression analysis. \
The dataset of the 2 subjects is compiled from the source without distinguishing labels, so the final annotations & plots cannot distinguish 
between sample sources. 

## Workflow
### ./Script/
**main.R**: the main script for this analysis. 
### ./Data/
The dataset is not uploaded to this repository due to GitHub size constraints. Details on the dataset can be found under 'Credits'. \
**Plots/UMAP_namedClusters.png** &  **/Plots/UMAP_namedClusters_color2.png**: the final annotated UMAP outputs in different colour-schemes. \
**Plots/Marker_FeaturePlots**: folder containing various cell-marker gene sets used to verify cluster identities.

## Credits
### Data
Data Source: *Aggregate of Human Whole Leukocytes Isolated from Fixed Whole Blood of a Diseased Donor and a Healthy Donor (Next GEM)*, Flex dataset analyzed using Cell Ranger 8.0.0, 10x Genomics, (2024, March 27). \
Data Link: https://www.10xgenomics.com/datasets/20k_WholeBlood_Fixation_Leukocyte_Isolation_2donor_aggr

### Cluster Annotation Resources
CellMarker2.0 Database:\
Congxue Hu, Tengyue Li, Yingqi Xu, Xinxin Zhang, Feng Li, Jing Bai, Jing Chen, Wenqi Jiang, Kaiyue Yang, Qi Ou, Xia Li, Peng Wang, Yunpeng Zhang, CellMarker 2.0: an updated database of manually curated cell markers in human/mouse and web tools based on scRNA-seq data, Nucleic Acids Research, Volume 51, Issue D1, 6 January 2023, Pages D870–D876, https://doi.org/10.1093/nar/gkac947

PanglaoDB: \
Oscar Franzén, Li-Ming Gan, Johan L M Björkegren, PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data, Database, Volume 2019, 2019, baz046, doi:10.1093/database/baz046

Review paper on T cell annotation: \
Mullan KA, de Vrij N, Valkiers S, Meysman P. Current annotation strategies for T cell phenotyping of single-cell RNA-seq data. Front Immunol. 2023 Dec 21;14:1306169. doi: 10.3389/fimmu.2023.1306169. PMID: 38187377; PMCID: PMC10768068.

ChatGPT: \
OpenAI. (2025). ChatGPT (GPT 4.5) [Large language model]. https://chat.openai.com/chat
