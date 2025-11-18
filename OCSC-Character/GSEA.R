#GSEA script

#libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)  
library(tidyverse)

#Get DEA results, run first if not
dge_table <- read.csv('DEA_results.csv')

#remove NA gene name rows
dge_table <- dge_table %>% filter(!is.na(genes))

#Remove duplicates
dge_table <- dge_table %>%
  group_by(genes) %>%
  slice_min(FDR, with_ties = FALSE) %>%
  ungroup()

#Map symbols to EntrezID
gene_ids <- bitr(dge_table$genes, fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = org.Hs.eg.db)

#Merge IDs to DEA results table
dge_mapped <- merge(dge_table, gene_ids, by.x = "genes", by.y = "SYMBOL")

#Vector for GSEA
gene_list <- dge_mapped$logFC
names(gene_list) <- dge_mapped$ENTREZID

#Sort by decreasing
gene_list <- sort(gene_list, decreasing = TRUE)

#Run GSEA with GO
gsea_go <- gseGO(geneList = gene_list,
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 ont = "BP",
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,  
                 verbose = FALSE)

#Visualize GSEA results
if (nrow(gsea_go@result) == 0) {
  print("No significant GO terms found under the specified p-value cutoff.")
} else {
  print("Significant GO terms found. Displaying plots...")
  #show top GO results
  print(head(gsea_go@result, 5))
  
  #Save Dotplot
  ggsave("./Visualizations/dotplot.png",
         plot = dotplot(gsea_go, showCategory = 10),
         width = 8, height = 6)
  
  #Save GSEA plot
  ggsave("./Visualizations/gsea_plot.png",
         plot = gseaplot2(gsea_go, geneSetID = gsea_go@result$ID[1]),
         width = 8, height = 6)

  #Save Ridgeplot
  ridge_plot <- ridgeplot(gsea_go, showCategory = 15) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8)) +
    ggplot2::labs(x = "Gene rank (Day 6 → Day 0)")
  
  ggsave("./Visualizations/ridgeplot.png",
         plot = ridge_plot,
         width = 8, height = 6)}
