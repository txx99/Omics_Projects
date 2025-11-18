# BINF5503 Final Project
# Authors: Usman Ahmed, Tazmeen Gill, Noha El-Haj
# Dataset: 1 - Cancer Transcriptomics

# ----- import required libraries -----
library(tidyverse)
# BiocManager::install("edgeR")
library(edgeR)
library(limma) #dependency of edgeR
library(tidyverse)
library(ggplot2)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)  
library(tidyverse)

# ----- import data -----
# paths to files
tsv_files <- list.files(path = "./Read_Count_Data/", pattern = "\\.tsv$", full.names = TRUE)
print(tsv_files)

# extract clean df names (patient ID + day)
df_names <- tsv_files %>% 
  str_replace_all(c("./.*/" = "", "_gene.tsv" = ""))

# list all the dfs (as tibbles) in one variable
all_dfs <- list()
for (file in unique(tsv_files)) {
  all_dfs[[file]] <- read_tsv(file = file, 
                              col_names = c("ensembl_gene_rec", 
                                            "gene_name", 
                                            "read_count"), 
                              na = "0")
}

# rename dfs from file path to patient + day 
all_dfs <- all_dfs %>% set_names(df_names)

# turn all listed dfs into individual objects
invisible(list2env(all_dfs, globalenv()))

# ----- create individual df for each patient -----
patient_ids <- unique(str_replace_all(string = df_names,
                                      pattern = "_D[06]$",
                                      replacement = ""))

for (id in patient_ids) {
  df_d0 <- get(paste0(id, "_D0")) # get data for id_D0 variable 
  df_d6 <- get(paste0(id, "_D6"))
  
  colnames(df_d0)[3] <- "d0_rc" # rename data col to include day
  colnames(df_d6)[3] <- "d6_rc"
  
  merged_df <- merge(df_d0, # merge d0 + d6 dfs 
                     df_d6, 
                     by = c("ensembl_gene_rec", "gene_name"))
  
  assign(paste0(id, "_df"), merged_df) #assign the id's merged df to id_df variable
}

# intermediary dfs (df_d0, df_d6, merged_df) currently hold final id's data (A899_df) == duplicate objects
rm(df_d0, df_d6, merged_df)

# create list of all patients for ease of function application
patients_list <- mget(ls(pattern = "^A.*_df$"))


# -------------- make inter-patient dfs by timepoint -----------
d0_dfs_list <- mget(ls(pattern = "D0$"))
d0_df<- purrr::reduce(.x=d0_dfs_list, merge, by=c("ensembl_gene_rec", "gene_name"), all=TRUE)
colnames(d0_df)[3:6]<- names(d0_dfs_list)

d6_dfs_list <- mget(ls(pattern = "D6$"))
d6_df<- purrr::reduce(.x=d6_dfs_list, merge, by=c("ensembl_gene_rec", "gene_name"), all=TRUE)
colnames(d6_df)[3:6]<- names(d6_dfs_list)

all_data<- merge(d0_df, d6_df, by= c("ensembl_gene_rec", "gene_name"), all=TRUE)


# ----- data exploration -----
lapply(patients_list, head)
lapply(patients_list, summary) 
lapply(patients_list, str)
lapply(patients_list[[1]], typeof)
lapply(patients_list[[1]], class)

# check for OCSC genes in datasets
raw_relevant_genes <- list()
for (i in 1:length(patients_list)) {
  raw_relevant_genes[[i]] <- patients_list[[i]] %>% filter(gene_name == "ALDH1A1"| 
                                                             gene_name =="CD44"| 
                                                             # gene_name =="CD133"| 
                                                             gene_name == "CD24"| 
                                                             # gene_name =="CD117"| 
                                                             gene_name == "EPCAM"| 
                                                             gene_name ==  "THY1")
  names(raw_relevant_genes)[i] <- names(patients_list)[i]
}


# ----- data cleaning & exporting -----
dir.create('./Clean_Data')
#intrapatient 
clean_patients <- list()
for (i in 1:length(patients_list)) {
  clean_df <- patients_list[[i]] %>% 
    filter(!grepl("^_", ensembl_gene_rec) & #removes first five rows relaying NA info
             !(is.na(d0_rc) & is.na(d6_rc)))
  
  clean_patients[[i]] <- clean_df 
  print(sum(duplicated(clean_df$gene_name)))
  write.csv(clean_df, paste0("./Clean_Data/clean_", names(patients_list)[i], '.csv'), row.names = FALSE)
}
rm(clean_df)

#interpatient
clean_all_data <- all_data[!apply(is.na(all_data[, 3:10]), 1, all), ] %>% 
  # looking at cols 3:10 by row, if all() cols TRUE for is.na(), exclude that row 
  filter(!grepl("^_", ensembl_gene_rec))
# remove the initial rows giving summary stats

sum(duplicated(clean_all_data$gene_name)) #TRUE

write.csv(clean_all_data, './Clean_Data/clean_all_data.csv', row.names = FALSE)

# ==============================

# edgeR DEA reference:
# https://github.com/rnnh/bioinfo-notebook/blob/master/docs/DE_analysis_edgeR_script.md

# --- load cleaned data ----

all_data_file<- "./Clean_Data/clean_all_data.csv"

# assign NAs as 0s
all_data<- read.csv(all_data_file)
all_data[is.na(all_data)]<- 0

#-----boxplots + scatterplots? ------
all_data_long <- all_data %>% 
  pivot_longer(!c(ensembl_gene_rec, gene_name), 
               names_to = c('patient', 'day'), 
               names_sep = '_',
               values_to = 'read_count')

all_data_wide <- all_data_long %>% 
  pivot_wider(names_from = day, 
              values_from = read_count)

all_data_wide %>%
  mutate(D0 = as.numeric(unlist(D0)),
         D6 = as.numeric(unlist(D6)),
         patient = as.factor(patient)) %>%
  ggplot(aes(x = D0, 
             y = D6)) +
  geom_point() +
  labs(x = 'D0 Read Count', 
       y = 'D6 Read Count') +
  facet_wrap(vars(patient)) +
  theme(legend.position = 'none', 
        panel.background = element_rect(fill = 'white', 
                                        colour = 'grey'))
ggsave("./Visualizations/exploratory_scatter.png", 
       width = 8, 
       height = 6)


# ----DGEList------
# create DGEList class object 
DGE_counts <- DGEList(counts = all_data,
                      genes = all_data$gene_name)
DGE_counts$genes$gene_name<-NULL

# add condition and patient grouping (factored) 
design_table<- data.frame(samples = colnames(all_data[2:9]))
design_table$patient<- as.factor(str_replace(design_table$samples, "_D[06]$",''))
design_table$day<- as.factor(c(rep("day 0",4), rep("day 6",4)))

# Confirm order matches 
summary(colnames(all_data[2:9]) == design_table$samples)

# Add grouping information to DGEList object
DGE_counts$samples$group<- as.factor(design_table$day)
DGE_counts$samples$group2<- as.factor(design_table$patient)


#------filter low expression genes---------
# filterByExpr keeps only genes meeting the min.count ( = 10) reads 
# will automatically look in your DGEList object for group info
genes_to_keep <- filterByExpr(DGE_counts)
summary(genes_to_keep)

# filter data 
DGE_counts <- DGE_counts[genes_to_keep, keep.lib.sizes = FALSE] #recalculate library sizes
dim(DGE_counts)[1]==summary(genes_to_keep)[3]

# ----------normalisation-----------------
# comparing gene expr bw samples, but different samples have different library sizes (total read count) --> can skew analysis
# normalise so that each sample has relatively similar impact on DEA analysis, reduce bias for high expr genes
DGE_counts$samples$norm.factors #currently all samples weighed equally

# calcNormFactors() normalizes across all samples based on lib sizes (assumes genes *not* differentially expressed)
# by scaling to minimise LogFC bw samples
DGE_counts <- calcNormFactors(DGE_counts)
DGE_counts$samples$norm.factors # see normalisation factors have been adjusted

#---------dispersion-------------------
# estimate gene dispersion == estimate relative variability of true expr bw replicates
# create design matrix (good practice)
condition_<- design_table$day
patient_<- design_table$patient

# ---- glm test (blocking intrapatient covariates) -------------

# define design = model matrix of what samples are being compared
design_2 <- model.matrix(~patient_+condition_)
DGE_counts<- estimateDisp(DGE_counts, design = design_2)

fit <- glmFit(DGE_counts, design_2)

# Likelihood ratio test
lrt <- glmLRT(fit, coef = "condition_day 6")
# looking at: change in gene expression at Day 6 after accounting for patient-specific differences

DEA_results <- topTags(lrt, n = Inf)$table
write.csv(DEA_results, 'DEA_results.csv', row.names = TRUE)


# ==================VISUALISATIONS ====================

# ---- find OCSC marker genes-------

all_genes<- topTags(lrt, n=Inf)
OCSC_markers<- list()
OCSC_markers <- all_genes$table %>% filter(genes == "ALDH1A1"| 
                                             genes =="CD44"| 
                                             genes == "CD24"| 
                                             genes == "EPCAM"| 
                                             genes ==  "THY1")
OCSC_markers # stats for all OCSC genes
rm(all_genes)


OCSC_markers %>% ggplot(aes(x = genes, 
                            y = logFC)) + 
  geom_col() + 
  geom_text(aes(y = 2, 
                label = round(FDR, digits = 4))) +
  labs(x = 'gene', 
       y = 'logFC Value') +
  theme_classic()
ggsave("./Visualizations/inter_expr.png", 
       width = 8, 
       height = 6)

OCSC_markers %>% EnhancedVolcano(lab = 'genes', 
                                 x = 'logFC' , 
                                 y = 'FDR',
                                 title = '',
                                 xlab = 'Log2 FC',
                                 ylab = 'FDR (p-adj)',
                                 pCutoff = 0.05, 
                                 FCcutoff = 2)
ggsave("./Visualizations/OCSC_volcano.png", 
       width = 8, 
       height = 6)

DEA_results %>% EnhancedVolcano(lab = DEA_results$genes, 
                                x = 'logFC' , 
                                y = 'FDR',
                                title = '',
                                xlab = 'Log2 FC',
                                ylab = 'FDR (p-adj)',
                                pCutoff = 0.05, 
                                FCcutoff = 2)

ggsave("./Visualizations/DEA_volcano.png", 
       width = 8, 
       height = 6)


# ---GSEA script ----------------

#Get DEA results, run first if not
dge_table <- DEA_results

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
