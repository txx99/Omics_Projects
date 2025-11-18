# edgeR DEA reference:
# https://github.com/rnnh/bioinfo-notebook/blob/master/docs/DE_analysis_edgeR_script.md

# -----library -------
# BiocManager::install("edgeR")
library(edgeR)
library(limma) #dependency of edgeR
library(tidyverse)
library(ggplot2)
library(EnhancedVolcano)

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




