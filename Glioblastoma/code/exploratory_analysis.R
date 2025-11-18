<<<<<<< HEAD
#Exploratory Analysis

library(GEOquery)
library(dplyr)
library(tidyverse)
library(ggplot2)

gset <- getGEO("GSE231994", GSEMatrix=TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL27956", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

gset 

#Pheno Data Exploratory Analysis
pheno_data <- pData(gset)
head(pheno_data)
str(pheno_data) #too messy maybe exclude
dim(pheno_data)

colnames(pheno_data)
table(pheno_data$characteristics_ch1.1)

#Expression Data Exploratory Analysis
expr_data <- exprs(gset) 
dim(expr_data)
hist(expr_data)
qqnorm(expr_data) #looks quite normal

sum(is.na(expr_data)) #no NAs

#Quality check of exprs data
boxplot(expr_data, main = "Boxplot of Expression Values", las=2, col="lightblue", outline=FALSE)
#Interpretation: want a similar distribution across our samples; medians are aligned
#Suggests data is comparable across our samples; no technical bias

pca <- prcomp(t(expr_data), scale. = TRUE)
pca_df <- data.frame(pca$x, Sample = rownames(pData(gset)))

pca_plot <- pca_df %>% 
            ggplot(aes(PC1, PC2, label = Sample)) + 
              geom_point(aes(color=pheno_data$characteristics_ch1.1)) + 
              geom_text(size=2, vjust = 1.5) + 
              theme_classic() + 
              labs(title = "PCA of Samples")
pca_plot









=======
#Exploratory Analysis

library(GEOquery)
library(dplyr)
library(tidyverse)
library(ggplot2)

gset <- getGEO("GSE231994", GSEMatrix=TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL27956", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

gset 

#Pheno Data Exploratory Analysis
pheno_data <- pData(gset)
head(pheno_data)
str(pheno_data) #too messy maybe exclude
dim(pheno_data)

colnames(pheno_data)
table(pheno_data$characteristics_ch1.1)

#Expression Data Exploratory Analysis
expr_data <- exprs(gset) 
dim(expr_data)
hist(expr_data)
qqnorm(expr_data) #looks quite normal

sum(is.na(expr_data)) #no NAs

#Quality check of exprs data
boxplot(expr_data, main = "Boxplot of Expression Values", las=2, col="lightblue", outline=FALSE)
#Interpretation: want a similar distribution across our samples; medians are aligned
#Suggests data is comparable across our samples; no technical bias

pca <- prcomp(t(expr_data), scale. = TRUE)
pca_df <- data.frame(pca$x, Sample = rownames(pData(gset)))

pca_plot <- pca_df %>% 
            ggplot(aes(PC1, PC2, label = Sample)) + 
              geom_point(aes(color=pheno_data$characteristics_ch1.1)) + 
              geom_text(size=2, vjust = 1.5) + 
              theme_classic() + 
              labs(title = "PCA of Samples")
pca_plot









>>>>>>> 5db9c95cae8c34e183c8b1a27d0cd955f1fd60d4
