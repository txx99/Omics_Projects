# BINF5503 Final Project
# Authors: Usman Ahmed, Tazmeen Gill, Noha El-Haj
# Dataset: 1 - Cancer Transcriptomics

# ----- import required libraries -----
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

