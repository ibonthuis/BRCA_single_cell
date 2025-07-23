source("/storage/kuijjerarea/ine/breast_met/single_cell/pre_processing_for_networks/run_scorpion_fn.R")

required_libraries <- c(
    "data.table",    
    "dplyr",
    "stringr",
    "Matrix",
    "SCORPION", # V1.02 because I don't want the 0 rows to be filtered.....
    "Seurat",
    "parallel",
    "furrr"
    )
    
for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
}

main_dir <- "/storage/kuijjerarea/ine/breast_met/single_cell"
setwd(main_dir)

TF <- fread("Input_network/motif_prior_names_2024.tsv")
PPI <- fread('Input_network/ppi_prior_2024.tsv')

ALL <- get(load("Data/LN_METS_merged_ANNOTATED.RData"))

count_matrix <- LayerData(ALL)
#colnames(count_matrix) <- ALL$patient_sampleID
count_matrix <- count_matrix[rowSums(count_matrix != 0) > 5,] 

TF <- TF %>% 
    filter(V1 %in% rownames(count_matrix)) %>%
    filter(V2 %in% rownames(count_matrix))

PPI <- PPI %>% 
    filter(V1 %in% rownames(count_matrix)) %>%
    filter(V2 %in% rownames(count_matrix))


count_matrix <- count_matrix[rownames(count_matrix) %in% TF$V2 , ]
list_of_subsamples <- purrr::map(IDs_of_interest, ~ subsample(.x, ALL, count_matrix))

list_of_subsamples_unlist <- unlist(list_of_subsamples)
names(list_of_subsamples) <- gsub("/", "", names(list_of_subsamples) )
save(list_of_subsamples, file = 'Results/Subsamples/all_subsamples.RData')
