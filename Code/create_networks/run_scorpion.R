# Be careful, it took 10 days to run 1555 networks

required_libraries <- c(
    "data.table",    
    "dplyr",
 #   "optparse",
  #  "rlang",
   # "ggplot2",
    #"purrr",
    #"tibble",
    "stringr",
    #"limma",
    #"tidyr",
    "Matrix",
    "SCORPION", # V1.02 because I don't want the 0 rows to be filtered.....
    "Seurat",
    "parallel",
    "furrr",
    "future"
    )

for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
}



main_dir <- "/storage/kuijjerarea/ine/breast_met/single_cell"
setwd(main_dir)

# File structure: one main folder containing folders "Data", "Input_network" and "Results". "Results" should contain a subfolder "Subsamples"

# Make sure TF and PPI are unfiltered
TF <- fread("Input_network/motif_prior_names_2024.tsv")
PPI <- fread('Input_network/ppi_prior_2024.tsv')
list_of_subsamples <- get(load("Results/Subsamples/all_subsamples.RData"))
print("Read in TF, PPI and gex done")

TF <- TF %>% 
    filter(V1 %in% rownames(list_of_subsamples[[1]])) %>%
    filter(V2 %in% rownames(list_of_subsamples[[1]]))

PPI <- PPI %>% 
    filter(V1 %in% rownames(list_of_subsamples[[1]])) %>%
    filter(V2 %in% rownames(list_of_subsamples[[1]]))



#list_of_subsamples <- unlist(list_of_subsamples)

options(future.globals.maxSize= 4891289600)

# length(list_of_subsamples)
# names(list_of_subsamples[1])
# dim(list_of_subsamples[[1]])

plan(multisession, workers = 4)
print("Plannned the multisession")

#1555 networks --> 20 cores
furrr::future_map2(list_of_subsamples, names(list_of_subsamples), function(df, name){
    #print(class(df))
    print(name)
    #print(paste0("Started future_map for sample ", names(df)))
    outFile <- paste0('Results/Networks/', name, ".RData")
    if(!file.exists(outFile)){
    N <- scorpion(tfMotifs = TF, gexMatrix = df, ppiNet = PPI, filterExpr = FALSE)[[1]]  #, nCores = 20
    N <- round(N, 3)
    gc()
    # i <- 1
    save(N, file = paste(outFile))
    }
}
)

plan(sequential)
