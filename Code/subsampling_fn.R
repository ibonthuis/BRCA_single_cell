# Note, maybe I can remove the print statements.
sample_cells_recursively <- function(nr_of_cells_to_sample, df) {
  # Initialize the list to store sampled dataframes
  sampled_dfs <- list()
  
  # While there are enough columns to sample
   number_of_cols <- ncol(df)
   print(number_of_cols)
  while (ncol(df) >= nr_of_cells_to_sample) {
    # Sample columns
    idx <- sample(ncol(df), nr_of_cells_to_sample)
    sampled_df <- df[, idx]
    sampled_dfs <- append(sampled_dfs, list(sampled_df))
    
    # Remove sampled columns from the original dataframe
    if(ncol(df) == 101){
      print("broke")
      break
    } else {
        df <- df[, -idx]
        print("done")
    }
    
  } 
    return(sampled_dfs)
}


subsample <- function(ID, ALL, count_matrix){
  # ALL = SEURATobject
  subsamples <- list()
  outFile_subsample <- paste0('Results/Subsamples/all_subsamples')
  if(!file.exists(outFile_subsample)){
    idData <- count_matrix[, ALL$patient_sampleID %in% ID]
   # colnames(idData) <- names(colnames(idData))
    if(ncol(idData) >= 100){
      list_dfs <- sample_cells_recursively(100, idData)
#length(list_dfs)
      for(i in 1:length(list_dfs)){
        print(i)
       # print(head(list_dfs[[i]]))
        #names(list_dfs)[i] <- paste0(ID, "_s", i)
        subsamples <- append(subsamples, list_dfs[i])
        names(subsamples)[length(subsamples)] <- paste0(ID, "_s", i)
      }
       # length(subsamples)
        #save(subsamples, file = paste(outFile_subsample, ".RData"))

    }
  } 
  return(subsamples)
}
