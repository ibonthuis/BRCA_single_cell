library("Seurat")
library("dplyr")
library("BiocManager")
library("glmGamPoi")
#BiocManager::install("glmGamPoi")

sc_breastmet <- readRDS("raw/BRCA.rds")





# In line with prior workflows, you can also into split your object into a list of multiple objects based on a metadata
# column creates a list of two objects
sc_breastmet <- SplitObject(sc_breastmet, split.by = "FlatFormID")

LN_METS_only <- merge(sc_breastmet$GSE167036, list(sc_breastmet$GSE195861, sc_breastmet$GSE225600))

LN_METS_only <- SCTransform(LN_METS_only) %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters() %>%
    RunUMAP(dims = 1:30)

save(
    LN_METS_only,
    file = "Data/LN_METS_plus_PRIMARY.RData"
)



LN_METS_only[["SCT"]]$scale.data <- NULL

save(
  LN_METS_only, 
  file = "Data/LN_METS_plus_PRIMARY_scaledata_removed.RData"
)


PatientID <- rep(paste("patient", seq(1,18), sep=""))

# Create a dataframe with sample pairs
sample_pairs <- data.frame(
Sample1 = c("CA1", "CA2", "CA3", "CA4", "CA5", "CA6", "CA7", "CA8", "GSM5852276_P2", "GSM5852278_P1", "GSM5852280_P6", "GSM5852282_P5", "GSM5852284_P4", "GSM5852286_P3", "L2", "L3", "L6", "L7", ""),
  Sample2 = c("LN1", "LN2", "LN3", "LN4", "LN5", "LN6", "LN7", "LN8", "GSM5852277_M2", "GSM5852279_M1", "GSM5852281_M6", "GSM5852283_M5", "GSM5852285_M4", "GSM5852287_M3", "T2", "T3", "T6", "T7"),
  PatientID = PatientID
)

LN_METS_only$patient_ID <- NA

# Assign patient IDs based on the sample pairs
for (i in seq_len(nrow(sample_pairs))) {
  LN_METS_only$patient_ID[grepl(sample_pairs$Sample1[i], LN_METS_only$orig.ident) |
                          grepl(sample_pairs$Sample2[i], LN_METS_only$orig.ident)] <- sample_pairs$PatientID[i]
}

# Create the patient_sampleID column by concatenating orig.ident, new_clusters, and patient_ID
LN_METS_only$patient_sampleID <- paste0(LN_METS_only$orig.ident, '-', LN_METS_only$new_clusters, '-', LN_METS_only$patient_ID)

save(
  LN_METS_only, 
  file = "Data/LN_METS_plus_PRIMARY_scaledata_removed.RData"
)

fwrite(ALL@meta.data, "metadata_ALL_annotated.tsv", sep = "\t", row.names = TRUE)