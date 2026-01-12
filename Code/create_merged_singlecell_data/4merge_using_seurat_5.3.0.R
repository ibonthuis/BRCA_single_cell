packageVersion("Seurat") # 5.3.0 necessary (5.0.1 is not enough)
library("Seurat")
packageVersion("sp")



LN_METS_only <- get(load("Data/LN_METS_plus_PRIMARY_scaledata_removed.RData"))

breastMETData <- get(load("Data/GSE180286_SCTransform_integrated_BREAST_MET.RData"))


mets_merged <- merge(LN_METS_only, breastMETData)

save(mets_merged,
    file = "Data/LN_METS_merged.RData")




save(
    LN_METS,
    "Data/LN_METS_merged_processed.RData"
)
