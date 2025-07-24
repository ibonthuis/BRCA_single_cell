packageVersion("Seurat") # 5.3.0 necessary (5.0.1 is not enough)
library("Seurat")

LN_METS_only <- get(load("/storage/kuijjerarea/ine/breast_met/single_cell/Data/LN_METS_plus_PRIMARY_scaledata_removed.RData"))

breastMETData <- get(load("/storage/kuijjerarea/ine/breast_met/single_cell/Data/GSE180286_SCTransform_integrated_BREAST_MET.RData"))


mets_merged <- merge(LN_METS_only, breastMETData)

save(mets_merged,
    file = "/storage/kuijjerarea/ine/breast_met/single_cell/Data/LN_METS_merged.RData")




save(
    LN_METS,
    "/storage/kuijjerarea/ine/breast_met/single_cell/Data/LN_METS_merged_processed.RData"
)
