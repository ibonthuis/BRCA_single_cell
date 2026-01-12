library(Seurat)
library("stringr")
#packageVersion("Seurat") # v 5.0.1




source("create_seurat_object_fn.R")

matricesList <- list.files(path = 'raw/GSE180286_extracted', pattern = '.txt', full.names = TRUE)
breastMETData1 <- readSample(matricesList[1])
breastMET_list <- lapply(matricesList[-1], readSample)

breastMETData <-  merge(breastMETData1, y = breastMET_list, project = "BRCAMET")
# all_combined_join_layers <- JoinLayers(breastMETData)
# breastMETData <- scQC(all_combined_join_layers)



breastMETData <- NormalizeData(breastMETData)
breastMETData <- FindVariableFeatures(breastMETData) # This step takes a very long time.
str(breastMETData)
breastMETData <- ScaleData(breastMETData)
breastMETData <- RunPCA(breastMETData)
#breastMETData <- RunHarmony(breastMETData, group.by.vars = 'orig.ident', max.iter.harmony = 100)

breastMETData <- RunUMAP(breastMETData, reduction = 'pca', dims = 1:20) 
# Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
#To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'

DimPlot(breastMETData)

save(
    breastMETData,
    file = "Data/GSE180286_BREAST_MET.RData"
)

# pdf("Figures/first_umap_gse180286.pdf")
# UMAPPlot(breastMETData)
# dev.off()


breastMETData <- SCTransform(breastMETData) # duurde best een tijdje
breastMETData <- scQC(breastMETData) 


save(
    breastMETData,
    file = "Data/GSE180286_SCTransform_integrated_BREAST_MET.RData"
)

