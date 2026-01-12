required_libraries <- c(
    "data.table",    
    "dplyr",
    "Seurat",
    "harmony",
    "Nebulosa" #BiocManager::install("Nebulosa")
    )
    
for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
}

LN_METS <- get(load("Data/LN_METS_merged.RData")) 

LN_METS_only <- get(load("Data/LN_METS_plus_PRIMARY_scaledata_removed.RData"))
breastMETData <- get(load("Data/GSE180286_SCTransform_integrated_BREAST_MET.RData"))

# addition since August 6th
LN_METS <- get(load("Data/LN_METS_merged_prepped_qc.RData"))
unique(LN_METS$celltype)

var_features_all <- c(VariableFeatures(breastMETData), VariableFeatures(LN_METS_only))
var_features_all <- unique(var_features_all) #4628 features 

VariableFeatures(LN_METS) <- var_features_all
LN_METS <- NormalizeData(LN_METS)
#LN_METS <- FindVariableFeatures(LN_METS) # This step takes a very long time.
LN_METS <- ScaleData(LN_METS)
LN_METS <- RunPCA(LN_METS)

LN_METS <- RunHarmony(LN_METS, "orig.ident")
LN_METS <- RunUMAP(LN_METS, reduction = 'harmony', dims = 1:50)
LN_METS <- FindNeighbors(LN_METS, reduction = 'harmony', dims = 1:50)
LN_METS <- FindClusters(LN_METS)

Seurat::DimPlot(LN_METS, reduction = "harmony", group.by = "orig.ident")
Seurat::DimPlot(LN_METS, reduction = "umap", group.by = "new_clusters")
DimPlot(LN_METS, reduction = "harmony")

markerList <- c('EPCAM', 'CDH1', 'COL1A1', 'COL3A1', 'MS4A1', 
                'CDH5', 'PECAM1', 'S100B', 'CDH2', 'PTPRC', 
                'CD3E', 'CD14', 'IL1RL1', 'MZB1', 'MKI67', 'ADIPOQ')

options(bitmapType='cairo')

png('Figures/atlasMarkers_adipo.png', width = 4000*.85, height = 3000*.85, res = 300)
plot_density(LN_METS, markerList) & 
  theme_light() & 
  theme(plot.title = element_text(face = 4), legend.key.width = unit(0.2, 'cm')) &
  scale_color_gradient(low = '#e5e5e5', high = '#d90429')
dev.off()


pdf("Figures/merged_LN_MET_featuresplot_extra.pdf", h = 12, w = 10)
FeaturePlot(LN_METS, c('EPCAM', 'CDH1', 'COL1A1', 'COL3A1', 'MS4A1', 'MZB1', 'MKI67',
                'CDH5', 'PECAM1', 'S100B', 'CDH2', 'PTPRC', 
                'CD3E', 'CD14', 'IL1RL1',  'PLD4',
                'CD8A', 'NKG7', 'FCER1A', 'CST3', 'CCR7', 'NRP1', 'CLEC4C' )) # CLEC4C (also BCDA-2), CCR7 and NRP1 mark plasmacytoid dendritic cells
dev.off()


pdf("Figures/merged_LN_MET_first_umap.pdf") 
DimPlot(LN_METS, reduction = "umap", group.by = "orig.ident")
dev.off()

pdf("Figures/merged_LN_MET_umap_new_clusters.pdf") 
DimPlot(LN_METS, reduction = "umap", group.by = "new_clusters")
dev.off()

pdf("Figures/merged_LN_MET_umap_seurat_clusters.pdf") 
DimPlot(LN_METS, label = TRUE)
dev.off()

pdf("Figures/merged_LN_MET_umap_seurat_clusters_harmony.pdf") 
DimPlot(LN_METS, reduction = "harmony")
dev.off()



save(
    LN_METS,
    file = "Data/LN_METS_merged.RData"
)


newID <- rep(NA, 30)
newID[c(4)] <- "Endothelial"
newID[c(2, 8, 19, 21, 28)] <- "B"
newID[c(5, 15, 20, 25, 22, 16)] <- "Epithelial"
newID[c(9, 12, 26, 27, 29)] <- "Fibroblast"
newID[c(11, 23)] <- "Myeloid"
newID[c(14)] <- "PlasmaB"
newID[c(10)] <- "Plasmablast"
newID[c(0, 1, 3, 6, 7, 13, 17, 18)] <- "T"
newID[c(24)] <- "PlasmacytoidDendritic"
unique(newID)
newID <- c("T", newID)
newID <- newID[c(1:30)]
newID[c(2)] <- "CD8/NK"

Idents(LN_METS) <- LN_METS$seurat_clusters
levels(Idents(LN_METS)) <- newID
LN_METS$celltype <- Idents(LN_METS)
#UMAPPlot(LN_METS, label = TRUE)

pdf("Figures/merged_LN_MET_umap_celltype_new_annotation_labelled.pdf") 
DimPlot(LN_METS, reduction = "umap", group.by = "celltype", label = TRUE)
dev.off()


metadata$patient_name <- paste(metadata$patient_ID, metadata$GroupID, metadata$celltype, sep = "_")
  print(unique(metadata$patient_name))
numberofcells <- as.data.frame(table(metadata$patient_name))

metadata <- LN_METS@meta.data
head(metadata)

