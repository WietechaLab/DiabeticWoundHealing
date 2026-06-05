# General
library(ggplot2)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(ggpubr)
library(Seurat)
library(CellChat)

options(future.globals.maxSize = 1e24)


#### Figure 1 ####

# Read in full dataset
MouseDataset <- readRDS("MouseDataset_final.rds")

# List the metadata labels
names(MouseDataset@meta.data)

# Figure S1B
ElbowPlot(MouseDataset, ndims = 200)

# Figure S1C
DimPlot(MouseDataset, group.by = "seurat_clusters", label = T)

# Figure S1D
Idents(MouseDataset) <- "orig.ident"
p1 <- DimPlot(
  object = MouseDataset,     
  cells.highlight = WhichCells(MouseDataset, idents = "db_D3"), 
  cols.highlight = "red",     
  cols = "gray",
  sizes.highlight = 0.5,
  pt.size = 0.5    
)
p2 <- DimPlot(
  object = MouseDataset,     
  cells.highlight = WhichCells(MouseDataset, idents = "db_D6"), 
  cols.highlight = "red",     
  cols = "gray",
  sizes.highlight = 0.5,
  pt.size = 0.5    
)
p3 <- DimPlot(
  object = MouseDataset,     
  cells.highlight = WhichCells(MouseDataset, idents = "db_D10"), 
  cols.highlight = "red",     
  cols = "gray",
  sizes.highlight = 0.5,
  pt.size = 0.5    
)
p4 <- DimPlot(
  object = MouseDataset,     
  cells.highlight = WhichCells(MouseDataset, idents = "WT_D3"), 
  cols.highlight = "red",     
  cols = "gray",
  sizes.highlight = 0.5,
  pt.size = 0.5    
)
p5 <- DimPlot(
  object = MouseDataset,     
  cells.highlight = WhichCells(MouseDataset, idents = "WT_D6"), 
  cols.highlight = "red",     
  cols = "gray",
  sizes.highlight = 0.5,
  pt.size = 0.5    
)
p6 <- DimPlot(
  object = MouseDataset,     
  cells.highlight = WhichCells(MouseDataset, idents = "WT_D10"), 
  cols.highlight = "red",     
  cols = "gray",
  sizes.highlight = 0.5,
  pt.size = 0.5    
)
p1 + p2 + p3 + p4 + p5 + p6


## Figure 1B
DimPlot(MouseDataset, group.by = "Protein_maxID")

## Figure 1C
DimPlot(MouseDataset, group.by = "timepoint")

## Figure 1D
DimPlot(MouseDataset, group.by = "genotype")

# Figure S1E
DimPlot(MouseDataset, group.by = "Protein_maxID")
DimPlot(MouseDataset, group.by = "Protein_secondID")

# Figure S2A
AssayedCombined <- c("protein_CD115", "rna_Csf1r",
                     "protein_CD117", "rna_Kit",
                     "protein_CD11b", "rna_Itgam",
                     "protein_CD11c", "rna_Itgax",
                     "protein_CD140a", "rna_Pdgfra",
                     "protein_CD202b", "rna_Tek",
                     "protein_CD26", "rna_Dpp4",
                     "protein_CD3", "rna_Cd3e",
                     "protein_CD31", "rna_Pecam1",
                     "protein_CD326", "rna_Epcam",
                     "protein_CD335", "rna_Ncr1",
                     "protein_CD4", "rna_Cd4",
                     "protein_CD45", "rna_Ptprc",
                     "protein_CD8a", "rna_Cd8a",
                     "protein_F4_80", "rna_Adgre1",
                     "protein_Ly6C", "rna_Ly6c2",
                     "protein_Ly6G", "rna_Ly6g",
                     "protein_NK1.1", "rna_Klrb1c"
)
DotPlot(MouseDataset, group.by = "Protein_maxID", features = AssayedCombined) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")
DotPlot(MouseDataset, group.by = "Protein_secondID", features = AssayedCombined) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

# Figure S2B
DotPlot(MouseDataset, group.by = "seurat_clusters", features = AssayedCombined) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

## Figure 1E
DimPlot(MouseDataset, group.by = "major_cells", label = T)

# Figure S2G
VlnPlot(MouseDataset, group.by = "major_cells", features = "nCount_RNA", raster = FALSE, pt.size = 0)
VlnPlot(MouseDataset, group.by = "major_cells", features = "nFeature_RNA", raster = FALSE, pt.size = 0)

# Figure S5
# CD45
FeaturePlot(MouseDataset, c("protein_CD45", "rna_Ptprc"), cols = c("lightgrey", "darkblue"))
FeatureScatter(MouseDataset, feature1 = "protein_CD45", feature2 = "rna_Ptprc")
# CD11b
FeaturePlot(MouseDataset, c("protein_CD11b", "rna_Itgam"), cols = c("lightgrey", "darkblue"))
FeatureScatter(MouseDataset, feature1 = "protein_CD11b", feature2 = "rna_Itgam")
# F4_80
FeaturePlot(MouseDataset, c("protein_F4_80", "rna_Adgre1"), cols = c("lightgrey", "darkblue"))
FeatureScatter(MouseDataset, feature1 = "protein_F4_80", feature2 = "rna_Adgre1")
# Ly6C
FeaturePlot(MouseDataset, c("protein_Ly6C", "rna_Ly6c2"), cols = c("lightgrey", "darkblue"))
FeatureScatter(MouseDataset, feature1 = "protein_Ly6C", feature2 = "rna_Ly6c2")
# CD11c
FeaturePlot(MouseDataset, c("protein_CD11c", "rna_Itgax"), cols = c("lightgrey", "darkblue"))
FeatureScatter(MouseDataset, feature1 = "protein_CD11c", feature2 = "rna_Itgax")
# Ly6G
FeaturePlot(MouseDataset, c("protein_Ly6G", "rna_Ly6g"), cols = c("lightgrey", "darkblue"))
FeatureScatter(MouseDataset, feature1 = "protein_Ly6G", feature2 = "rna_Ly6g")
# CD31
FeaturePlot(MouseDataset, c("protein_CD31", "rna_Pecam1"), cols = c("lightgrey", "darkblue"))
FeatureScatter(MouseDataset, feature1 = "protein_CD31", feature2 = "rna_Pecam1")
# CD26
FeaturePlot(MouseDataset, c("protein_CD26", "rna_Dpp4"), cols = c("lightgrey", "darkblue"))
FeatureScatter(MouseDataset, feature1 = "protein_CD26", feature2 = "rna_Dpp4")
# CD140a
FeaturePlot(MouseDataset, c("protein_CD140a", "rna_Pdgfra"), cols = c("lightgrey", "darkblue"))
FeatureScatter(MouseDataset, feature1 = "protein_CD140a", feature2 = "rna_Pdgfra")
# CD3
FeaturePlot(MouseDataset, c("protein_CD3", "rna_Cd3e"), cols = c("lightgrey", "darkblue"))
FeatureScatter(MouseDataset, feature1 = "protein_CD3", feature2 = "rna_Cd3e")
# CD4
FeaturePlot(MouseDataset, c("protein_CD4", "rna_Cd4"), cols = c("lightgrey", "darkblue"))
FeatureScatter(MouseDataset, feature1 = "protein_CD4", feature2 = "rna_Cd4")
# NK1.1
FeaturePlot(MouseDataset, c("protein_NK1.1", "rna_Klrb1c"), cols = c("lightgrey", "darkblue"))
FeatureScatter(MouseDataset, feature1 = "protein_NK1.1", feature2 = "rna_Klrb1c")
# CD335
FeaturePlot(MouseDataset, c("protein_CD335", "rna_Ncr1"), cols = c("lightgrey", "darkblue"))
FeatureScatter(MouseDataset, feature1 = "protein_CD335", feature2 = "rna_Ncr1")
# CD202b
FeaturePlot(MouseDataset, c("protein_CD202b", "rna_Tek"), cols = c("lightgrey", "darkblue"))
FeatureScatter(MouseDataset, feature1 = "protein_CD202b", feature2 = "rna_Tek")
# CD117
FeaturePlot(MouseDataset, c("protein_CD117", "rna_Kit"), cols = c("lightgrey", "darkblue"))
FeatureScatter(MouseDataset, feature1 = "protein_CD117", feature2 = "rna_Kit")
# CD115
FeaturePlot(MouseDataset, c("protein_CD115", "rna_Csf1r"), cols = c("lightgrey", "darkblue"))
FeatureScatter(MouseDataset, feature1 = "protein_CD115", feature2 = "rna_Csf1r")
# CD326
FeaturePlot(MouseDataset, c("protein_CD326", "rna_Epcam"), cols = c("lightgrey", "darkblue"))
FeatureScatter(MouseDataset, feature1 = "protein_CD326", feature2 = "rna_Epcam")
# CD8a
FeaturePlot(MouseDataset, c("protein_CD8a", "rna_Cd8a"), cols = c("lightgrey", "darkblue"))
FeatureScatter(MouseDataset, feature1 = "protein_CD8a", feature2 = "rna_Cd8a")

## Figure 1F
FinalAssayed <- c("protein_CD140a","rna_Pdgfra",
                             "protein_CD26", "rna_Dpp4",
                             "protein_CD45", "rna_Ptprc",
                             "protein_CD11b", "rna_Itgam",
                             "protein_F4_80", "rna_Adgre1",
                             "protein_Ly6C", "rna_Ly6c2",
                             "protein_Ly6G", "rna_Ly6g", "rna_Csf3r",
                             "protein_CD31", "rna_Pecam1",
                             "protein_CD202b", "rna_Tek",
                             "rna_Acta2", "rna_Myh11",
                             "protein_CD3", "rna_Cd3e",
                             "protein_CD4", "rna_Cd4",
                             "rna_Gpx3", "rna_Gnas",
                             "rna_Des", "rna_Tpm2"
)
DotPlot(MouseDataset, group.by = "major_cells", features = FinalAssayed) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

## Figure 1G
table(MouseDataset$major_cells, MouseDataset$orig.ident)

## Figure 1H
DefaultAssay(MouseDataset) <- "RNA"
av.exp <- AggregateExpression(
  object = MouseDataset,
  assays = "RNA",
  features = NULL,
  return.seurat = FALSE,
  group.by = "orig.ident",
  add.ident = NULL,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  margin = 1,
  verbose = TRUE
)
av.exp <- as.data.frame(av.exp)
cor.exp <- cor(av.exp, method = "pearson")
Heatmap(cor.exp)

# Save the mouse dataset as RDS file
saveRDS(MouseDataset, "MouseDataset_final.rds")


#### CellChat Analysis of Mouse Dataset ####

# Subset the dataset to 3 DPW only for CellChat analysis
Idents(MouseDataset) <- "timepoint"
MouseDataset_3 <- subset(MouseDataset, idents = "3")
# Subset the dataset according to genotype
Idents(MouseDataset_3) <- "genotype"
MouseDataset_3_WT <- subset(MouseDataset_3, idents = "WT")
MouseDataset_3_db <- subset(MouseDataset_3, idents = "db")
# Set cell subtypes as the default identity for analysis of cell-cell interactions
Idents(MouseDataset_3_WT) <- "major_cells"
Idents(MouseDataset_3_db) <- "major_cells"
# Create CellChat Objects
cellchat_3_WT <- createCellChat(MouseDataset_3_WT, group.by = "major_cells")
cellchat_3_db <- createCellChat(MouseDataset_3_db, group.by = "major_cells")
# Set cell-cell interaction database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB 
# Set the used database in the object
cellchat_3_WT@DB <- CellChatDB.use
cellchat_3_db@DB <- CellChatDB.use
# Process the datasets for CellChat
cellchat_3_WT <- subsetData(cellchat_3_WT) 
cellchat_3_db <- subsetData(cellchat_3_db) 
# Run the standard CellChat workflow
future::plan("multisession", workers = 4)
cellchat_3_WT <- identifyOverExpressedGenes(cellchat_3_WT)
cellchat_3_db <- identifyOverExpressedGenes(cellchat_3_db)
cellchat_3_WT <- identifyOverExpressedInteractions(cellchat_3_WT)
cellchat_3_db <- identifyOverExpressedInteractions(cellchat_3_db)
cellchat_3_WT <- computeCommunProb(cellchat_3_WT, type = "triMean", population.size = TRUE) 
cellchat_3_db <- computeCommunProb(cellchat_3_db, type = "triMean", population.size = TRUE) 
cellchat_3_WT <- filterCommunication(cellchat_3_WT, min.cells = 10)
cellchat_3_db <- filterCommunication(cellchat_3_db, min.cells = 10)
cellchat_3_WT <- computeCommunProbPathway(cellchat_3_WT)
cellchat_3_db <- computeCommunProbPathway(cellchat_3_db)
cellchat_3_WT <- aggregateNet(cellchat_3_WT)
cellchat_3_db <- aggregateNet(cellchat_3_db)
cellchat_3_WT <- netAnalysis_computeCentrality(cellchat_3_WT, slot.name = "netP")
cellchat_3_db <- netAnalysis_computeCentrality(cellchat_3_db, slot.name = "netP")

# Save CellChat objects as RDS files
saveRDS(cellchat_3_WT, "CellChat_MajorCells_WT_3DPW.rds")
saveRDS(cellchat_3_db, "CellChat_MajorCells_db_3DPW.rds")
# Load CellChat objects from RDS files
cellchat_3_WT <- readRDS("CellChat_MajorCells_WT_3DPW.rds")
cellchat_3_db <- readRDS("CellChat_MajorCells_WT_3DPW.rds")

# Figure S6A
netAnalysis_signalingRole_scatter(cellchat_3_WT)
netAnalysis_signalingRole_scatter(cellchat_3_db)

# Figure S6B left
groupSize <- as.numeric(table(cellchat_3_WT@idents))
netVisual_circle(cellchat_3_WT@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_3_WT@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Figure S6B right
groupSize <- as.numeric(table(cellchat_3_db@idents))
netVisual_circle(cellchat_3_db@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_3_db@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Differential CellChat analysis
object.list_db_v_WT <- list(ND = cellchat_3_WT, db = cellchat_3_db)
cellchat_db_v_WT <- mergeCellChat(object.list_db_v_WT, add.names = names(object.list_db_v_WT))

# Figure S6C
gg1 <- compareInteractions(cellchat_db_v_WT, show.legend = F)
gg2 <- compareInteractions(cellchat_db_v_WT, show.legend = F, measure = "weight")
gg1 + gg2

# Figure S6D
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat_db_v_WT, weight.scale = T)
netVisual_diffInteraction(cellchat_db_v_WT, weight.scale = T, measure = "weight")

# Figure S6E
netVisual_heatmap(cellchat_db_v_WT)

## Figure 1I
netVisual_heatmap(cellchat_db_v_WT, measure = "weight")

# Figure S6F
netAnalysis_signalingChanges_scatter(cellchat_db_v_WT, idents.use = "Neutrophil", top.label = 1, color.use = "black")
netAnalysis_signalingChanges_scatter(cellchat_db_v_WT, idents.use = "Macrophage", top.label = 1, color.use = "black")
netAnalysis_signalingChanges_scatter(cellchat_db_v_WT, idents.use = "Fibroblast", top.label = 1, color.use = "black")
netAnalysis_signalingChanges_scatter(cellchat_db_v_WT, idents.use = "Smooth muscle", top.label = 1, color.use = "black")
netAnalysis_signalingChanges_scatter(cellchat_db_v_WT, idents.use = "Endothelial", top.label = 1, color.use = "black")


#### Cell subtype general analysis ####

# Subsetting the full dataset to fibroblasts, macrophages and neutrophils
Idents(MouseDataset) <- "major_cells"
dataset_fibroblast <- subset(MouseDataset, idents = "Fibroblast")
dataset_macrophage <- subset(MouseDataset, idents = "Macrophage")
dataset_neutrophil <- subset(MouseDataset, idents = "Neutrophil")

DefaultAssay(dataset_fibroblast) <- "RNA"
DefaultAssay(dataset_macrophage) <- "RNA"
DefaultAssay(dataset_neutrophil) <- "RNA"

# PCA
dataset_fibroblast <- RunPCA(dataset_fibroblast, verbose = TRUE)
dataset_macrophage <- RunPCA(dataset_macrophage, verbose = TRUE)
dataset_neutrophil <- RunPCA(dataset_neutrophil, verbose = TRUE)

# Figure S7A
ElbowPlot(dataset_fibroblast)

# Figure S7B
ElbowPlot(dataset_macrophage)

# Figure S7C
ElbowPlot(dataset_neutrophil)

# Seurat workflow to determine fibroblast subtypes
dataset_fibroblast <- RunUMAP(dataset_fibroblast, dims = 1:6)
dataset_fibroblast <- FindNeighbors(dataset_fibroblast, dims = 1:6)
dataset_fibroblast <- FindClusters(dataset_fibroblast, resolution = .1)

# Seurat workflow to determine macrophage subtypes
dataset_macrophage <- RunUMAP(dataset_macrophage, dims = 1:9)
dataset_macrophage <- FindNeighbors(dataset_macrophage, dims = 1:9)
dataset_macrophage <- FindClusters(dataset_macrophage, resolution = .1)

# Seurat workflow to determine neutrophil subtypes
dataset_neutrophil <- RunUMAP(dataset_neutrophil, dims = 1:5)
dataset_neutrophil <- FindNeighbors(dataset_neutrophil, dims = 1:5)
dataset_neutrophil <- FindClusters(dataset_neutrophil, resolution = .1)


#### Figure 2 ####

# Read in the fibroblast dataset
dataset_fibroblast <- readRDS("MouseDataset_fibroblast_final.rds")

# List the metadata labels
names(dataset_fibroblast@meta.data)

# Figure 2A
DimPlot(dataset_fibroblast, group.by = "timepoint")

# Figure 2B
DimPlot(dataset_fibroblast, group.by = "genotype")

# Figure 2C
DimPlot(dataset_fibroblast, group.by = "seurat_clusters", label = T)

# Figure 2D
table(dataset_fibroblast$orig.ident, dataset_fibroblast$seurat_clusters)

# Figure 2E
fibroblast_markers <- FindAllMarkers(dataset_fibroblast, assay = "RNA", only.pos = TRUE, min.pct = 0.50, logfc.threshold = 1)

# Figure 2F
VlnPlot(dataset_fibroblast, features =c("protein_CD140a","rna_Pdgfra","protein_CD26","rna_Dpp4","rna_Pi16","rna_Lrrc15","rna_Tnc","rna_Cilp","rna_Mgp"), pt.size = 0, stack = F, combine = T) 

# Extract from Table S18 and read in module scoring markers for mouse wound fibroblast subtypes
mmFB_Markers <- read_delim("mmFB_markers.txt", delim="\t", col_names = T)
# Perform module scoring
dataset_fibroblast <- AddModuleScore(
  object = dataset_fibroblast,
  features = mmFB_Markers[1],
  ctrl = 100,
  name = 'FIB_II_'
)
dataset_fibroblast <- AddModuleScore(
  object = dataset_fibroblast,
  features = mmFB_Markers[2],
  ctrl = 100,
  name = 'FIB_I_'
)
dataset_fibroblast <- AddModuleScore(
  object = dataset_fibroblast,
  features = mmFB_Markers[3],
  ctrl = 100,
  name = 'FIB_V_'
)
dataset_fibroblast <- AddModuleScore(
  object = dataset_fibroblast,
  features = mmFB_Markers[4],
  ctrl = 100,
  name = 'FIB_VII_'
)
dataset_fibroblast <- AddModuleScore(
  object = dataset_fibroblast,
  features = mmFB_Markers[5],
  ctrl = 100,
  name = 'FIB_IV_'
)
dataset_fibroblast <- AddModuleScore(
  object = dataset_fibroblast,
  features = mmFB_Markers[6],
  ctrl = 100,
  name = 'FIB_VIII_'
)
dataset_fibroblast <- AddModuleScore(
  object = dataset_fibroblast,
  features = mmFB_Markers[7],
  ctrl = 100,
  name = 'FIB_III_'
)
dataset_fibroblast <- AddModuleScore(
  object = dataset_fibroblast,
  features = mmFB_Markers[8],
  ctrl = 100,
  name = 'FIB_DP_'
)
dataset_fibroblast <- AddModuleScore(
  object = dataset_fibroblast,
  features = mmFB_Markers[9],
  ctrl = 100,
  name = 'FIB_IX_'
)
dataset_fibroblast <- AddModuleScore(
  object = dataset_fibroblast,
  features = mmFB_Markers[10],
  ctrl = 100,
  name = 'FIB_X_'
)
dataset_fibroblast <- AddModuleScore(
  object = dataset_fibroblast,
  features = mmFB_Markers[11],
  ctrl = 100,
  name = 'FIB_VI_'
)

# Figure 2H
DotPlot(dataset_fibroblast, group.by="orig.ident", features = c("FIB_I_1","FIB_II_1","FIB_III_1","FIB_IV_1","FIB_V_1","FIB_VI_1","FIB_VII_1","FIB_VIII_1","FIB_IX_1","FIB_X_1","FIB_DP_1")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

# Figure 2I
DotPlot(dataset_fibroblast, group.by="seurat_clusters", split.by = "genotype", cols = "Set3", features = c("FIB_I_1","FIB_II_1","FIB_III_1","FIB_IV_1","FIB_V_1","FIB_VI_1","FIB_VII_1","FIB_VIII_1","FIB_IX_1","FIB_X_1","FIB_DP_1")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

# Extract from Table S18 and read in module scoring markers for wound fibroblast phenotypes
wFB_pheno_Markers <- read_delim("wFB_pheno_Markers.txt", delim="\t", col_names = T)
# Perform module scoring
dataset_fibroblast <- AddModuleScore(
  object = dataset_fibroblast,
  features = wFB_pheno_Markers[1],
  ctrl = 100,
  name = 'CollagenScore'
)
dataset_fibroblast <- AddModuleScore(
  object = dataset_fibroblast,
  features = wFB_pheno_Markers[2],
  ctrl = 100,
  name = 'ElastinScore'
)
dataset_fibroblast <- AddModuleScore(
  object = dataset_fibroblast,
  features = wFB_pheno_Markers[3],
  ctrl = 100,
  name = 'ContractionScore'
)
dataset_fibroblast <- AddModuleScore(
  object = dataset_fibroblast,
  features = wFB_pheno_Markers[4],
  ctrl = 100,
  name = 'AdipogenesisScore'
)
dataset_fibroblast <- AddModuleScore(
  object = dataset_fibroblast,
  features = wFB_pheno_Markers[5],
  ctrl = 100,
  name = 'EarlyFBScore'
)
dataset_fibroblast <- AddModuleScore(
  object = dataset_fibroblast,
  features = wFB_pheno_Markers[6],
  ctrl = 100,
  name = 'LateFBScore'
)

# Figure 2J
DotPlot(dataset_fibroblast, group.by="orig.ident", features = c("CollagenScore1","ElastinScore1","ContractionScore1","AdipogenesisScore1","EarlyFBScore1","LateFBScore1")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

# Figure 2K
DotPlot(dataset_fibroblast, group.by="seurat_clusters", split.by = "genotype", cols = "Set3", features = c("CollagenScore1","ElastinScore1","ContractionScore1","AdipogenesisScore1","EarlyFBScore1","LateFBScore1")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

# Save the fibroblast dataset as RDS file
saveRDS(dataset_fibroblast, "MouseDataset_fibroblast_final.rds")


#### Figure 3 ####

# Read in the macrophage dataset
dataset_macrophage <- readRDS("MouseDataset_macrophage_final.rds")

# List the metadata labels
names(dataset_macrophage@meta.data)

## Figure 3A
DimPlot(dataset_macrophage, group.by = "timepoint")

## Figure 3B
DimPlot(dataset_macrophage, group.by = "genotype")

## Figure 3C
DimPlot(dataset_macrophage, group.by = "seurat_clusters", label = T)

## Figure 3D
plot1 <- FeatureScatter(dataset_macrophage, feature1 = "protein_F4_80", feature2 = "protein_Ly6C")
plot2 <- FeatureScatter(dataset_macrophage, feature1 = "protein_F4_80", feature2 = "protein_CD11b")
plot3 <- FeatureScatter(dataset_macrophage, feature1 = "protein_F4_80", feature2 = "protein_CD11c")
plot1 + plot2 + plot3

## Figure 3E
VlnPlot(dataset_macrophage, group.by = "seurat_clusters", pt.size = 0, features = c("protein_F4_80", "rna_Adgre1","protein_CD11b", "rna_Itgam","protein_Ly6C", "rna_Ly6c2","protein_CD11c", "rna_Itgax","protein_CD45"))
VlnPlot(dataset_macrophage, group.by = "seurat_clusters", pt.size = 0, features = c("rna_Vcan","rna_Pf4","rna_Il1b","rna_Cd74","rna_Lsp1"))

## Figure 3F
table(dataset_macrophage$orig.ident, dataset_macrophage$seurat_clusters)

## Figure 3G
dataset_macrophage <- FindAllMarkers(dataset_macrophage, assay = "RNA", only.pos = TRUE, min.pct = 0.50, logfc.threshold = 1)

# Extract from Table S18 and read in module scoring markers for macrophage signatues
CommonMoMarkers <- read_delim("CommonMoMarkers.txt", delim="\t", col_names = T)
mouse_M1_features <- CommonMoMarkers[1]
mouse_M2a_features <- CommonMoMarkers[2]
mouse_M2c_features <- CommonMoMarkers[3]
Polarization_list <- list(M1 = mouse_M1_features, 
                        M2a = mouse_M2a_features, 
                        M2c = mouse_M2c_features
                        )
# Perform module scoring
dataset_macrophage <- AddModuleScore(object = dataset_macrophage, features = Polarization_list, name = "MacPolarization")

## Figure 3H
DotPlot(dataset_macrophage, group.by="orig.ident", features = c("MacPolarization1","MacPolarization2","MacPolarization3")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

## Figure 3I
DotPlot(dataset_macrophage, group.by="seurat_clusters", split.by = "genotype", cols = "Set3", features = c("MacPolarization1","MacPolarization2","MacPolarization3")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

# Extract from Table S18 and read in module scoring markers for monocyte-derived macrophage phenotypes
CommonMoMarkers <- read_delim("CommonMoMarkers.txt", delim="\t", col_names = T)
Phagocytic_markers <- CommonMoMarkers[4]
Oxidative_markers <- CommonMoMarkers[5]
Inflammatory_markers <- CommonMoMarkers[6]
Remodeling_markers <- CommonMoMarkers[7]
# Perform module scoring
dataset_macrophage <- AddModuleScore(
  object = dataset_macrophage,
  assay = "RNA",
  features = Phagocytic_markers,
  ctrl = 100,
  name = 'Phagocytic_Mo'
)
dataset_macrophage <- AddModuleScore(
  object = dataset_macrophage,
  assay = "RNA",
  features = Oxidative_markers,
  ctrl = 100,
  name = 'Oxidative_Mo'
)
dataset_macrophage <- AddModuleScore(
  object = dataset_macrophage,
  assay = "RNA",
  features = Inflammatory_markers,
  ctrl = 100,
  name = 'Inflammatory_Mo'
)
dataset_macrophage <- AddModuleScore(
  object = dataset_macrophage,
  assay = "RNA",
  features = Remodeling_markers,
  ctrl = 100,
  name = 'Remodeling_Mo'
)

## Figure 3J
DotPlot(dataset_macrophage, group.by="orig.ident", features = c("Phagocytic_Mo1","Oxidative_Mo1","Inflammatory_Mo1","Remodeling_Mo1")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

## Figure 3K
DotPlot(dataset_macrophage, group.by="seurat_clusters", split.by = "genotype", cols = "Set3", features = c("Phagocytic_Mo1","Oxidative_Mo1","Inflammatory_Mo1","Remodeling_Mo1")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

# Save the macrophage dataset as RDS file
saveRDS(dataset_macrophage, "MouseDataset_macrophage_final.rds")


#### Figure 4 ####

# Read in the neutrophil dataset
dataset_neutrophil <- readRDS("MouseDataset_neutrophil_final.rds")

# List the metadata labels
names(dataset_neutrophil@meta.data)

## Figure 4A
DimPlot(dataset_neutrophil, group.by = "timepoint")

## Figure 4B
DimPlot(dataset_neutrophil, group.by = "genotype")

## Figure 4C
DimPlot(dataset_neutrophil, group.by = "seurat_clusters", label = T)

## Figure 4D
table(dataset_neutrophil$orig.ident, dataset_neutrophil$seurat_clusters)

## Figure 4E
neutrophil_markers <- FindAllMarkers(dataset_neutrophil, assay = "RNA", only.pos = TRUE, min.pct = 0.50, logfc.threshold = 0.50)

## Figure 4F
VlnPlot(dataset_neutrophil, features =c("protein_Ly6G","protein_CD11b","rna_Csf3r","rna_Fos","rna_Tnfaip3","rna_S100a8","rna_Spp1","rna_Sod2","rna_Vmp1"), pt.size = 0)

# Extract from Table S18 and read in module scoring markers for monocyte-derived macrophage phenotypes
CommonNeutMarkers <- read_delim("CommonNeutMarkers.txt", delim="\t", col_names = T)
CommonNeutMarkers_N1 <- CommonNeutMarkers[1]
CommonNeutMarkers_N2 <- CommonNeutMarkers[2]
scNeutMarkers_early <- CommonNeutMarkers[3]
scNeutMarkers_late <- CommonNeutMarkers[4]
scNeutMarkers_Arthritis <- CommonNeutMarkers[5]
scNeutMarkers_IL1_peritoneum <- CommonNeutMarkers[6]
scNeutMarkers_IL1_lung <- CommonNeutMarkers[7]
# Perform module scoring
dataset_neutrophil <- AddModuleScore(
  object = dataset_neutrophil,
  assay = "RNA",
  features = CommonNeutMarkers_N1,
  ctrl = 100,
  name = 'N1_Neut'
)
dataset_neutrophil <- AddModuleScore(
  object = dataset_neutrophil,
  assay = "RNA",
  features = CommonNeutMarkers_N2,
  ctrl = 100,
  name = 'N2_Neut'
)
dataset_neutrophil <- AddModuleScore(
  object = dataset_neutrophil,
  assay = "RNA",
  features = scNeutMarkers_early,
  ctrl = 100,
  name = 'scNeutMarkers_early'
)
dataset_neutrophil <- AddModuleScore(
  object = dataset_neutrophil,
  assay = "RNA",
  features = scNeutMarkers_late,
  ctrl = 100,
  name = 'scNeutMarkers_late'
)
dataset_neutrophil <- AddModuleScore(
  object = dataset_neutrophil,
  assay = "RNA",
  features = scNeutMarkers_Arthritis,
  ctrl = 100,
  name = 'scNeutMarkers_Arthritis'
)
dataset_neutrophil <- AddModuleScore(
  object = dataset_neutrophil,
  assay = "RNA",
  features = scNeutMarkers_IL1_peritoneum,
  ctrl = 50,
  name = 'scNeutMarkers_IL1_peritoneum'
)
dataset_neutrophil <- AddModuleScore(
  object = dataset_neutrophil,
  assay = "RNA",
  features = scNeutMarkers_IL1_lung,
  ctrl = 50,
  name = 'scNeutMarkers_IL1_lung'
)

## Figure 4H
DotPlot(dataset_neutrophil, group.by = "orig.ident", features = c("N1_Neut1","N2_Neut1","scNeutMarkers_early1","scNeutMarkers_late1","scNeutMarkers_Arthritis1","scNeutMarkers_IL1_peritoneum1","scNeutMarkers_IL1_lung1")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

## Figure 4I
DotPlot(dataset_neutrophil, group.by = "seurat_clusters", split.by = "genotype",cols = "Set3", features = c("N1_Neut1","N2_Neut1","scNeutMarkers_early1","scNeutMarkers_late1","scNeutMarkers_Arthritis1","scNeutMarkers_IL1_peritoneum1","scNeutMarkers_IL1_lung1")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

# Save the neutrophil dataset as RDS file
saveRDS(dataset_neutrophil, "MouseDataset_neutrophil_final.rds")


#### Figure 5 ####

# Naming of all fibroblast, macrophage and neutrophil subtypes prior to merging

# Naming of fibroblast subtypes
Idents(dataset_fibroblast) <- "seurat_clusters"
dataset_fibroblast[["final_subtypes"]] <- Idents(dataset_fibroblast)
Idents(dataset_fibroblast) <- "final_subtypes"
new.cluster.ids <- c("Fb_Pi16+_Dpp4+", "Fb_Lrrc15+_Tnc+", "Fb_Cilp+_Mgp+", "Fb_Prolif")
names(new.cluster.ids) <- levels(dataset_fibroblast)
dataset_fibroblast <- RenameIdents(dataset_fibroblast, new.cluster.ids)
dataset_fibroblast[["subtype"]] <- Idents(dataset_fibroblast)
Idents(dataset_fibroblast) <- "subtype"
levels(dataset_fibroblast)

# Naming of macrohage subtypes
Idents(dataset_macrophage) <- "seurat_clusters"
dataset_macrophage[["subtype"]] <- Idents(dataset_macrophage)
Idents(dataset_macrophage) <- "subtype"
new.cluster.ids <- c("Mo_F4/80+_CD11b+", "Mo_Ly6c+_CD11b+", "Mo_CD11c+_CD11b-", "Mo_Prolif")
names(new.cluster.ids) <- levels(dataset_macrophage)
dataset_macrophage <- RenameIdents(dataset_macrophage, new.cluster.ids)
dataset_macrophage[["subtype"]] <- Idents(dataset_macrophage)
Idents(dataset_macrophage) <- "subtype"
levels(dataset_macrophage)

# Naming of neutrophil subtypes
Idents(dataset_neutrophil) <- "seurat_clusters"
dataset_neutrophil[["subtype"]] <- Idents(dataset_neutrophil)
Idents(dataset_neutrophil) <- "subtype"
new.cluster.ids <- c("Neu_Csf3r+_Fos+", "Neu_Tnfaip3+Sod2+")
names(new.cluster.ids) <- levels(dataset_neutrophil)
dataset_neutrophil <- RenameIdents(dataset_neutrophil, new.cluster.ids)
dataset_neutrophil[["subtype"]] <- Idents(dataset_neutrophil)
Idents(dataset_neutrophil) <- "subtype"
levels(dataset_neutrophil)

# Merging of cell subtype datasets
MouseDataset_subtypes <- merge(dataset_neutrophil, y = c(dataset_macrophage, dataset_fibroblast))
# OR load in the saved RDS file
MouseDataset_subtypes <- readRDS("MouseDataset_subtypes_final.rds")

# Reorder cell subtypes
unique(MouseDataset_subtypes@meta.data[["subtype"]])
mm_levels <- c("Neu_Tnfaip3+Sod2+", "Neu_Csf3r+_Fos+", "Mo_Prolif", "Mo_F4/80+_CD11b+", "Mo_Ly6c+_CD11b+", "Mo_CD11c+_CD11b-", "Fb_Prolif", "Fb_Pi16+_Dpp4+", "Fb_Lrrc15+_Tnc+", "Fb_Cilp+_Mgp+")
MouseDataset_subtypes@meta.data[["subtype"]] <- factor(MouseDataset_subtypes@meta.data[["subtype"]], levels = mm_levels)
levels(MouseDataset_subtypes@meta.data[["subtype"]])

## Figure 5A
DotPlot(MouseDataset_subtypes, group.by = "subtype",
        features = c("protein_CD45","rna_Ptprc","protein_Ly6G", "rna_Ly6g","Tnfaip3","Sod2","Csf3r","Fos","protein_CD11b","rna_Itgam","protein_F4-80","rna_Adgre1","protein_Ly6C", "rna_Ly6c2","protein_CD11c", "rna_Itgax","protein_CD140a", "rna_Pdgfra","protein_CD26", "rna_Dpp4","Pi16","Lrrc15","Tnc","Cilp","Mgp")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

# Save the merged subtype dataset as RDS file
saveRDS(MouseDataset_subtypes, "MouseDataset_subtypes_final.rds")

# Subset the merged dataset to 3 DPW only for CellChat analysis
Idents(MouseDataset_subtypes) <- "timepoint"
MouseDataset_subtypes_3 <- subset(MouseDataset_subtypes, idents = "3")
# Subset the merged dataset according to genotype
Idents(MouseDataset_subtypes_3) <- "genotype"
MouseDataset_subtypes_WT <- subset(MouseDataset_subtypes_3, idents = "WT")
MouseDataset_subtypes_db <- subset(MouseDataset_subtypes_3, idents = "db")
# Set cell subtypes as the default identity for analysis of cell-cell interactions
Idents(MouseDataset_subtypes_WT) <- "subtype"
Idents(MouseDataset_subtypes_db) <- "subtype"
# Create CellChat Objects
cellchat_3_WT <- createCellChat(MouseDataset_subtypes_WT, group.by = "subtype")
cellchat_3_db <- createCellChat(MouseDataset_subtypes_db, group.by = "subtype")
# Set cell-cell interaction database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB 
# Set the used database in the object
cellchat_3_WT@DB <- CellChatDB.use
cellchat_3_db@DB <- CellChatDB.use
# Process the datasets for CellChat
cellchat_3_WT <- subsetData(cellchat_3_WT) 
cellchat_3_db <- subsetData(cellchat_3_db) 
# Run the standard CellChat workflow
future::plan("multisession", workers = 4)
cellchat_3_WT <- identifyOverExpressedGenes(cellchat_3_WT)
cellchat_3_db <- identifyOverExpressedGenes(cellchat_3_db)
cellchat_3_WT <- identifyOverExpressedInteractions(cellchat_3_WT)
cellchat_3_db <- identifyOverExpressedInteractions(cellchat_3_db)
cellchat_3_WT <- computeCommunProb(cellchat_3_WT, type = "triMean", population.size = TRUE) 
cellchat_3_db <- computeCommunProb(cellchat_3_db, type = "triMean", population.size = TRUE) 
cellchat_3_WT <- filterCommunication(cellchat_3_WT, min.cells = 10)
cellchat_3_db <- filterCommunication(cellchat_3_db, min.cells = 10)
cellchat_3_WT <- computeCommunProbPathway(cellchat_3_WT)
cellchat_3_db <- computeCommunProbPathway(cellchat_3_db)
cellchat_3_WT <- aggregateNet(cellchat_3_WT)
cellchat_3_db <- aggregateNet(cellchat_3_db)
cellchat_3_WT <- netAnalysis_computeCentrality(cellchat_3_WT, slot.name = "netP")
cellchat_3_db <- netAnalysis_computeCentrality(cellchat_3_db, slot.name = "netP")

# Save CellChat objects as RDS files
saveRDS(cellchat_3_WT, "CellChat_subtypes_WT_3DPW.rds")
saveRDS(cellchat_3_db, "CellChat_subtypes_db_3DPW.rds")
# Load CellChat objects from RDS files
cellchat_3_WT <- readRDS("CellChat_subtypes_WT_3DPW.rds")
cellchat_3_db <- readRDS("CellChat_subtypes_db_3DPW.rds")

# Figure S8A
netAnalysis_signalingRole_scatter(cellchat_3_WT)
netAnalysis_signalingRole_scatter(cellchat_3_db)

# Figure S8B left
groupSize <- as.numeric(table(cellchat_3_WT@idents))
netVisual_circle(cellchat_3_WT@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_3_WT@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Figure S8B right
groupSize <- as.numeric(table(cellchat_3_db@idents))
netVisual_circle(cellchat_3_db@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_3_db@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Differential CellChat analysis
object.list_db_v_WT <- list(ND = cellchat_3_WT, db = cellchat_3_db)
cellchat_db_v_WT <- mergeCellChat(object.list_db_v_WT, add.names = names(object.list_db_v_WT))

# Figure S8C
gg1 <- compareInteractions(cellchat_db_v_WT, show.legend = F)
gg2 <- compareInteractions(cellchat_db_v_WT, show.legend = F, measure = "weight")
gg1 + gg2

# Figure S8D
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat_db_v_WT, weight.scale = T)
netVisual_diffInteraction(cellchat_db_v_WT, weight.scale = T, measure = "weight")

## Figure 5B
netVisual_heatmap(cellchat_db_v_WT)

## Figure 5C
netVisual_heatmap(cellchat_db_v_WT, measure = "weight")

## Figure 5D
netAnalysis_signalingChanges_scatter(cellchat_db_v_WT, idents.use = "Neu_Tnfaip3+Sod2+", top.label = 1, color.use = "black")

## Figure 5E
netAnalysis_signalingChanges_scatter(cellchat_db_v_WT, idents.use = "Mo_F4/80+_CD11b+", top.label = 1, color.use = "black")

## Figure 5F
netAnalysis_signalingChanges_scatter(cellchat_db_v_WT, idents.use = "Mo_Ly6c+_CD11b+", top.label = 1, color.use = "black")

# Figure S8E
netAnalysis_signalingChanges_scatter(cellchat_db_v_WT, idents.use = "Neu_Csf3r+_Fos+", top.label = 1, color.use = "black")
netAnalysis_signalingChanges_scatter(cellchat_db_v_WT, idents.use = "Mo_Prolif", top.label = 1, color.use = "black")
netAnalysis_signalingChanges_scatter(cellchat_db_v_WT, idents.use = "Mo_CD11c+_CD11b-", top.label = 1, color.use = "black")

# Figure S8F
netAnalysis_signalingChanges_scatter(cellchat_db_v_WT, idents.use = "Fb_Pi16+_Dpp4+", top.label = 1, color.use = "black")
netAnalysis_signalingChanges_scatter(cellchat_db_v_WT, idents.use = "Fb_Lrrc15+_Tnc+", top.label = 1, color.use = "black")
netAnalysis_signalingChanges_scatter(cellchat_db_v_WT, idents.use = "Fb_Cilp+_Mgp+", top.label = 1, color.use = "black")
netAnalysis_signalingChanges_scatter(cellchat_db_v_WT, idents.use = "Fb_Prolif", top.label = 1, color.use = "black")

# Code to perform differential expression analysis of ligands and receptors in db vs ND
pos.dataset = "db"
cellchat_db_v_WT@var.features$features
cellchat_db_v_WT@var.features[[features.name]]
features.name = paste0(pos.dataset, ".merged")
cellchat_db_v_WT <- identifyOverExpressedGenes(cellchat_db_v_WT, group.dataset = "genotype", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 0.05) 
net <- netMappingDEG(cellchat_db_v_WT, features.name = features.name, thresh = 0.05)
write.csv(net, file = file.path(getwd(), "CellChat_mouse_DEG_LRpairs.txt"))
net.up <- subsetCommunication(cellchat_db_v_WT, net = net, datasets = "db", ligand.logFC = 0.20, receptor.logFC = 0.20)
net.down <- subsetCommunication(cellchat_db_v_WT, net = net, datasets = "db", ligand.logFC = -0.20, receptor.logFC = -0.20)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat_db_v_WT)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_db_v_WT)
pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]
sort.pairLR.use.up <- arrange(pairLR.use.up, interaction_name)
sort.pairLR.use.down <- arrange(pairLR.use.down, interaction_name)

## Figure 5G
netVisual_bubble(cellchat_db_v_WT, pairLR.use = pairLR.use.up, sources.use = "Neu_Tnfaip3+Sod2+", comparison = 2,  angle.x = 90, remove.isolate = T, title.name = paste0("Up-regulated signaling in ", names(object.list_db_v_WT)[2]))

## Figure 5H
netVisual_bubble(cellchat_db_v_WT, pairLR.use = pairLR.use.up, sources.use = "Mo_F4/80+_CD11b+", comparison = 2,  angle.x = 90, remove.isolate = T, title.name = paste0("Up-regulated signaling in ", names(object.list_db_v_WT)[2]))

## Figure 5I
netVisual_bubble(cellchat_db_v_WT, pairLR.use = pairLR.use.up, sources.use = "Mo_Ly6c+_CD11b+", comparison = 2,  angle.x = 90, remove.isolate = T, title.name = paste0("Up-regulated signaling in ", names(object.list_db_v_WT)[2]))

# CD44 pathway expression in cell subtypes
Cd44_pathway <- c("Cd44","Spp1","Cxcl2","Cxcr2","Icam1","Itgal","Itgam","Itgb2","Vcam1","Fn1","Col1a1")
Idents(MouseDataset_subtypes_3) <- "subtype"
## Figure 5J bottom
DotPlot(MouseDataset_subtypes_3, ident = c("Neu_Tnfaip3+Sod2+","Neu_Csf3r+_Fos+"), group.by = "subtype",split.by = "genotype", cols = "Set3", 
        features = Cd44_pathway) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")
## Figure 5J middle
DotPlot(MouseDataset_subtypes_3, ident = c("Mo_F4/80+_CD11b+","Mo_Ly6c+_CD11b+","Mo_CD11c+_CD11b-"), group.by = "subtype",split.by = "genotype", cols = "Set3", 
        features = Cd44_pathway) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")
## Figure 5J top
DotPlot(MouseDataset_subtypes_3, ident = c("Fb_Pi16+_Dpp4+","Fb_Lrrc15+_Tnc+","Fb_Cilp+_Mgp+"), group.by = "subtype",split.by = "genotype", cols = "Set3", 
        features = Cd44_pathway) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")


#### Figure 6 ####

## Figure 6C
Idents(MouseDataset_subtypes) <- "subtype"
p <- DotPlot(MouseDataset_subtypes, ident = c("Mo_F4/80+_CD11b+","Mo_Ly6c+_CD11b+","Mo_CD11c+_CD11b-"), group.by = "subtype", split.by = "orig.ident", cols = "Set3", 
             features = c("MacPolarization1","MacPolarization2","MacPolarization3","Cd44","Spp1","Icam1","S100a9","Ncf4","Sod2","Il1b","Tnf","Mrc1","Cd36")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")
p$data$id <- factor(x = p$data$id, levels = c("Mo_F4/80+_CD11b+_db_D3", "Mo_F4/80+_CD11b+_db_D6", "Mo_F4/80+_CD11b+_db_D10", 
                                              "Mo_F4/80+_CD11b+_WT_D3", "Mo_F4/80+_CD11b+_WT_D6", "Mo_F4/80+_CD11b+_WT_D10",
                                              "Mo_Ly6c+_CD11b+_db_D3", "Mo_Ly6c+_CD11b+_db_D6", "Mo_Ly6c+_CD11b+_db_D10", 
                                              "Mo_Ly6c+_CD11b+_WT_D3", "Mo_Ly6c+_CD11b+_WT_D6", "Mo_Ly6c+_CD11b+_WT_D10",
                                              "Mo_CD11c+_CD11b-_db_D3", "Mo_CD11c+_CD11b-_db_D6", "Mo_CD11c+_CD11b-_db_D10",
                                              "Mo_CD11c+_CD11b-_WT_D3", "Mo_CD11c+_CD11b-_WT_D6", "Mo_CD11c+_CD11b-_WT_D10"))
print(p)


#### Figure 7 ####

## Figure 7J
table(MouseDataset_subtypes$subtype, MouseDataset_subtypes$orig.ident)
