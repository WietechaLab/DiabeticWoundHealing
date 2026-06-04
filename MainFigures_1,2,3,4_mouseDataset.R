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

# Figure 1B
DimPlot(MouseDataset, group.by = "Protein_maxID")

# Figure 1C
DimPlot(MouseDataset, group.by = "timepoint")

# Figure 1D
DimPlot(MouseDataset, group.by = "genotype")

# Figure 1E
DimPlot(MouseDataset, group.by = "major_cells", label = T)

# Figure 1F
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

# Figure 1G
table(MouseDataset$major_cells, MouseDataset$orig.ident)

# Figure 1H
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

# Figure 1I --- please refer to CellChat script

# Save the mouse dataset as RDS file
saveRDS(MouseDataset, "MouseDataset_final.rds")


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

# Figure 3A
DimPlot(dataset_macrophage, group.by = "timepoint")

# Figure 3B
DimPlot(dataset_macrophage, group.by = "genotype")

# Figure 3C
DimPlot(dataset_macrophage, group.by = "seurat_clusters", label = T)

# Figure 3D
plot1 <- FeatureScatter(dataset_macrophage, feature1 = "protein_F4_80", feature2 = "protein_Ly6C")
plot2 <- FeatureScatter(dataset_macrophage, feature1 = "protein_F4_80", feature2 = "protein_CD11b")
plot3 <- FeatureScatter(dataset_macrophage, feature1 = "protein_F4_80", feature2 = "protein_CD11c")
plot1 + plot2 + plot3

# Figure 3E
VlnPlot(dataset_macrophage, group.by = "seurat_clusters", pt.size = 0, features = c("protein_F4_80", "rna_Adgre1","protein_CD11b", "rna_Itgam","protein_Ly6C", "rna_Ly6c2","protein_CD11c", "rna_Itgax","protein_CD45"))
VlnPlot(dataset_macrophage, group.by = "seurat_clusters", pt.size = 0, features = c("rna_Vcan","rna_Pf4","rna_Il1b","rna_Cd74","rna_Lsp1"))

# Figure 3F
table(dataset_macrophage$orig.ident, dataset_macrophage$seurat_clusters)

# Figure 3G
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

# Figure 3H
DotPlot(dataset_macrophage, group.by="orig.ident", features = c("MacPolarization1","MacPolarization2","MacPolarization3")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

# Figure 3I
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

# Figure 3J
DotPlot(dataset_macrophage, group.by="orig.ident", features = c("Phagocytic_Mo1","Oxidative_Mo1","Inflammatory_Mo1","Remodeling_Mo1")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

# Figure 3K
DotPlot(dataset_macrophage, group.by="seurat_clusters", split.by = "genotype", cols = "Set3", features = c("Phagocytic_Mo1","Oxidative_Mo1","Inflammatory_Mo1","Remodeling_Mo1")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

# Save the macrophage dataset as RDS file
saveRDS(dataset_macrophage, "MouseDataset_macrophage_final.rds")


#### Figure 4 ####

# Read in the neutrophil dataset
dataset_neutrophil <- readRDS("MouseDataset_neutrophil_final.rds")

# List the metadata labels
names(dataset_neutrophil@meta.data)

# Figure 4A
DimPlot(dataset_neutrophil, group.by = "timepoint")

# Figure 4B
DimPlot(dataset_neutrophil, group.by = "genotype")

# Figure 4C
DimPlot(dataset_neutrophil, group.by = "seurat_clusters", label = T)

# Figure 4D
table(dataset_neutrophil$orig.ident, dataset_neutrophil$seurat_clusters)

# Figure 4E
neutrophil_markers <- FindAllMarkers(dataset_neutrophil, assay = "RNA", only.pos = TRUE, min.pct = 0.50, logfc.threshold = 0.50)

# Figure 4F
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

# Figure 4H
DotPlot(dataset_neutrophil, group.by = "orig.ident", features = c("N1_Neut1","N2_Neut1","scNeutMarkers_early1","scNeutMarkers_late1","scNeutMarkers_Arthritis1","scNeutMarkers_IL1_peritoneum1","scNeutMarkers_IL1_lung1")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

# Figure 4I
DotPlot(dataset_neutrophil, group.by = "seurat_clusters", split.by = "genotype",cols = "Set3", features = c("N1_Neut1","N2_Neut1","scNeutMarkers_early1","scNeutMarkers_late1","scNeutMarkers_Arthritis1","scNeutMarkers_IL1_peritoneum1","scNeutMarkers_IL1_lung1")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

# Save the neutrophil dataset as RDS file
saveRDS(dataset_neutrophil, "MouseDataset_neutrophil_final.rds")



### Human subtypes ####

MouseDataset_subtypes <- readRDS("MouseDataset_subtypes_final.rds")

# List the metadata
names(MouseDataset_subtypes@meta.data)

# List the annotated cell subtypes
Idents(MouseDataset_subtypes) <- "final_subtypes"
levels(MouseDataset_subtypes)


levels(MouseDataset_subtypes$orig.ident)


levels(MouseDataset_subtypes$timepoint)

# Show the breakdown of cell subtypes per wound type
table(MouseDataset_subtypes$final_subtypes, MouseDataset_subtypes$orig.ident)


### Human subtypes ####

HumanChronicWounds_subtypes <- readRDS("HumanChronicWounds_subtypes_final.rds")

JoinLayers(HumanChronicWounds_subtypes)

names(HumanChronicWounds_subtypes@meta.data)

Idents(HumanChronicWounds_subtypes) <- "final_subtypes"

# List the annotated cell subtypes
levels(HumanChronicWounds_subtypes$final_subtypes)

# List the annotated wound types
levels(HumanChronicWounds_subtypes$wound_type)

# Subset to chronic wounds and healthy skin only
Idents(HumanChronicWounds_subtypes) <- "wound_type"
HumanChronicWounds_subtypes <- subset(HumanChronicWounds_subtypes, idents = c("Diabetic foot ulcer","Non-diabetic foot ulcer","Healthy skin"))

# Subset to chronic wounds only
Idents(HumanChronicWounds_subtypes) <- "wound_type"
HumanChronicWounds_subtypes <- subset(HumanChronicWounds_subtypes, idents = c("Diabetic foot ulcer","Non-diabetic foot ulcer"))

# Show the breakdown of cell subtypes per wound type
table(HumanChronicWounds_subtypes$final_subtypes, HumanChronicWounds_subtypes$wound_type)

