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

#### Figure 7 ####

# Code for re-analysis of human datasets.
# For this code to work, the datasets must be downloaded from the relevant GEO accessions and placed into respective named folders

#### GSE223964 ####

GSE223964_data_dir <- file.path(getwd(),"GSE223964")
GSE223964_data_samples <- list.files(GSE223964_data_dir, full.names=TRUE, recursive=FALSE)
GSE223964_data_names <- list.files(GSE223964_data_dir, full.names=FALSE, recursive=FALSE)
listSeurat <- list()
for (i in 1:length(GSE223964_data_samples)) {
  sample <- Read10X_h5(GSE223964_data_samples[i], use.names = TRUE, unique.features = TRUE)
  seuratObject = toString(GSE223964_data_names[i])
  seuratObject <- CreateSeuratObject(counts = sample)
  seuratObject@meta.data$sampleName <- GSE223964_data_names[i]
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^MT-", col.name = "percent.mt")
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^RPL", col.name = "percent.rpl")
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^RPS", col.name = "percent.rps")
  listSeurat[i] <- seuratObject
}

j = 9 # cycle from j = 1 to 9
listSeurat[[j]]
VlnPlot(listSeurat[[j]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rpl", "percent.rps"), ncol = 5)
FeatureScatter(listSeurat[[j]], feature1 = "percent.mt", feature2 = "nCount_RNA")
FeatureScatter(listSeurat[[j]], feature1 = "percent.rps", feature2 = "nFeature_RNA")
FeatureScatter(listSeurat[[j]], feature1 = "percent.rpl", feature2 = "nCount_RNA")
listSeurat[[j]] <- subset(listSeurat[[j]], subset = nFeature_RNA > 100 & nCount_RNA < 75000 & percent.mt < 15 & percent.rps < 30 & percent.rpl < 40)
DefaultAssay(listSeurat[[j]]) <- "RNA"
sce <- scDblFinder(GetAssayData(listSeurat[[j]], assay="RNA", slot="counts"), samples=Idents(listSeurat[[j]]))
listSeurat[[j]]$scDblFinder.score <- sce$scDblFinder.score
VlnPlot(listSeurat[[j]], features = "scDblFinder.score", raster=FALSE, pt.size=0.5)
listSeurat[[j]] <- subset(listSeurat[[j]], scDblFinder.score < 0.5)

GSE223964.merged <- merge(x = listSeurat[[1]], y = list(listSeurat[[2]], listSeurat[[3]], listSeurat[[4]], listSeurat[[5]], 
                                                        listSeurat[[6]], listSeurat[[7]], listSeurat[[8]]),
                          project = "GSE223964_merged", merge.data = TRUE)
DefaultAssay(GSE223964.merged) <- 'RNA'


#### GSE265972 ####

GSE265972_data_dir <- file.path(getwd(),"GSE265972")
GSE265972_data_dirs <- list.dirs(GSE265972_data_dir, full.names=TRUE, recursive=FALSE)
GSE265972_data_samples <- list.dirs(GSE265972_data_dir, full.names=FALSE, recursive=FALSE)
listSeurat <- list()
for (i in 1:length(GSE265972_data_dirs)) {
  sample <- Read10X_GEO(data_dir = GSE265972_data_dirs[i], gene.column = 2)
  seuratObject = toString(GSE265972_data_samples[i])
  seuratObject <- CreateSeuratObject(counts = sample)
  seuratObject@meta.data$sampleName <- GSE265972_data_samples[i]
  seuratObject@meta.data$sampleType <- substr(GSE265972_data_samples[i], start = 1, stop = 2)
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^MT-", col.name = "percent.mt")
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^RPL", col.name = "percent.rpl")
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^RPS", col.name = "percent.rps")
  listSeurat[i] <- seuratObject
}

j = 9 # cycle from j = 1 to 9
listSeurat[[j]]
VlnPlot(listSeurat[[j]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rpl", "percent.rps"), ncol = 5)
FeatureScatter(listSeurat[[j]], feature1 = "percent.mt", feature2 = "nCount_RNA")
FeatureScatter(listSeurat[[j]], feature1 = "percent.rps", feature2 = "nFeature_RNA")
FeatureScatter(listSeurat[[j]], feature1 = "percent.rpl", feature2 = "nCount_RNA")
listSeurat[[j]] <- subset(listSeurat[[j]], subset = nFeature_RNA > 50 & nCount_RNA < 75000 & percent.mt < 20 & percent.rps < 20 & percent.rpl < 30)
DefaultAssay(listSeurat[[j]]) <- "RNA"
sce <- scDblFinder(GetAssayData(listSeurat[[j]], assay="RNA", slot="counts"), samples=Idents(listSeurat[[j]]))
listSeurat[[j]]$scDblFinder.score <- sce$scDblFinder.score
VlnPlot(listSeurat[[j]], features = "scDblFinder.score", raster=FALSE, pt.size=0.5)
listSeurat[[j]] <- subset(listSeurat[[j]], scDblFinder.score < 0.5)

GSE265972.merged <- merge(x = listSeurat[[1]], y = list(listSeurat[[2]], listSeurat[[3]], listSeurat[[4]], listSeurat[[5]], 
                                                        listSeurat[[6]], listSeurat[[7]], listSeurat[[8]], listSeurat[[9]]),
                          project = "GSE265972_merged", merge.data = TRUE)
DefaultAssay(GSE265972.merged) <- 'RNA'


#### GSE241132 ####

GSE241132_data_dir <- file.path(getwd(),"GSE241132")
GSE241132_data_dirs <- list.dirs(GSE241132_data_dir, full.names=TRUE, recursive=FALSE)
GSE241132_data_samples <- list.dirs(GSE241132_data_dir, full.names=FALSE, recursive=FALSE)
listSeurat <- list()
for (i in 1:length(GSE241132_data_dirs)) {
  sample <- Read10X_GEO(data_dir = GSE241132_data_dirs[i], gene.column = 2)
  seuratObject = toString(GSE241132_data_samples[i])
  seuratObject <- CreateSeuratObject(counts = sample)
  seuratObject@meta.data$sampleName <- GSE241132_data_samples[i]
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^MT-", col.name = "percent.mt")
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^RPL", col.name = "percent.rpl")
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^RPS", col.name = "percent.rps")
  listSeurat[i] <- seuratObject
}

j = 1 # cycle from j = 1 to 12
listSeurat[[j]]
VlnPlot(listSeurat[[j]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rpl", "percent.rps"), ncol = 5)
FeatureScatter(listSeurat[[j]], feature1 = "percent.mt", feature2 = "nCount_RNA")
FeatureScatter(listSeurat[[j]], feature1 = "percent.rps", feature2 = "nFeature_RNA")
FeatureScatter(listSeurat[[j]], feature1 = "percent.rpl", feature2 = "nCount_RNA")
VlnPlot(listSeurat[[j]], features = c("CSF3R","FCGR3B","CMTM2","PROK2","CXCR2"))
listSeurat[[j]] <- subset(listSeurat[[j]], subset = nFeature_RNA > 50 & nCount_RNA < 75000 & percent.mt < 20 & percent.rps < 20 & percent.rpl < 30)
DefaultAssay(listSeurat[[j]]) <- "RNA"
sce <- scDblFinder(GetAssayData(listSeurat[[j]], assay="RNA", slot="counts"), samples=Idents(listSeurat[[j]]))
listSeurat[[j]]$scDblFinder.score <- sce$scDblFinder.score
VlnPlot(listSeurat[[j]], features = "scDblFinder.score", raster=FALSE, pt.size=0.5)
listSeurat[[j]] <- subset(listSeurat[[j]], scDblFinder.score < 0.5)

GSE241132.merged <- merge(x = listSeurat[[1]], y = list(listSeurat[[2]], listSeurat[[3]], listSeurat[[4]], listSeurat[[5]], 
                                                        listSeurat[[6]], listSeurat[[7]], listSeurat[[8]], listSeurat[[9]],
                                                        listSeurat[[10]], listSeurat[[11]], listSeurat[[12]]), 
                          project = "GSE241132_merged", merge.data = TRUE)
DefaultAssay(GSE241132.merged) <- 'RNA'


#### GSE231643 ####

GSE231643_data_dir <- file.path(getwd(),"GSE231643")
GSE231643_data_dirs <- list.dirs(GSE231643_data_dir, full.names=TRUE, recursive=FALSE)
GSE231643_data_samples <- list.dirs(GSE231643_data_dir, full.names=FALSE, recursive=FALSE)
listSeurat <- list()
for (i in 1:length(GSE231643_data_dirs)) {
  sample <- Read10X_GEO(data_dir = GSE231643_data_dirs[i], gene.column = 2)
  seuratObject = toString(GSE231643_data_samples[i])
  seuratObject <- CreateSeuratObject(counts = sample)
  seuratObject@meta.data$sampleName <- GSE231643_data_samples[i]
  seuratObject@meta.data$sampleType <- substr(GSE231643_data_samples[i], start = 1, stop = 2)
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^MT-", col.name = "percent.mt")
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^RPL", col.name = "percent.rpl")
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^RPS", col.name = "percent.rps")
  listSeurat[i] <- seuratObject
}

j = 1 # cycle from j = 1 to 8
listSeurat[[j]]
VlnPlot(listSeurat[[j]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rpl", "percent.rps"), ncol = 5)
FeatureScatter(listSeurat[[j]], feature1 = "percent.mt", feature2 = "nCount_RNA")
FeatureScatter(listSeurat[[j]], feature1 = "percent.rps", feature2 = "nFeature_RNA")
FeatureScatter(listSeurat[[j]], feature1 = "percent.rpl", feature2 = "nCount_RNA")
listSeurat[[j]] <- subset(listSeurat[[j]], subset = nFeature_RNA > 100 & nCount_RNA < 50000 & percent.mt < 15 & percent.rps < 20 & percent.rpl < 30)
DefaultAssay(listSeurat[[j]]) <- "RNA"
sce <- scDblFinder(GetAssayData(listSeurat[[j]], assay="RNA", slot="counts"), samples=Idents(listSeurat[[j]]))
listSeurat[[j]]$scDblFinder.score <- sce$scDblFinder.score
VlnPlot(listSeurat[[j]], features = "scDblFinder.score", raster=FALSE, pt.size=0.5)
listSeurat[[j]] <- subset(listSeurat[[j]], scDblFinder.score < 0.5)

GSE231643.merged <- merge(x = listSeurat[[1]], y = list(listSeurat[[2]], listSeurat[[3]], listSeurat[[4]], listSeurat[[5]], 
                                                        listSeurat[[6]], listSeurat[[7]], listSeurat[[8]]),
                          project = "GSE231643_merged", merge.data = TRUE)
DefaultAssay(GSE231643.merged) <- 'RNA'


#### GSE268834 ####

GSE268834_data_dir <- file.path(getwd(),"GSE268834_RAW")
GSE268834_data_dirs <- list.dirs(GSE268834_data_dir, full.names=TRUE, recursive=FALSE)
GSE268834_data_samples <- list.dirs(GSE268834_data_dir, full.names=FALSE, recursive=FALSE)
listSeurat <- list()
for (i in 1:length(GSE268834_data_dirs)) {
  sample <- Read10X_GEO(data_dir = GSE268834_data_dirs[i], gene.column = 2)
  seuratObject = toString(GSE268834_data_samples[i])
  seuratObject <- CreateSeuratObject(counts = sample)
  seuratObject@meta.data$sampleName <- GSE268834_data_samples[i]
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^MT-", col.name = "percent.mt")
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^RPL", col.name = "percent.rpl")
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^RPS", col.name = "percent.rps")
  listSeurat[i] <- seuratObject
}

j = 1 # cycle from j = 1 to 8
listSeurat[[j]]
VlnPlot(listSeurat[[j]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rpl", "percent.rps"), ncol = 5)
FeatureScatter(listSeurat[[j]], feature1 = "percent.mt", feature2 = "nFeature_RNA")
FeatureScatter(listSeurat[[j]], feature1 = "percent.rps", feature2 = "nFeature_RNA")
FeatureScatter(listSeurat[[j]], feature1 = "percent.rpl", feature2 = "nCount_RNA")
listSeurat[[j]] <- subset(listSeurat[[j]], subset = nFeature_RNA > 50 & nCount_RNA < 75000 & percent.mt < 20 & percent.rps < 20 & percent.rpl < 25)
DefaultAssay(listSeurat[[j]]) <- "RNA"
sce <- scDblFinder(GetAssayData(listSeurat[[j]], assay="RNA", slot="counts"), samples=Idents(listSeurat[[j]]))
listSeurat[[j]]$scDblFinder.score <- sce$scDblFinder.score
VlnPlot(listSeurat[[j]], features = "scDblFinder.score", raster=FALSE, pt.size=0.5)
listSeurat[[j]] <- subset(listSeurat[[j]], scDblFinder.score < 0.5)

GSE268834.merged <- merge(x = listSeurat[[1]], y = list(listSeurat[[2]], listSeurat[[3]], listSeurat[[4]], listSeurat[[5]], 
                                                        listSeurat[[6]], listSeurat[[7]], listSeurat[[8]]),
                          project = "GSE268834_merged", merge.data = TRUE)
DefaultAssay(GSE268834.merged) <- 'RNA'


#### GSE248247 ####

GSE248247_data_dir <- file.path(getwd(),"Human diabetic healing/GSE248247_RAW")
GSE248247_data_dirs <- list.dirs(GSE248247_data_dir, full.names=TRUE, recursive=FALSE)
GSE248247_data_samples <- list.dirs(GSE248247_data_dir, full.names=FALSE, recursive=FALSE)
listSeurat <- list()
for (i in 1:length(GSE248247_data_dirs)) {
  sample <- Read10X_GEO(data_dir = GSE248247_data_dirs[i], gene.column = 2)
  seuratObject = toString(GSE248247_data_samples[i])
  seuratObject <- CreateSeuratObject(counts = sample)
  seuratObject@meta.data$sampleName <- GSE248247_data_samples[i]
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^MT-", col.name = "percent.mt")
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^RPL", col.name = "percent.rpl")
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^RPS", col.name = "percent.rps")
  listSeurat[i] <- seuratObject
}
j = 1 # cycle from j = 1 to 2
listSeurat[[j]]
VlnPlot(listSeurat[[j]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rpl", "percent.rps"), ncol = 5)
FeatureScatter(listSeurat[[j]], feature1 = "percent.mt", feature2 = "nCount_RNA")
FeatureScatter(listSeurat[[j]], feature1 = "percent.rps", feature2 = "nFeature_RNA")
FeatureScatter(listSeurat[[j]], feature1 = "percent.rpl", feature2 = "nCount_RNA")
listSeurat[[j]] <- subset(listSeurat[[j]], subset = nFeature_RNA > 50 & nCount_RNA < 75000 & percent.mt < 25 & percent.rps < 20 & percent.rpl < 30)
DefaultAssay(listSeurat[[j]]) <- "RNA"
sce <- scDblFinder(GetAssayData(listSeurat[[j]], assay="RNA", slot="counts"), samples=Idents(listSeurat[[j]]))
listSeurat[[j]]$scDblFinder.score <- sce$scDblFinder.score
VlnPlot(listSeurat[[j]], features = "scDblFinder.score", raster=FALSE, pt.size=0.5)
listSeurat[[j]] <- subset(listSeurat[[j]], scDblFinder.score < 0.5)

GSE248247.merged <- merge(x = listSeurat[[1]], y = list(listSeurat[[2]]), project = "GSE248247_merged", merge.data = TRUE)
DefaultAssay(GSE248247.merged) <- 'RNA'


#### INTEGRATIVE ANALYSIS OF HUMAN CHRONIC WOUNDS ####

GSE231643.merged[["orig.dataset"]] <- "GSE231643"
GSE223964.merged[["orig.dataset"]] <- "GSE223964"
GSE265972.merged[["orig.dataset"]] <- "GSE265972"
GSE241132.merged[["orig.dataset"]] <- "GSE241132"
GSE248247.merged[["orig.dataset"]] <- "GSE248247"
GSE268834.merged[["orig.dataset"]] <- "GSE268834"

HumanChronicWounds <- merge(x = GSE231643.merged, y = list(GSE223964.merged, GSE265972.merged, GSE241132.merged, GSE248247.merged, GSE268834.merged), 
                            project = "HumanChronicWounds", merge.data = TRUE)
DefaultAssay(HumanChronicWounds)<-'RNA'

HumanChronicWounds <- NormalizeData(HumanChronicWounds)
HumanChronicWounds <- FindVariableFeatures(HumanChronicWounds)

HumanChronicWounds[["RNA"]] <- JoinLayers(HumanChronicWounds[['RNA']])
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
HumanChronicWounds <- CellCycleScoring(HumanChronicWounds, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

HumanChronicWounds <- ScaleData(HumanChronicWounds, vars.to.regress = c("S.Score", "G2M.Score"))
HumanChronicWounds <- ScaleData(HumanChronicWounds)
HumanChronicWounds <- RunPCA(HumanChronicWounds)

HumanChronicWounds[['RNA']] <- split(x = HumanChronicWounds[['RNA']], f = HumanChronicWounds$orig.dataset)

HumanChronicWounds <- IntegrateLayers(object = HumanChronicWounds, 
                                      method = RPCAIntegration, 
                                      orig.reduction = "pca", 
                                      new.reduction = "integrated.rpca",
                                      verbose = FALSE)
HumanChronicWounds <- JoinLayers(HumanChronicWounds)
DefaultAssay(HumanChronicWounds) <- "RNA"

# Seurat workflow without integration
HumanChronicWounds <- FindNeighbors(HumanChronicWounds, reduction = "pca", dims = 1:15)
HumanChronicWounds <- FindClusters(HumanChronicWounds, resolution = .3, cluster.name = "unintegrated_clusters")
HumanChronicWounds <- RunUMAP(HumanChronicWounds, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")

# Figure S11A
DimPlot(
  HumanChronicWounds, 
  reduction = "umap.unintegrated", 
  group.by = "orig.dataset",
  shuffle = TRUE,
  raster = FALSE
)

# Figure S11B
DimPlot(
  HumanChronicWounds, 
  reduction = "umap.unintegrated", 
  group.by = "sampleName",
  shuffle = TRUE,
  raster = FALSE
)

# Seurat workflow with RPCA integration
HumanChronicWounds <- FindNeighbors(HumanChronicWounds, reduction = "integrated.rpca", dims = 1:15)
HumanChronicWounds <- FindClusters(HumanChronicWounds, resolution = .3, cluster.name = "rpca_clusters")
HumanChronicWounds <- RunUMAP(HumanChronicWounds, dims = 1:15, reduction = "integrated.rpca", reduction.name = "umap.rpca")

# Figure S11C
DimPlot(
  HumanChronicWounds,
  reduction = "umap.rpca",
  group.by = "wound_type",
  shuffle = TRUE,
  raster = FALSE
)

# Figure S11D
DimPlot(
  HumanChronicWounds, 
  reduction = "umap.rpca", 
  group.by = "sampleName",
  shuffle = TRUE,
  raster = FALSE
)

Idents(HumanChronicWounds) <- "rpca_clusters"
levels(HumanChronicWounds)
HumanChronicWounds[["int_cell_types"]] <- Idents(HumanChronicWounds)
Idents(HumanChronicWounds) <- "int_cell_types"
levels(HumanChronicWounds)
HumanChronicWounds <- RenameIdents(HumanChronicWounds, 
                                   "0"="Keratinocyte", 
                                   "1"="Fibroblast", 
                                   "2"="Endothelial cell", 
                                   "3"="Macrophage", 
                                   "4"="T cell", 
                                   "5"="Keratinocyte", 
                                   "6"="Smooth muscle cell", 
                                   "7"="Fibroblast",
                                   "8"="Fibroblast",
                                   "9"="Keratinocyte", 
                                   "10"="B cell",
                                   "11"="Keratinocyte",
                                   "12"="Plasma cell",
                                   "13"="Macrophage",
                                   "14"="Neutrophil",
                                   "15"="Mast cell",
                                   "16"="Pericyte",
                                   "17"="Melanocyte",
                                   "18"="Keratinocyte"
)
levels(HumanChronicWounds)
HumanChronicWounds[["int_cell_types"]] <- Idents(HumanChronicWounds)

# Subset to chronic wounds and healthy skin only
Idents(HumanChronicWounds) <- "wound_type"
HumanChronicWounds <- subset(HumanChronicWounds, idents = c("Diabetic foot ulcer", "Non-diabetic foot ulcer", "Healthy skin"))

## Figure 7C
DimPlot(
  HumanChronicWounds,
  reduction = "umap.rpca",
  group.by = "int_cell_types",
  shuffle = TRUE,
  label = TRUE,
  raster = FALSE
)

## Figure 7D
DotPlot(HumanChronicWounds, group.by = "int_cell_types", 
        features = c("KRT14","PDGFRA","PECAM1","MRC1","CD3D","ACTA2","CPA3","CD38","CSF3R","MS4A1","MLANA")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "darkblue", mid = "white", high = "darkred")

## Figure 7E
DimPlot(
  HumanChronicWounds,
  reduction = "umap.rpca",
  group.by = "wound_type",
  shuffle = TRUE,
  label = FALSE,
  raster = FALSE
)

## Figure 7F
DimPlot(
  HumanChronicWounds,
  reduction = "umap.rpca",
  group.by = "sampleName",
  shuffle = TRUE,
  label = FALSE,
  raster = FALSE
)

## Figure 7G
table(HumanChronicWounds$int_cell_types, HumanChronicWounds$wound_type)

# Figure S11E
table(HumanChronicWounds$int_cell_types, HumanChronicWounds$sampleName)

# # Save the integrated human wound  cell dataset as RDS file
saveRDS(HumanChronicWounds, "HumanChronicWounds_integrated_final.rds")
# Read in the integrated human wound  cell dataset from RDS file
HumanChronicWounds <- readRDS("HumanChronicWounds_integrated_final.rds")


#### Human wound cell subtypes ####

# Subset fibroblasts, macrophages and neutrophils
Idents(HumanChronicWounds) <- "int_cell_types"
HumanChronicWounds_fibroblast <- subset(HumanChronicWounds, idents = "Fibroblast")
HumanChronicWounds_macrophage <- subset(HumanChronicWounds, idents = "Macrophage")
HumanChronicWounds_neutrophil <- subset(HumanChronicWounds, idents = "Neutrophil")

DefaultAssay(HumanChronicWounds_fibroblast) <- "RNA"
HumanChronicWounds_fibroblast <- RunPCA(HumanChronicWounds_fibroblast)
DefaultAssay(HumanChronicWounds_macrophage) <- "RNA"
HumanChronicWounds_macrophage <- RunPCA(HumanChronicWounds_macrophage)
DefaultAssay(HumanChronicWounds_neutrophil) <- "RNA"
HumanChronicWounds_neutrophil <- RunPCA(HumanChronicWounds_neutrophil)

# Figure S12A
ElbowPlot(HumanChronicWounds_fibroblast)
ElbowPlot(HumanChronicWounds_macrophage)
ElbowPlot(HumanChronicWounds_neutrophil)

# Fibroblasts
HumanChronicWounds_fibroblast <- FindNeighbors(HumanChronicWounds_fibroblast, reduction = "integrated.rpca", dims = 1:6)
HumanChronicWounds_fibroblast <- FindClusters(HumanChronicWounds_fibroblast, resolution = .15)
HumanChronicWounds_fibroblast <- RunUMAP(HumanChronicWounds_fibroblast, dims = 1:6, reduction = "integrated.rpca")

# Figure S12B
DimPlot(
  HumanChronicWounds_fibroblast,
  reduction = "umap",
  group.by = c("wound_type", "seurat_clusters", "cell_subtypes"),
  combine = TRUE,
  label = FALSE,
  raster = FALSE
)

# Macrophages
HumanChronicWounds_macrophage <- FindNeighbors(HumanChronicWounds_macrophage, reduction = "integrated.rpca", dims = 1:10)
HumanChronicWounds_macrophage <- FindClusters(HumanChronicWounds_macrophage, resolution = .2)
HumanChronicWounds_macrophage <- RunUMAP(HumanChronicWounds_macrophage, dims = 1:10, reduction = "integrated.rpca")

# Figure S12C
DimPlot(
  HumanChronicWounds_macrophage,
  reduction = "umap",
  group.by = c("wound_type", "seurat_clusters", "cell_subtypes"),
  combine = TRUE,
  label = FALSE,
  raster = FALSE
)

# Neutrophils
HumanChronicWounds_neutrophil <- FindNeighbors(HumanChronicWounds_neutrophil, reduction = "integrated.rpca", dims = 1:8)
HumanChronicWounds_neutrophil <- FindClusters(HumanChronicWounds_neutrophil, resolution = .1)
HumanChronicWounds_neutrophil <- RunUMAP(HumanChronicWounds_neutrophil, dims = 1:8, reduction = "integrated.rpca")

# Figure S12D
DimPlot(
  HumanChronicWounds_neutrophil,
  reduction = "umap",
  group.by = c("wound_type", "seurat_clusters", "cell_subtypes"),
  combine = TRUE,
  label = FALSE,
  raster = FALSE
)

# Naming of all fibroblast, macrophage and neutrophil subtypes prior to merging

# Naming of fibroblast subtypes
Idents(HumanChronicWounds_fibroblast) <- "seurat_clusters"
HumanChronicWounds_fibroblast[["final_subtypes"]] <- Idents(HumanChronicWounds_fibroblast)
Idents(HumanChronicWounds_fibroblast) <- "final_subtypes"
HumanChronicWounds_fibroblast <- RenameIdents(HumanChronicWounds_fibroblast, 
                                      "0"="Fb_LRRC15+_TNC+", 
                                      "1"="Fb_CILP+_MGP+", 
                                      "2"="Fb_PRG4+_DPP4+",
                                      "3"="Fb_LRRC15+_TNC+", 
                                      "4"="Fb_LRRC15+_TNC+"
)
levels(HumanChronicWounds_fibroblast)
HumanChronicWounds_fibroblast[["final_subtypes"]] <- Idents(HumanChronicWounds_fibroblast)

# Naming of macrophage subtypes
Idents(HumanChronicWounds_macrophage) <- "seurat_clusters"
HumanChronicWounds_macrophage[["final_subtypes"]] <- Idents(HumanChronicWounds_macrophage)
Idents(HumanChronicWounds_macrophage) <- "final_subtypes"
HumanChronicWounds_macrophage <- RenameIdents(HumanChronicWounds_macrophage, 
                                      "0"="Mo_CD163+_MRC1+", 
                                      "1"="Mo_CD163+_MRC1+", 
                                      "2"="Mo_VCAN+_IL1B+", 
                                      "3"="Mo_CD86+_CCR7+",
                                      "4"="Mo_CD86+_CCR7+"
)
levels(HumanChronicWounds_macrophage)
HumanChronicWounds_macrophage[["final_subtypes"]] <- Idents(HumanChronicWounds_macrophage)

# Naming of neutrophil subtypes
Idents(HumanChronicWounds_neutrophil) <- "seurat_clusters"
HumanChronicWounds_neutrophil[["final_subtypes"]] <- Idents(HumanChronicWounds_neutrophil)
Idents(HumanChronicWounds_neutrophil) <- "final_subtypes"
HumanChronicWounds_neutrophil <- RenameIdents(HumanChronicWounds_neutrophil, 
                                         "0"="Neu_TNFAIP3+_SOD2+", 
                                         "1"="Neu_CSF3R+_FOS+"
)
levels(HumanChronicWounds_neutrophil)
HumanChronicWounds_neutrophil[["final_subtypes"]] <- Idents(HumanChronicWounds_neutrophil)

# Merging of cell subtype datasets
HumanChronicWounds_subtypes <- merge(HumanChronicWounds_neutrophil, y = c(HumanChronicWounds_macrophage, HumanChronicWounds_fibroblast))

# Read in human cell subtypes dataset from RDS file
HumanChronicWounds_subtypes <- readRDS("HumanChronicWounds_subtypes_final.rds")
# Save the human cell subtypes dataset as RDS file
saveRDS(HumanChronicWounds_subtypes, "HumanChronicWounds_subtypes_final.rds")

# List the annotated cell subtypes
Idents(HumanChronicWounds_subtypes) <- "final_subtypes"
levels(HumanChronicWounds_subtypes$final_subtypes)

# List the annotated wound types
levels(HumanChronicWounds_subtypes$wound_type)

## Figure 7I
DotPlot(HumanChronicWounds_subtypes, group.by = "final_subtypes",
        features = c("PTPRC","ITGAM","TNFAIP3","SOD2","CSF3R","FOS","CD163","MRC1","VCAN","IL1B","CD86","CCR7","PDGFRA","DPP4","PRG4","LRRC15","TNC","CILP","MGP")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

## Figure 7J
table(HumanChronicWounds_subtypes$final_subtypes, HumanChronicWounds_subtypes$wound_type)

# Figure S11F
table(HumanChronicWounds_subtypes$final_subtypes, HumanChronicWounds_subtypes$sampleName)

## Figure 7K
DotPlot(HumanChronicWounds_subtypes, group.by = "final_subtypes",
        features = c("mm_neutDB_score1","mm_neutWT_score1","mm_F480_score1","mm_Ly6c_score1","mm_Cd11c_score1","mm_Cd26_score1","mm_Inhba_score1","mm_Cilp_score1")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")


#### Figure 8 ####

# Subset to chronic wounds only
Idents(HumanChronicWounds_subtypes) <- "wound_type"
HumanChronicWounds_subtypes <- subset(HumanChronicWounds_subtypes, idents = c("Diabetic foot ulcer","Non-diabetic foot ulcer"))

## Figure 8A
Idents(HumanChronicWounds_subtypes) <- "final_subtypes"
CD44_pathway_hs <- c("CD44","SPP1","CXCL2","CXCR2","ICAM1","VCAM1","FN1")
DotPlot(HumanChronicWounds_subtypes, group.by = "final_subtypes", split.by = "wound_type", cols = "Set3", 
        features = CD44_pathway_hs) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")

# Load in gene sets for module scoring
human_M1_features <- c("IL1B", "TNF", "IL6", "CXCL9", "CXCL10", "CXCL11", "CD86", "CD80", 
                       "CD38", "MARCO", "S100A8", "S100A9", "ACOD1", "GBP2", "IRF5", 
                       "STAT1", "SOCS1", "IFIT1", "CMPK2") 
human_M2a_features <- c("MRC1", "TGM2", "CCL17", "CCL22", "CCL24", "CD200R1", "CLEC10A", 
                        "PPARG", "STAT6", "DAB2", "PDCD1LG2", "ALOX15", "CD36", "FABP4", 
                        "MAFB", "MS4A4A", "MYC")
human_M2c_features <- c("CD163", "MERTK", "IL10", "TGFB1", "STAT3", "CCL18", "CCL16", 
                        "PGF", "MMP9", "MMP7", "MMP8", "TREM2", "C1QA", "C1QB", 
                        "C1QC", "GAS6", "SLCO2B1", "LYZ", "MPO")
human_mac_lists <- list(M1 = human_M1_features, 
                        M2a = human_M2a_features, 
                        M2c = human_M2c_features)
# Perform module scoring
HumanChronicWounds_subtypes <- AddModuleScore(object = HumanChronicWounds_subtypes, features = human_mac_lists, name = "MacPolarization")


## Figure 8B
Idents(HumanChronicWounds_subtypes) <- "final_subtypes"
DotPlot(HumanChronicWounds_subtypes, ident = c("Mo_CD163+_MRC1+", "Mo_VCAN+_IL1B+", "Mo_CD86+_CCR7+"), group.by = "final_subtypes", split.by = "wound_type", cols = "Set3",     
        features = c("MacPolarization1","MacPolarization2","MacPolarization3","CD44","SPP1","ICAM1","S100A9","FPR1","SOD2","NFKB1","IL1B","TNF","MERTK","CD36")) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "white", high = "darkred")


#### CellChat analysis of human cell subtypes ####

# Subset according to diabetic vs non-diabetic wound type
Idents(HumanChronicWounds_subtypes) <- "wound_type"
HumanChronicWounds_subtypes_D <- subset(HumanChronicWounds_subtypes, idents = "Diabetic foot ulcer")
HumanChronicWounds_subtypes_ND <- subset(HumanChronicWounds_subtypes, idents = "Non-diabetic foot ulcer")
# Set cell subtypes as the default identity for analysis of cell-cell interactions
Idents(HumanChronicWounds_subtypes_D) <- "final_subtypes"
Idents(HumanChronicWounds_subtypes_ND) <- "final_subtypes"
# Create CellChat Objects
cellchat_HumanChronicWounds_D <- createCellChat(HumanChronicWounds_subtypes_D, group.by = "final_subtypes")
cellchat_HumanChronicWounds_ND <- createCellChat(HumanChronicWounds_subtypes_ND, group.by = "final_subtypes")
# Set cell-cell interaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
# set the used database in the object
cellchat_HumanChronicWounds_D@DB <- CellChatDB.use
cellchat_HumanChronicWounds_ND@DB <- CellChatDB.use
# Process the datasets for CellChat
cellchat_HumanChronicWounds_D <- subsetData(cellchat_HumanChronicWounds_D) 
cellchat_HumanChronicWounds_ND <- subsetData(cellchat_HumanChronicWounds_ND) 
# Run the standard CellChat workflow
future::plan("multisession", workers = 4)
cellchat_HumanChronicWounds_D <- identifyOverExpressedGenes(cellchat_HumanChronicWounds_D)
cellchat_HumanChronicWounds_ND <- identifyOverExpressedGenes(cellchat_HumanChronicWounds_ND)
cellchat_HumanChronicWounds_D <- identifyOverExpressedInteractions(cellchat_HumanChronicWounds_D)
cellchat_HumanChronicWounds_ND <- identifyOverExpressedInteractions(cellchat_HumanChronicWounds_ND)
cellchat_HumanChronicWounds_D <- computeCommunProb(cellchat_HumanChronicWounds_D, type = "triMean", population.size = TRUE) 
cellchat_HumanChronicWounds_ND <- computeCommunProb(cellchat_HumanChronicWounds_ND, type = "triMean", population.size = TRUE) 
cellchat_HumanChronicWounds_D <- filterCommunication(cellchat_HumanChronicWounds_D, min.cells = 10)
cellchat_HumanChronicWounds_ND <- filterCommunication(cellchat_HumanChronicWounds_ND, min.cells = 10)
cellchat_HumanChronicWounds_D <- computeCommunProbPathway(cellchat_HumanChronicWounds_D)
cellchat_HumanChronicWounds_ND <- computeCommunProbPathway(cellchat_HumanChronicWounds_ND)
cellchat_HumanChronicWounds_D <- aggregateNet(cellchat_HumanChronicWounds_D)
cellchat_HumanChronicWounds_ND <- aggregateNet(cellchat_HumanChronicWounds_ND)
cellchat_HumanChronicWounds_D <- netAnalysis_computeCentrality(cellchat_HumanChronicWounds_D, slot.name = "netP")
cellchat_HumanChronicWounds_ND <- netAnalysis_computeCentrality(cellchat_HumanChronicWounds_ND, slot.name = "netP")

# Save CellChat objects as RDS files
saveRDS(cellchat_HumanChronicWounds_ND, "CellChat_human_NonDiabetic.rds")
saveRDS(cellchat_HumanChronicWounds_D, "CellChat_human_Diabetic.rds")
# Load CellChat objects from RDS files
cellchat_HumanChronicWounds_ND <- readRDS("CellChat_human_NonDiabetic.rds")
cellchat_HumanChronicWounds_D <- readRDS("CellChat_human_Diabetic.rds")

# Figure S13A
netAnalysis_signalingRole_scatter(cellchat_HumanChronicWounds_ND)
netAnalysis_signalingRole_scatter(cellchat_HumanChronicWounds_D)

# Figure S13B left
groupSize <- as.numeric(table(cellchat_HumanChronicWounds_ND@idents))
netVisual_circle(cellchat_HumanChronicWounds_ND@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_HumanChronicWounds_ND@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Figure S13B right
groupSize <- as.numeric(table(cellchat_HumanChronicWounds_D@idents))
netVisual_circle(cellchat_HumanChronicWounds_D@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_HumanChronicWounds_D@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Differential CellChat analysis
object.list_D_v_ND <- list(ND = cellchat_HumanChronicWounds_ND, D = cellchat_HumanChronicWounds_D)
cellchat_D_v_ND <- mergeCellChat(object.list_D_v_ND, add.names = names(object.list_D_v_ND))

# Figure S13C
gg1 <- compareInteractions(cellchat_D_v_ND, show.legend = F)
gg2 <- compareInteractions(cellchat_D_v_ND, show.legend = F, measure = "weight")
gg1 + gg2

# Figure S13D
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat_D_v_ND, weight.scale = T)
netVisual_diffInteraction(cellchat_D_v_ND, weight.scale = T, measure = "weight")

## Figure 8C left
netVisual_heatmap(cellchat_D_v_ND)

## Figure 8C right
netVisual_heatmap(cellchat_D_v_ND, measure = "weight")

## Figure 8D
netAnalysis_signalingChanges_scatter(cellchat_D_v_ND, idents.use = "Mo_CD163+_MRC1+", top.label = 1)

## Figure 8E
netAnalysis_signalingChanges_scatter(cellchat_D_v_ND, idents.use = "Mo_VCAN+_IL1B+", top.label = 1)

# Figure S13E
netAnalysis_signalingChanges_scatter(cellchat_D_v_ND, idents.use = "Neu_TNFAIP3+_SOD2+", top.label = 1)
netAnalysis_signalingChanges_scatter(cellchat_D_v_ND, idents.use = "Neu_CSF3R+_FOS+", top.label = 1)
netAnalysis_signalingChanges_scatter(cellchat_D_v_ND, idents.use = "Mo_CD86+_CCR7+", top.label = 1)
netAnalysis_signalingChanges_scatter(cellchat_D_v_ND, idents.use = "Fb_PRG4+_DPP4+", top.label = 1)
netAnalysis_signalingChanges_scatter(cellchat_D_v_ND, idents.use = "Fb_LRRC15+_TNC+", top.label = 1)
netAnalysis_signalingChanges_scatter(cellchat_D_v_ND, idents.use = "Fb_CILP+_MGP+", top.label = 1)

# Code to perform differential expression analysis of ligands and receptors in D vs ND
pos.dataset = "Diabetic foot ulcer"
cellchat_D_v_ND@var.features$features
cellchat_D_v_ND@var.features[[features.name]]
features.name = paste0(pos.dataset, ".merged")
cellchat_D_v_ND <- identifyOverExpressedGenes(cellchat_D_v_ND, group.dataset = "wound_type", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 0.05) 
net <- netMappingDEG(cellchat_D_v_ND, features.name = features.name, thresh = 0.05)
write.csv(net, file = file.path(getwd(), "CellChat_huamn_DEG_LRpairs.txt"))
net.up <- subsetCommunication(cellchat_D_v_ND, net = net, datasets = "D", ligand.logFC = 0.20, receptor.logFC = 0.20)
net.down <- subsetCommunication(cellchat_D_v_ND, net = net, datasets = "D", ligand.logFC = -0.20, receptor.logFC = -0.20)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat_D_v_ND)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_D_v_ND)
pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]
sort.pairLR.use.up <- arrange(pairLR.use.up, interaction_name)
sort.pairLR.use.down <- arrange(pairLR.use.down, interaction_name)

## Figure 8F
netVisual_bubble(cellchat_D_v_ND, pairLR.use = pairLR.use.up, sources.use = "Neu_TNFAIP3+_SOD2+", comparison = 2,  angle.x = 90, remove.isolate = T, title.name = paste0("Up-regulated signaling in ", names(object.list_D_v_ND)[2]))

## Figure 8G
netVisual_bubble(cellchat_D_v_ND, pairLR.use = pairLR.use.up, sources.use = "Mo_CD163+_MRC1+", comparison = 2,  angle.x = 90, remove.isolate = T, title.name = paste0("Up-regulated signaling in ", names(object.list_D_v_ND)[2]))

## Figure 8H
netVisual_bubble(cellchat_D_v_ND, pairLR.use = pairLR.use.up, sources.use = "Mo_VCAN+_IL1B+", comparison = 2,  angle.x = 90, remove.isolate = T, title.name = paste0("Up-regulated signaling in ", names(object.list_D_v_ND)[2]))
