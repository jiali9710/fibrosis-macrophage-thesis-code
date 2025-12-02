# 04_Subset_MP.R
# Description: Subsetting mononuclear phagocytes (MPs) from liver, lung, and kidney immune & blood datasets for integration and clustering

# ===============================
# 0. Load Required Packages
# ===============================
library(Seurat)
library(future)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(Rcpp)
library(harmony)
library(clustree)
library(cowplot)
library(data.table)
library(stringr)

ulimit::memory_limit(120000)

# ===============================
# 1. Set Output Paths
# ===============================
dir.create("output/Subset_MP_Liver/", recursive = TRUE)
dir.create("output/Subset_MP_Lung/", recursive = TRUE)
dir.create("output/Subset_MP_Kidney/", recursive = TRUE)

# ===============================
# 2. LIVER: Subset MPs and Integrate
# ===============================
# Load annotated liver object
liver_all <- readRDS("output/Harmony/liver_seurat_annotated.rds")
liver_mps <- subset(liver_all, subset = lineage == "mononuclearphagocyte")

# Subset raw merged object using barcodes
write.csv(rownames(liver_mps), "output/Subset_MP_Liver/liver_MP_cellbarcodes.csv", row.names = FALSE)
barcodes <- read.csv("output/Subset_MP_Liver/liver_MP_cellbarcodes.csv", stringsAsFactors = FALSE)[,1]

liver_merged <- readRDS("output/Harmony/liver_immuneblood_afterharmony.rds")
liver_mps <- subset(liver_merged, cells = barcodes)

# Normalize and integrate
DefaultAssay(liver_mps) <- "RNA"
liver_mps <- NormalizeData(liver_mps)
liver_mps <- FindVariableFeatures(liver_mps, selection.method = "vst", nfeatures = 2000)
liver_mps <- ScaleData(liver_mps)
liver_mps <- RunPCA(liver_mps, features = VariableFeatures(liver_mps))
liver_mps <- RunHarmony(liver_mps, group.by.vars = "sampleid", plot_convergence = TRUE, lambda = 2)
liver_mps <- RunUMAP(liver_mps, dims = 1:50, reduction = "harmony")
liver_mps <- FindNeighbors(liver_mps, dims = 1:50, reduction = "harmony")
liver_mps <- FindClusters(liver_mps, resolution = 0.7)

# Save
saveRDS(liver_mps, file = "output/Subset_MP_Liver/liver_MP_seurat_harmony.rds")

# ===============================
# 3. LUNG: Subset MPs and Integrate
# ===============================
lung_all <- readRDS("output/Harmony/lung_seurat_annotated.rds")
lung_mps <- subset(lung_all, subset = lineage == "mononuclearphagocyte")

write.csv(rownames(lung_mps), "output/Subset_MP_Lung/lung_MP_cellbarcodes.csv", row.names = FALSE)
barcodes <- read.csv("output/Subset_MP_Lung/lung_MP_cellbarcodes.csv", stringsAsFactors = FALSE)[,1]

lung_merged <- readRDS("output/Harmony/lung_immuneblood_afterharmony.rds")
lung_mps <- subset(lung_merged, cells = barcodes)

DefaultAssay(lung_mps) <- "RNA"
lung_mps <- NormalizeData(lung_mps)
lung_mps <- FindVariableFeatures(lung_mps, selection.method = "vst", nfeatures = 2000)
lung_mps <- ScaleData(lung_mps)
lung_mps <- RunPCA(lung_mps, features = VariableFeatures(lung_mps))
lung_mps <- RunHarmony(lung_mps, group.by.vars = "sampleid", plot_convergence = TRUE, lambda = 2)
lung_mps <- RunUMAP(lung_mps, dims = 1:50, reduction = "harmony")
lung_mps <- FindNeighbors(lung_mps, dims = 1:50, reduction = "harmony")
lung_mps <- FindClusters(lung_mps, resolution = 0.7)

saveRDS(lung_mps, file = "output/Subset_MP_Lung/lung_MP_seurat_harmony.rds")

# ===============================
# 4. KIDNEY: Subset MPs and Integrate
# ===============================
kidney_all <- readRDS("output/Harmony/kidney_seurat_annotated.rds")
kidney_mps <- subset(kidney_all, subset = lineage == "mononuclearphagocyte")

write.csv(rownames(kidney_mps), "output/Subset_MP_Kidney/kidney_MP_cellbarcodes.csv", row.names = FALSE)
barcodes <- read.csv("output/Subset_MP_Kidney/kidney_MP_cellbarcodes.csv", stringsAsFactors = FALSE)[,1]

kidney_merged <- readRDS("output/Harmony/kidney_immuneblood_afterharmony.rds")
kidney_mps <- subset(kidney_merged, cells = barcodes)

DefaultAssay(kidney_mps) <- "RNA"
kidney_mps <- NormalizeData(kidney_mps)
kidney_mps <- FindVariableFeatures(kidney_mps, selection.method = "vst", nfeatures = 2000)
kidney_mps <- ScaleData(kidney_mps)
kidney_mps <- RunPCA(kidney_mps, features = VariableFeatures(kidney_mps))
kidney_mps <- RunHarmony(kidney_mps, group.by.vars = "sampleid", plot_convergence = TRUE, lambda = 2)
kidney_mps <- RunUMAP(kidney_mps, dims = 1:50, reduction = "harmony")
kidney_mps <- FindNeighbors(kidney_mps, dims = 1:50, reduction = "harmony")
kidney_mps <- FindClusters(kidney_mps, resolution = 0.7)

saveRDS(kidney_mps, file = "output/Subset_MP_Kidney/kidney_MP_seurat_harmony.rds")