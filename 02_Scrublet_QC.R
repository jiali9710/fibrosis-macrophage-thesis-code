# 02_Scrublet_QC.R - Doublet Detection, QC, and Harmony Integration
# Author: Lily Jia
# Description:
# This script performs doublet detection using Scrublet (via Python),
# quality control filtering, and Harmony-based integration for immune
# and blood single-cell RNA-seq samples from liver, lung, and kidney.

# =============================
# 1. Load required packages
# =============================
library(Seurat)
library(reticulate)
library(dplyr)
library(patchwork)
library(cowplot)
library(harmony)
library(ggplot2)

use_python(Sys.which("python3"), required = TRUE)
scr <- import("scrublet")
plt <- import("matplotlib.pyplot")

# =============================
# 2. Define directories
# =============================
data_dir <- "../data/Soupx_corrected"
output_dir <- "../data/Scrublet_labeled"
harmony_dir <- "../data/harmony"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(harmony_dir, showWarnings = FALSE, recursive = TRUE)

# =============================
# 3. Example: Run Scrublet on one sample (ml_immune_ccl4)
# =============================
# Replace this with a loop or apply to all samples
input_file <- file.path(data_dir, "ml_immune_ccl4.rds")
seurat_obj <- readRDS(input_file)

# Run Scrublet on raw RNA matrix
doublet_obj <- scr$Scrublet(t(as.matrix(seurat_obj@assays$RNA@data)), 
                            expected_doublet_rate=0.1,
                            sim_doublet_ratio=5)
doublet_out <- doublet_obj$scrub_doublets()
doublet_obj$plot_histogram()
plt$savefig("hist_doublet_threshold.png")

# Add doublet label
seurat_obj$scrublet_label <- doublet_out[[2]]
DimPlot(seurat_obj, group.by="scrublet_label", reduction = "umap")
dev.print(png, "UMAP_doublet_labeled.png", width=5,height=5,res=300,units="in")

# Save labeled object
saveRDS(seurat_obj, file=file.path(output_dir, "ml_immune_ccl4_scrublet_labeled.rds"))

# =============================
# 4. Merge LIVER samples
# =============================
ml_immune_ccl4 <- readRDS(file.path(output_dir, "ml_immune_ccl4_scrublet_labeled.rds"))
ml_immune_uninjC <- readRDS(file.path(output_dir, "ml_immune_uninjC_scrublet_labeled.rds"))
ml_immune_uninjD <- readRDS(file.path(output_dir, "ml_immune_uninjD_scrublet_labeled.rds"))
ml_blood_ccl4 <- readRDS(file.path(output_dir, "ml_blood_ccl4_scrublet_labeled.rds"))
ml_blood_uninj <- readRDS(file.path(output_dir, "ml_blood_uninj_scrublet_labeled.rds"))
mk_blood_uninj <- readRDS(file.path(output_dir, "mk_blood_uninj_scrublet_labeled.rds"))

liver_seurat <- merge(ml_immune_ccl4, ml_immune_uninjC, add.cell.ids = c("ml_immune_ccl4","ml_immune_uninjC"))
liver_seurat <- merge(liver_seurat, ml_immune_uninjD, add.cell.ids = c("","ml_immune_uninjD"))
liver_seurat <- merge(liver_seurat, ml_blood_ccl4, add.cell.ids = c("","ml_blood_ccl4"))
liver_seurat <- merge(liver_seurat, ml_blood_uninj, add.cell.ids = c("","ml_blood_uninj"))
liver_seurat <- merge(liver_seurat, mk_blood_uninj, add.cell.ids = c("","mk_blood_uninj"))

# QC filtering
VlnPlot(liver_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), group.by="sampleid")
liver_seurat <- subset(liver_seurat, subset = nFeature_RNA > 400 & percent.mito < 15)

# Normalization and variable gene selection
DefaultAssay(liver_seurat) <- "RNA"
liver_seurat <- NormalizeData(liver_seurat)
liver_seurat <- FindVariableFeatures(liver_seurat)
liver_seurat <- ScaleData(liver_seurat)
liver_seurat <- RunPCA(liver_seurat)

# Harmony integration (LIVER)
pc.num <- 1:50
liver_seurat.harmony <- RunHarmony(liver_seurat, group.by.vars="orig.ident", plot_convergence=TRUE, lambda=2)
liver_harmony_embeddings <- Embeddings(liver_seurat.harmony, 'harmony')
liver_harmony_embeddings[1:5,1:5]
liver_seurat.harmony <- liver_seurat.harmony %>% FindNeighbors(dims=pc.num) %>% FindClusters(resolution=0.6)
liver_seurat.harmony <- liver_seurat.harmony %>% RunUMAP(dims=pc.num, reduction="harmony")

p1 <- DimPlot(liver_seurat.harmony, reduction="harmony", group.by="sampleid")
p2 <- VlnPlot(liver_seurat.harmony, features="harmony_1", group.by="sampleid")
plot_grid(p1, p2)

FeaturePlot(liver_seurat.harmony, c("Spp1","Trem2"))
DimPlot(liver_seurat.harmony, group.by="orig.ident") + plot_annotation(title="liver after integration (Harmony)")
DimPlot(liver_seurat.harmony, group.by="orig.ident", split.by="orig.ident") + plot_annotation(title="liver after integration (Harmony)")

saveRDS(liver_seurat.harmony, file=file.path(harmony_dir, "liver_immuneblood_afterharmony.rds"))

# =============================
# 5. Merge and integrate LUNG samples (full code retained)
# =============================
mlung_immune_bleo <- readRDS(file.path(output_dir, "mlung_immune_bleo_scrublet_labeled.rds"))
mlung_immune_uninj <- readRDS(file.path(output_dir, "mlung_immune_uninj_scrublet_labeled.rds"))
mlung_blood_bleo <- readRDS(file.path(output_dir, "mlung_blood_bleo_scrublet_labeled.rds"))
ml_blood_uninj <- readRDS(file.path(output_dir, "ml_blood_uninj_scrublet_labeled.rds"))
mk_blood_uninj <- readRDS(file.path(output_dir, "mk_blood_uninj_scrublet_labeled.rds"))

lung_seurat <- merge(mlung_immune_bleo, mlung_immune_uninj, add.cell.ids = c("mlung_immune_bleo","mlung_immune_uninj"))
lung_seurat <- merge(lung_seurat, mlung_blood_bleo, add.cell.ids = c("","mlung_blood_bleo"))
lung_seurat <- merge(lung_seurat, ml_blood_uninj, add.cell.ids = c("","ml_blood_uninj"))
lung_seurat <- merge(lung_seurat, mk_blood_uninj, add.cell.ids = c("","mk_blood_uninj"))

VlnPlot(lung_seurat, features=c("nFeature_RNA","nCount_RNA","percent.mito"), group.by="sampleid")
lung_seurat <- subset(lung_seurat, subset = nFeature_RNA > 400 & percent.mito < 15)

DefaultAssay(lung_seurat) <- "RNA"
lung_seurat <- NormalizeData(lung_seurat)
lung_seurat <- FindVariableFeatures(lung_seurat)
lung_seurat <- ScaleData(lung_seurat)
lung_seurat <- RunPCA(lung_seurat)

lung_seurat.harmony <- RunHarmony(lung_seurat, group.by.vars="orig.ident", plot_convergence=TRUE, lambda=2)
lung_harmony_embeddings <- Embeddings(lung_seurat.harmony, 'harmony')
lung_harmony_embeddings[1:5,1:5]
lung_seurat.harmony <- lung_seurat.harmony %>% RunUMAP(dims=pc.num, reduction="harmony")
lung_seurat.harmony <- lung_seurat.harmony %>% FindNeighbors(dims=pc.num) %>% FindClusters(resolution=0.6)

p1 <- DimPlot(lung_seurat.harmony, reduction="harmony", pt.size=0.1, group.by="sampleid")
p2 <- VlnPlot(lung_seurat.harmony, features="harmony_1", group.by="sampleid")
plot_grid(p1, p2)
FeaturePlot(lung_seurat.harmony, c("Spp1","Trem2"))
DimPlot(lung_seurat.harmony, group.by="orig.ident") + plot_annotation(title="lung after integration (Harmony)")
DimPlot(lung_seurat.harmony, group.by="orig.ident", split.by="orig.ident") + plot_annotation(title="lung after integration (Harmony)")

saveRDS(lung_seurat.harmony, file=file.path(harmony_dir, "lung_immuneblood_afterharmony.rds"))

# =============================
# 6. Merge and integrate KIDNEY samples (full code retained)
# =============================
mk_immune_uuo7 <- readRDS(file.path(output_dir, "mk_immune_uuo7_scrublet_labeled.rds"))
mk_immune_uninj <- readRDS(file.path(output_dir, "mk_immune_uninj_scrublet_labeled.rds"))
mk_blood_uuo7 <- readRDS(file.path(output_dir, "mk_blood_uuo7_scrublet_labeled.rds"))
ml_blood_uninj <- readRDS(file.path(output_dir, "ml_blood_uninj_scrublet_labeled.rds"))
mk_blood_uninj <- readRDS(file.path(output_dir, "mk_blood_uninj_scrublet_labeled.rds"))

kidney_seurat <- merge(mk_immune_uuo7, mk_immune_uninj, add.cell.ids = c("mk_immune_uuo7","mk_immune_uninj"))
kidney_seurat <- merge(kidney_seurat, mk_blood_uuo7, add.cell.ids = c("","mk_blood_uuo7"))
kidney_seurat <- merge(kidney_seurat, ml_blood_uninj, add.cell.ids = c("","ml_blood_uninj"))
kidney_seurat <- merge(kidney_seurat, mk_blood_uninj, add.cell.ids = c("","mk_blood_uninj"))

VlnPlot(kidney_seurat, features=c("nFeature_RNA","nCount_RNA","percent.mito"), group.by="sampleid")
kidney_seurat <- subset(kidney_seurat, subset = nFeature_RNA > 400 & percent.mito < 15)

DefaultAssay(kidney_seurat) <- "RNA"
kidney_seurat <- NormalizeData(kidney_seurat)
kidney_seurat <- FindVariableFeatures(kidney_seurat)
kidney_seurat <- ScaleData(kidney_seurat)
kidney_seurat <- RunPCA(kidney_seurat)

kidney_seurat.harmony <- RunHarmony(kidney_seurat, group.by.vars="orig.ident", plot_convergence=TRUE, lambda=6, theta=0)
kidney_harmony_embeddings <- Embeddings(kidney_seurat.harmony, 'harmony')
kidney_harmony_embeddings[1:5,1:5]
kidney_seurat.harmony <- kidney_seurat.harmony %>% RunUMAP(dims=pc.num, reduction="harmony")
kidney_seurat.harmony <- kidney_seurat.harmony %>% FindNeighbors(dims=pc.num) %>% FindClusters(resolution=0.6)

p1 <- DimPlot(kidney_seurat.harmony, reduction="harmony", pt.size=0.1, group.by="sampleid")
p2 <- VlnPlot(kidney_seurat.harmony, features="harmony_1", group.by="sampleid")
plot_grid(p1, p2)
FeaturePlot(kidney_seurat.harmony, c("Ccl8","Mrc1"), split.by="condition")
DimPlot(kidney_seurat.harmony, group.by="orig.ident") + plot_annotation(title="kidney after integration (Harmony)")
DimPlot(kidney_seurat.harmony, group.by="seurat_clusters", split.by="condition") + plot_annotation(title="kidney after integration (Harmony)")
DimPlot(kidney_seurat.harmony, group.by="orig.ident", split.by="orig.ident") + plot_annotation(title="kidney after integration (Harmony)")

saveRDS(kidney_seurat.harmony, file=file.path(harmony_dir, "kidney_immuneblood_afterharmony.rds"))
