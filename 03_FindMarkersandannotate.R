# Load packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(clustree)
library(scales)
library(SingleR)
library(celldex)
library(data.table)
library(stringr)

# Define input/output directories
input_dir <- "../data/harmony/"
output_dir <- "../results/markers/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ===================== #
#      LIVER ANALYSIS   #
# ===================== #

liver_harmony <- readRDS(file.path(input_dir, "liver_immuneblood_afterharmony.rds"))
ElbowPlot(liver_harmony, ndims = 50)
pc.num <- 1:50
liver_harmony <- FindNeighbors(liver_harmony, dims = pc.num)
liver_harmony <- FindClusters(liver_harmony, resolution = c(0.3, 0.5, 0.6, 0.7, 1.0))
clustree(liver_harmony@meta.data, prefix = "RNA_snn_res.")
liver_harmony <- FindClusters(liver_harmony, resolution = 0.6)

liver.markers <- FindAllMarkers(liver_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(liver.markers, file = file.path(output_dir, "liver_allmarkers_res0.6.csv"))

top3 <- liver.markers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
pdf(file.path(output_dir, "liver_heatmap.pdf"), width = 10, height = 12)
DoHeatmap(liver_harmony, features = as.character(top3$gene)) + NoLegend()
dev.off()

genes_to_check <- c("Pdgfrb", "Pdfgra", "Fbln2", "Cd34", "Mfap4", "Reln", "Ecm1", "Rgs5", "Vipr1", "Ngfr",
                    "Adamtsl2", "Col1a1", "Col3a1", "Timp1", "Tagln", "Myh11", "Pdpn", "Dmkn", "Egln3", "Ptprc",
                    "Lyz2", "Cd68", "Ly6c2", "Ccr2", "Treml4", "Ace", "Clec4f", "Timd4", "Marco", "Trem2", "Cd9",
                    "Spp1", "Cd63", "Cd14", "Mmp13", "Cxcl10", "Ifitm3", "Arg1", "Mmp12", "Cd209a", "Cd74", "Xcr1",
                    "Cd79a", "Cd79b", "Siglech", "Cd3d", "Cd3g", "Nkg7", "Hbb-bs", "Hba-a1", "Ly6g", "Fcer1a",
                    "Cpa3", "Pecam1", "Kdr", "Mki67", "Top2a", "Msln", "Cpe", "Slc22a1", "Trdn", "Clec4g", "Rspo3", "Gja5")
p_all_markers <- DotPlot(liver_harmony, features = genes_to_check) + coord_flip()
print(p_all_markers)


# ===================== #
#      LUNG ANALYSIS    #
# ===================== #

lung_harmony <- readRDS(file.path(input_dir, "lung_immuneblood_afterharmony.rds"))
ElbowPlot(lung_harmony, ndims = 50)
pc.num <- 1:50
lung_harmony <- FindNeighbors(lung_harmony, dims = pc.num)
lung_harmony <- FindClusters(lung_harmony, resolution = c(0.3, 0.5, 0.6, 0.7, 1.0))
clustree(lung_harmony@meta.data, prefix = "RNA_snn_res.")
lung_harmony <- FindClusters(lung_harmony, resolution = 0.6)

lung.markers <- FindAllMarkers(lung_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(lung.markers, file = file.path(output_dir, "lung_allmarkers_res0.6.csv"))

top3 <- lung.markers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
pdf(file.path(output_dir, "lung_heatmap.pdf"), width = 10, height = 12)
DoHeatmap(lung_harmony, features = as.character(top3$gene)) + NoLegend()
dev.off()

p_all_markers <- DotPlot(lung_harmony, features = genes_to_check) + coord_flip()
print(p_all_markers)


# ===================== #
#     KIDNEY ANALYSIS   #
# ===================== #

kidney_harmony <- readRDS(file.path(input_dir, "kidney_immuneblood_afterharmony.rds"))
ElbowPlot(kidney_harmony, ndims = 50)
pc.num <- 1:50
kidney_harmony <- FindNeighbors(kidney_harmony, dims = pc.num)
kidney_harmony <- FindClusters(kidney_harmony, resolution = c(0.3, 0.5, 0.6, 0.7, 1.0))
clustree(kidney_harmony@meta.data, prefix = "RNA_snn_res.")
kidney_harmony <- FindClusters(kidney_harmony, resolution = 0.6)

kidney.markers <- FindAllMarkers(kidney_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(kidney.markers, file = file.path(output_dir, "kidney_allmarkers_res0.6.csv"))

top3 <- kidney.markers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
pdf(file.path(output_dir, "kidney_heatmap.pdf"), width = 10, height = 12)
DoHeatmap(kidney_harmony, features = as.character(top3$gene)) + NoLegend()
dev.off()

p_all_markers <- DotPlot(kidney_harmony, features = genes_to_check) + coord_flip()
print(p_all_markers)

# =============================================== #
#      Organ-specific celltype & lineage annotation
# =============================================== #

# ============================================================
#                       LIVER ANNOTATION
# ============================================================

liver_harmony <- readRDS(file.path(input_dir, "liver_immuneblood_afterharmony.rds"))

# ---------------- Annotation labels ---------------- #
cluster.annotation <- c(
  "Liver_Ly6Clo_mono","Liver_KC1","Liver_Ly6Chi_mono1","Liver_Bcell1",
  "Liver_Neutrophil1","Liver_SAMP1","Liver_Tcell1","Liver_NKcell1",
  "Liver_pdc1","Liver_Tcell2","Liver_Neutrophil2","Liver_cDC1",
  "Liver_cDC2","Liver_Tcell3","Liver_Endo","Liver_Proliferating1",
  "Liver_NKcell2","Liver_Bcell2","Liver_Ly6Chi_mono2","Liver_MP/Bcell_doublet",
  "Liver_KC2","Liver_SAMP2","Liver_pdc2","Liver_MP1","Liver_Basophil",
  "Liver_Proliferating2"
)

object <- liver_harmony
names(cluster.annotation) <- levels(object)
object <- RenameIdents(object, cluster.annotation)

# ---------------- Add celltype metadata ---------------- #
cell.data <- data.table(
  barcode = colnames(object),
  celltype = Idents(object)
)
cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$barcode <- NULL
object <- AddMetaData(object, cell.data, col.name = "celltype")

# ---------------- Colors ---------------- #
my_cols <- c(
  'Liver_Ly6Clo_mono'='#D575FE','Liver_KC1'='#EB69F0','Liver_Ly6Chi_mono1'='#FF65AC',
  'Liver_Bcell1'='#24B700','Liver_Neutrophil1'='#BE9C00','Liver_SAMP1'='#8B93FF',
  'Liver_Tcell1'='#00BB49','Liver_NKcell1'='#00BE70','Liver_pdc1'='#F8766D',
  'Liver_Tcell2'='#00C090','Liver_Neutrophil2'='#8CAB00','Liver_cDC1'='#00ACFC',
  'Liver_cDC2'='#42A0FF','Liver_Tcell3'='#00C1AB','Liver_Endo'='#E18A00',
  'Liver_Proliferating1'='#68B100','Liver_NKcell2'='#00BFC4','Liver_Bcell2'='#00BBDA',
  'Liver_Ly6Chi_mono2'='#FF6C91','Liver_MP/Bcell_doublet'='#A8A400','Liver_KC2'='#F962DD',
  'Liver_SAMP2'='#B684FF','Liver_pdc2'='#EE8043','Liver_MP1'='#FF61C6',
  'Liver_Basophil'='#D19300','Liver_Proliferating2'='#00B5ED'
)

# ---------------- Plot celltypes ---------------- #
pdf(file.path(output_dir, "liver_celltype_umap.pdf"), width=10, height=8)
print(
  DimPlot(object, reduction="umap", group.by="celltype",
          label=TRUE, cols=my_cols, pt.size=0.5) +
    theme(plot.title=element_text(hjust=0.5))
)
dev.off()

# ---------------- Lineage annotation ---------------- #
lineage.annotation <- c(
  "mononuclearphagocyte","mononuclearphagocyte","mononuclearphagocyte","Lymphoid",
  "Granulocyte","mononuclearphagocyte","Lymphoid","Lymphoid","mononuclearphagocyte",
  "Lymphoid","Granulocyte","mononuclearphagocyte","mononuclearphagocyte","Lymphoid",
  "Endothelia","Proliferating","Lymphoid","Lymphoid","mononuclearphagocyte",
  "Doublet","mononuclearphagocyte","mononuclearphagocyte","mononuclearphagocyte",
  "mononuclearphagocyte","Granulocyte","Proliferating"
)

lineage.data <- data.table(
  celltype = cluster.annotation,
  lineage = lineage.annotation
)

meta.data <- merge(cell.data, lineage.data, by = "celltype")
meta.data <- data.frame(meta.data, row.names = meta.data$barcode)
meta.data$barcode <- NULL
meta.data$celltype <- NULL
object <- AddMetaData(object, meta.data, col.name = "lineage")

pdf(file.path(output_dir, "liver_lineage_umap.pdf"), width=10, height=8)
print(
  DimPlot(object, reduction="umap", group.by="lineage",
          pt.size=0.3, label=TRUE)
)
dev.off()

saveRDS(object, file.path(input_dir, "liver_seurat_annotated.rds"))


# ============================================================
#                       LUNG ANNOTATION
# ============================================================

lung_harmony <- readRDS(file.path(input_dir, "lung_immuneblood_afterharmony.rds"))

# SingleR reference
mouseImmu <- ImmGenData()
mouseRNA  <- MouseRNAseqData()

sce <- lung_harmony
sce_for_SingleR <- GetAssayData(sce, slot="data")
clusters <- sce@meta.data$seurat_clusters

pred.mouseImmu <- SingleR(
  test = sce_for_SingleR,
  ref = mouseImmu,
  labels = mouseImmu$label.main,
  method = "cluster",
  clusters = clusters
)

pred.mouseRNA <- SingleR(
  test = sce_for_SingleR,
  ref = mouseRNA,
  labels = mouseRNA$label.fine,
  method = "cluster",
  clusters = clusters
)

# ---------------- Annotation labels ---------------- #
cluster.annotation <- c(
  "Lung_Bcell1","Lung_Ly6Clo_mono","Lung_Ly6Chi_mono","Lung_Tcell1",
  "Lung_mo_dc","Lung_Ilc.NKcell","Lung_AlveolarMac","Lung_Tcell2",
  "Lung_Neutrophil","Lung_Bcell2","Lung_NKcell","Lung_Tcell3",
  "Lung_Bcell3","Lung_MP1","Lung_mono1","Lung_cDC2","Lung_cDC1",
  "Lung_mono2","Lung_Endo","Lung_Bcell4","Lung_SAMP","Lung_Tcell4",
  "Lung_mono3","Lung_Proliferating","Lung_RBC/Ilc_doublet","Lung_pDC",
  "Lung_Basophil","Lung_Bcell5","Lung_MP2"
)

object <- lung_harmony
names(cluster.annotation) <- levels(object)
object <- RenameIdents(object, cluster.annotation)

# ---------------- Add celltype metadata ---------------- #
cell.data <- data.table(
  barcode = colnames(object),
  celltype = Idents(object)
)
cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$barcode <- NULL
object <- AddMetaData(object, cell.data, col.name="celltype")

# ---------------- Colors ---------------- #
my_cols <- c(
  'Lung_Bcell1'='#F365E6','Lung_Bcell2'='#FC61D4','Lung_Bcell3'='#FF62BE',
  'Lung_Bcell4'='#FF67A6','Lung_Bcell5'='#FF6C91','Lung_Ly6Clo_mono'='#46A0FF',
  'Lung_Ly6Chi_mono'='#E46DF6','Lung_Tcell1'='#B4A000','Lung_Tcell2'='#00ABFD',
  'Lung_Tcell3'='#83AD00','Lung_Tcell4'='#60B200','Lung_Ilc.NKcell'='#F8766D',
  'Lung_mono3'='#00BE6B','Lung_RBC/Ilc_doublet'='#00C088','Lung_SAMP'='#19B700',
  'Lung_mono1'='#00BB48','Lung_mono2'='#C69900','Lung_mo_dc'='#9EA700',
  'Lung_MP1'='#00C1A2','Lung_MP2'='#EF7F48','Lung_Endo'='#E48800',
  'Lung_Proliferating'='#D69100','Lung_Neutrophil'='#00C0BA','Lung_cDC1'='#00BECF',
  'Lung_cDC2'='#00B9E1','Lung_pDC'='#8794FF','Lung_AlveolarMac'='#B086FF',
  'Lung_Basophil'='#CE79FF','Lung_NKcell'='#00B3F1'
)

# ---------------- UMAP plots ---------------- #
pdf(file.path(output_dir, "lung_celltype_umap.pdf"), width=10, height=8)
print(
  DimPlot(object, reduction="umap", group.by="celltype",
          label=TRUE, cols=my_cols, pt.size=0.5)
)
dev.off()

pdf(file.path(output_dir, "lung_celltype_split_condition.pdf"), width=12, height=6)
print(
  DimPlot(object, reduction="umap", group.by="celltype",
          cols=my_cols, pt.size=0.4, split.by="condition")
)
dev.off()

# ---------------- Lineage annotation ---------------- #
lineage.annotation <- c(
  "Lymphoid","mononuclearphagocyte","mononuclearphagocyte","Lymphoid",
  "mononuclearphagocyte","Lymphoid","mononuclearphagocyte","Lymphoid",
  "Granulocyte","Lymphoid","Lymphoid","Lymphoid","Lymphoid","mononuclearphagocyte",
  "mononuclearphagocyte","mononuclearphagocyte","mononuclearphagocyte",
  "mononuclearphagocyte","Endothelia","Lymphoid","mononuclearphagocyte",
  "Lymphoid","mononuclearphagocyte","Proliferating","Doublet","mononuclearphagocyte",
  "Granulocyte","Lymphoid","mononuclearphagocyte"
)

lineage.data <- data.table(celltype=cluster.annotation, lineage=lineage.annotation)
meta.data <- merge(cell.data, lineage.data, by="celltype")
meta.data <- data.frame(meta.data, row.names=meta.data$barcode)
meta.data$barcode <- NULL
meta.data$celltype <- NULL
object <- AddMetaData(object, meta.data, col.name="lineage")

pdf(file.path(output_dir, "lung_lineage_umap.pdf"), width=10, height=8)
print(
  DimPlot(object, reduction="umap", group.by="lineage",
          pt.size=0.3, label=TRUE)
)
dev.off()

saveRDS(object, file.path(input_dir, "lung_seurat_annotated.rds"))


# ============================================================
#                     KIDNEY ANNOTATION
# ============================================================

kidney_harmony <- readRDS(file.path(input_dir, "kidney_immuneblood_afterharmony.rds"))

mouseImmu <- ImmGenData()
mouseRNA  <- MouseRNAseqData()

sce <- kidney_harmony
sce_for_SingleR <- GetAssayData(sce, slot="data")
clusters <- sce@meta.data$seurat_clusters

pred.mouseImmu <- SingleR(
  test = sce_for_SingleR,
  ref = mouseImmu,
  labels = mouseImmu$label.main,
  method = "cluster",
  clusters = clusters
)

pred.mouseRNA <- SingleR(
  test = sce_for_SingleR,
  ref = mouseRNA,
  labels = mouseRNA$label.fine,
  method = "cluster",
  clusters = clusters
)

cluster.annotation <- c(
  "Kid_Bcell1","Kid_Ly6Chi_mono","Kid_Mrc1+Mac","Kid_Patrollingmono",
  "Kid_Mmp12+Mac","Kid_Neutrophil","Kid_Tcell1","Kid_iLC","Kid_Tcell2",
  "Kid_cDC2","Kid_pDC","Kid_Bcell2","Kid_Bcell3","Kid_mono1","Kid_NKcell",
  "Kid_Tcell3","Kid_Monocyte_drivedBcell","Kid_cDC1","Kid_Proliferating1",
  "Kid_IFNinduceMac","Kid_Endo","Kid_Tcell4","Kid_Basophil",
  "Kid_Proliferating2","Kid_Ccr7+DC","Kid_Tcell5","Kid_Ecm1+Mac",
  "Kid_mono2","Kid_RBC"
)

object <- kidney_harmony
names(cluster.annotation) <- levels(object)
object <- RenameIdents(object, cluster.annotation)

cell.data <- data.table(
  barcode = colnames(object),
  celltype = Idents(object)
)
cell.data <- data.frame(cell.data, row.names=cell.data$barcode)
cell.data$barcode <- NULL
object <- AddMetaData(object, cell.data, col.name="celltype")

my_cols <- c(
  'Kid_Tcell1'='#F365E6','Kid_Tcell2'='#FC61D4','Kid_Tcell3'='#FF62BE',
  'Kid_Tcell4'='#FF67A6','Kid_Tcell5'='#FF6C91','Kid_Ly6Chi_mono'='#46A0FF',
  'Kid_Patrollingmono'='#E46DF6','Kid_Monocyte_drivedBcell'='#B4A000',
  'Kid_Bcell1'='#19B700','Kid_Bcell2'='#83AD00','Kid_Bcell3'='#60B200',
  'Kid_iLC'='#F8766D','Kid_RBC'='#00BE6B','Kid_mono1'='#00C088',
  'Kid_Mmp12+Mac'='#00ABFD','Kid_Mrc1+Mac'='#00BB48','Kid_mono2'='#C69900',
  'Kid_cDC2'='#9EA700','Kid_pDC'='#00C1A2','Kid_cDC1'='#EF7F48',
  'Kid_Proliferating1'='#E48800','Kid_IFNinduceMac'='#D69100',
  'Kid_Proliferating2'='#00C0BA','Kid_Ccr7+DC'='#00BECF','Kid_Ecm1+Mac'='#00B9E1',
  'Kid_Neutrophil'='#8794FF','Kid_Endo'='#B086FF','Kid_Basophil'='#CE79FF',
  'Kid_NKcell'='#00B3F1'
)

pdf(file.path(output_dir, "kidney_celltype_umap.pdf"), width=10, height=8)
print(
  DimPlot(object, reduction="umap", group.by="celltype",
          label=TRUE, cols=my_cols, pt.size=0.5)
)
dev.off()

pdf(file.path(output_dir, "kidney_celltype_split_condition.pdf"), width=12, height=6)
print(
  DimPlot(object, reduction="umap", group.by="celltype",
          cols=my_cols, split.by="condition", pt.size=0.4)
)
dev.off()

lineage.annotation <- c(
  "Lymphoid","mononuclearphagocyte","mononuclearphagocyte","mononuclearphagocyte",
  "mononuclearphagocyte","Granulocyte","Lymphoid","Lymphoid","Lymphoid",
  "mononuclearphagocyte","mononuclearphagocyte","Lymphoid","Lymphoid",
  "mononuclearphagocyte","Lymphoid","Lymphoid","Lymphoid","mononuclearphagocyte",
  "Proliferating","mononuclearphagocyte","Endothelia","Lymphoid","Granulocyte",
  "Proliferating","mononuclearphagocyte","Lymphoid","mononuclearphagocyte",
  "mononuclearphagocyte","Redbloodcell"
)

lineage.data <- data.table(celltype=cluster.annotation, lineage=lineage.annotation)
meta.data <- merge(cell.data, lineage.data, by="celltype")
meta.data <- data.frame(meta.data, row.names=meta.data$barcode)
meta.data$barcode <- NULL
meta.data$celltype <- NULL
object <- AddMetaData(object, meta.data, col.name="lineage")

pdf(file.path(output_dir, "kidney_lineage_umap.pdf"), width=10, height=8)
print(
  DimPlot(object, reduction="umap", group.by="lineage",
          pt.size=0.3, label=TRUE)
)
dev.off()

saveRDS(object, file.path(input_dir, "kidney_seurat_annotated.rds"))
