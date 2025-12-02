# ============================================================
# Main_Part1.R – Ambient RNA Removal using SoupX
# Author: Lily Jia
# Description:
#     This script removes ambient RNA contamination from
#     raw 10x single-cell data using SoupX.
#     Applies to 12 immune & blood samples from liver, lung, kidney.
# ============================================================

# 1. Load required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
  library(dplyr)
  library(patchwork)
})

# 2. Set paths
data_dir <- "../data/raw_10x/"
output_dir <- "../data/SoupX_corrected/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 3. Define sample list
samples <- list(
  ml_immune_ccl4     = "ml_immune_ccl4",
  ml_immune_uninjC   = "ml_immune_uninjC",
  ml_immune_uninjD   = "ml_immune_uninjD",
  ml_blood_ccl4      = "ml_blood_ccl4",
  ml_blood_uninj     = "ml_blood_uninj",
  mlung_immune_bleo  = "mlung_immune_bleo",
  mlung_immune_uninj = "mlung_immune_uninj",
  mlung_blood_bleo   = "mlung_blood_bleo",
  mk_blood_uuo7      = "mk_blood_uuo7",
  mk_blood_uninj     = "mk_blood_uninj",
  mk_immune_uuo7     = "mk_immune_uuo7",
  mk_immune_uninj    = "mk_immune_uninj"
)

# 4. Process each sample
for (sample_name in names(samples)) {
  sample_id <- samples[[sample_name]]
  message("Processing: ", sample_id)
  
  raw_path <- file.path(data_dir, sample_id, "raw_feature_bc_matrix")
  filt_path <- file.path(data_dir, sample_id, "filtered_feature_bc_matrix")
  
  # Read raw and filtered matrices
  tod <- Read10X(data.dir = raw_path)
  toc <- Read10X(data.dir = filt_path)
  
  # Create SoupChannel object
  sc <- SoupChannel(tod = tod, toc = toc)
  
  # Create Seurat object for clustering
  srat <- CreateSeuratObject(toc)
  srat$percent.mito <- PercentageFeatureSet(srat, pattern = "^mt-")
  srat <- SCTransform(srat, method = "glmGamPoi", vars.to.regress = "percent.mito", verbose = FALSE)
  srat <- RunPCA(srat, verbose = FALSE)
  srat <- RunUMAP(srat, dims = 1:50, verbose = FALSE)
  srat <- FindNeighbors(srat, dims = 1:50, verbose = FALSE)
  srat <- FindClusters(srat, verbose = TRUE)
  
  # Provide clustering info to SoupX
  meta <- srat@meta.data
  umap <- Embeddings(srat, "umap")
  sc <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
  sc <- setDR(sc, umap)
  
  # Estimate contamination and adjust
  sc <- autoEstCont(sc)
  out <- adjustCounts(sc)
  
  # Assign object and save
  assign(paste0(sample_name, ".data"), out, envir = .GlobalEnv)
  saveRDS(out, file = file.path(output_dir, paste0(sample_id, "_soupx_corrected.rds")))
  
  message("Saved: ", sample_id, "_soupx_corrected.rds")
}
