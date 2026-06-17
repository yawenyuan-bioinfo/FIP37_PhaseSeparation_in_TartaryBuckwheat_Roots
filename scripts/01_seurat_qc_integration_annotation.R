#!/usr/bin/env Rscript

# Supplementary Code 1
# Seurat preprocessing, sample-level QC, doublet removal, integration, clustering,
# and marker-based annotation for tartary buckwheat root snRNA-seq.

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

source("config/config.R")
set.seed(config$random_seed)
setwd(config$seurat_outdir)

# Root marker genes used for manual annotation and dot-plot visualization.
root_markers <- list(
  Cap = c("FT01Gene00363", "FT01Gene12258", "FT01Gene28886", "FT01Gene33822", "FT01Gene26089"),
  QC = c("FT01Gene01386", "FT01Gene05103", "FT01Gene05160", "FT01Gene28016", "FT01Gene30434"),
  Epidermis = c("FT01Gene07238", "FT01Gene15393", "FT01Gene08384"),
  Cortex = c("FT01Gene06602", "FT01Gene29523"),
  Endodermis = c("FT01Gene00010", "FT01Gene00375", "FT01Gene16897", "FT01Gene28520"),
  Stele = c("FT01Gene27095", "FT01Gene32607", "FT01Gene14654"),
  Pericycle = c("FT01Gene03156", "FT01Gene12718", "FT01Gene16331", "FT01Gene21385"),
  Xylem = c("FT01Gene13753", "FT01Gene23546"),
  Phloem = c("FT01Gene00929", "FT01Gene01294", "FT01Gene10703", "FT01Gene18257", "FT01Gene12523"),
  CasparianStrip = c("FT01Gene16897", "FT01Gene00375"),
  Trichoblast = c("FT01Gene25819", "FT01Gene05408", "FT01Gene26575"),
  Atrichoblast = c("FT01Gene31029", "FT01Gene26329", "FT01Gene13506", "FT01Gene13504", "FT01Gene12209", "FT01Gene09698", "FT01Gene06645", "FT01Gene03056")
)

read_sample_metadata <- function(path) {
  meta <- read.csv(path, stringsAsFactors = FALSE)
  stopifnot(all(c("sample_id", "condition", "replicate", "cellranger_path") %in% colnames(meta)))
  meta
}

run_single_sample <- function(sample_row) {
  sample_name <- sample_row$sample_id
  message("Processing sample: ", sample_name)
  matrix_dir <- file.path(sample_row$cellranger_path, "filtered_feature_bc_matrix")

  counts <- Read10X(matrix_dir)
  obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
  obj$sample_id <- sample_name
  obj$condition <- sample_row$condition
  obj$replicate <- sample_row$replicate

  # QC filtering follows the thresholds used in the original analysis.
  pdf(paste0("QC_violin_before_filter_", sample_name, ".pdf"), width = 6, height = 4)
  print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0.1))
  dev.off()

  obj <- subset(obj, subset = nFeature_RNA > config$min_features &
                  nFeature_RNA < config$max_features &
                  nCount_RNA > config$min_counts &
                  nCount_RNA < config$max_counts)

  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = config$n_variable_features)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = config$sample_npcs)

  pdf(paste0("PCA_elbow_", sample_name, ".pdf"), width = 6, height = 4)
  print(ElbowPlot(obj, ndims = config$sample_npcs))
  dev.off()

  obj <- RunUMAP(obj, dims = seq_len(config$sample_dims))
  obj <- RunTSNE(obj, dims = seq_len(config$sample_dims))
  obj <- FindNeighbors(obj, reduction = "pca", dims = seq_len(config$sample_dims))
  obj <- FindClusters(obj, resolution = config$sample_resolution)

  # DoubletFinder: expected doublets are estimated using the empirical rule from the original workflow.
  sweep_res <- paramSweep_v3(obj, PCs = seq_len(config$sample_dims), sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$MeanBC)]))
  homotypic_prop <- modelHomotypic(obj$seurat_clusters)
  nExp_poi <- round(0.000008 * ncol(obj)^2)
  nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))

  obj <- doubletFinder_v3(obj, PCs = seq_len(config$sample_dims), pN = 0.25, pK = pK,
                          nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  obj <- doubletFinder_v3(obj, PCs = seq_len(config$sample_dims), pN = 0.25, pK = pK,
                          nExp = nExp_poi_adj,
                          reuse.pANN = paste0("pANN_0.25_", pK, "_", nExp_poi), sct = FALSE)

  pdf(paste0("DoubletFinder_summary_", sample_name, ".pdf"), width = 12, height = 4)
  df_cols <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
  plots <- c(list(DimPlot(obj, reduction = "tsne", group.by = "seurat_clusters", label = TRUE) + ggtitle("Clusters")),
             lapply(df_cols, function(col) DimPlot(obj, reduction = "tsne", group.by = col) + ggtitle(col)))
  print(wrap_plots(plots))
  dev.off()
  obj
}

sample_meta <- read_sample_metadata(config$sample_metadata)
obj_list <- lapply(seq_len(nrow(sample_meta)), function(i) run_single_sample(sample_meta[i, ]))
names(obj_list) <- sample_meta$sample_id

# Remove cells classified as doublets by any DoubletFinder classification column.
obj_list <- lapply(obj_list, function(obj) {
  df_cols <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
  doublet_cells <- rownames(obj@meta.data)[rowSums(obj@meta.data[, df_cols, drop = FALSE] == "Doublet") >= 1]
  obj[, !colnames(obj) %in% doublet_cells]
})

features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = config$n_variable_features)
anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features)
combined <- IntegrateData(anchorset = anchors)
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = config$integrated_npcs, verbose = FALSE)

pdf("PCA_elbow_integrated.pdf", width = 6, height = 4)
print(ElbowPlot(combined, ndims = config$integrated_npcs))
dev.off()

combined <- RunUMAP(combined, dims = seq_len(config$integrated_dims))
combined <- RunTSNE(combined, dims = seq_len(config$integrated_dims))
combined <- FindNeighbors(combined, dims = seq_len(config$integrated_dims))
combined <- FindClusters(combined, resolution = config$integrated_resolution)

DefaultAssay(combined) <- "RNA"
markers <- FindAllMarkers(combined, only.pos = FALSE)
write.csv(markers, "RNA_wilcox_all_cluster_markers.csv", row.names = FALSE)

# Add manual annotation from metadata/cluster_annotation.csv.
anno <- read.csv(config$cluster_annotation, stringsAsFactors = FALSE)
combined$CellAnno_Major <- anno$cell_type[match(as.character(combined$seurat_clusters), as.character(anno$cluster))]
combined$CellAnno_Major[is.na(combined$CellAnno_Major)] <- "Unknown"

pdf("UMAP_by_sample.pdf", width = 8, height = 5)
print(DimPlot(combined, reduction = "umap", split.by = "sample_id", label = TRUE))
dev.off()

pdf("UMAP_by_cell_type.pdf", width = 6, height = 5)
print(DimPlot(combined, reduction = "umap", group.by = "CellAnno_Major", label = TRUE, pt.size = 0.1))
dev.off()

pdf("DotPlot_root_marker_genes.pdf", width = 10, height = 12)
print(DotPlot(combined, features = root_markers, group.by = "CellAnno_Major") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
              strip.background = element_rect(fill = "white")) +
        scale_colour_gradient2(low = "white", mid = "lightgrey", high = "firebrick", midpoint = -0.5))
dev.off()

saveRDS(combined, config$seurat_rds)
writeLines(capture.output(sessionInfo()), "sessionInfo_01_seurat.txt")
