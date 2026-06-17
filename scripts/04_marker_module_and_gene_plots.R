#!/usr/bin/env Rscript

# Supplementary Code 4: Marker module scores and selected gene visualization

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
})

set.seed(42)

# -----------------------------
# 1. User configuration
# -----------------------------
source("config/config.R")
seurat_rds <- config$seurat_rds
feature_map_csv <- config$feature_map
outdir <- file.path(result_dir, "marker_plots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

obj <- readRDS(seurat_rds)
DefaultAssay(obj) <- "RNA"

# -----------------------------
# 2. Arabidopsis-to-Fagopyrum marker helper
# -----------------------------
find_mapped_genes <- function(seurat_obj, mapping_df, marker_list) {
  out <- list()
  for (group_name in names(marker_list)) {
    query_genes <- marker_list[[group_name]]
    valid_genes <- intersect(query_genes, mapping_df$V2)
    if (length(valid_genes) == 0) next

    mapped <- mapping_df$V1[mapping_df$V2 %in% valid_genes]
    mapped <- mapped[mapped %in% rownames(seurat_obj)]
    if (length(mapped) > 0) out[[group_name]] <- unique(mapped)
  }
  out
}

feature_map <- read.csv(feature_map_csv, header = FALSE)

plant_marker_at <- list(
  cortex_markers = c("AT5G57360", "AT5G47710", "AT2G22890", "AT3G26765"),
  endodermis_markers = c("AT2G20830", "AT5G13320", "AT4G34120", "AT4G34660"),
  xylem_markers = c("AT1G21960", "AT1G71930", "AT4G24290", "AT2G34510"),
  procambium_markers = c("AT4G37650", "AT3G20840", "AT1G51190", "AT5G59730", "AT3G25710", "AT2G31420"),
  root_cap_markers = c("AT2G37620", "AT1G05580", "AT4G22010", "AT1G72260"),
  epidermis_markers = c("AT2G46400", "AT3G47910", "AT5G17790", "AT1G74930"),
  elongation_markers = c("AT2G43750", "AT2G46790", "AT3G13420", "AT2G14710"),
  maturation_markers = c("AT5G15710", "AT3G21710", "AT2G37770", "AT1G02900")
)

plant_markers <- find_mapped_genes(obj, feature_map, plant_marker_at)
plant_markers <- plant_markers[sapply(plant_markers, length) > 0]

# -----------------------------
# 3. Module score plots
# -----------------------------
obj_score <- AddModuleScore(obj, features = plant_markers, name = "Marker")
score_cols <- grep("^Marker", colnames(obj_score@meta.data), value = TRUE)
colnames(obj_score@meta.data)[match(score_cols, colnames(obj_score@meta.data))] <- names(plant_markers)

pdf("03_marker_module_score_featureplots.pdf", width = 12, height = 10)
print(FeaturePlot(obj_score, features = names(plant_markers), ncol = 3, pt.size = 0.05, order = TRUE, reduction = "umap"))
dev.off()

pdf("03_marker_module_score_violin.pdf", width = 10, height = max(6, 0.4 * length(plant_markers)))
print(VlnPlot(obj_score, features = names(plant_markers), ncol = 1, pt.size = 0))
dev.off()

# -----------------------------
# 4. Selected FIP37-related genes
# -----------------------------
selected_genes <- c("FT01Gene07011", "FT01Gene08370")
selected_genes <- selected_genes[selected_genes %in% rownames(obj)]

for (gene in selected_genes) {
  p <- FeaturePlot(obj, features = gene, cols = c("lightgray", "#FFD700", "red"),
                   split.by = "group", ncol = 4, pt.size = 0.05, order = TRUE, reduction = "umap")
  ggsave(paste0("04_featureplot_", gene, "_by_group.pdf"), p, width = 16, height = 5)

  p <- VlnPlot(obj, features = gene, split.by = "group", group.by = "CellAnno_Major", pt.size = 0.1)
  ggsave(paste0("04_violin_", gene, "_by_celltype_group.pdf"), p, width = 18, height = 5)
}

if (length(selected_genes) > 0) {
  obj <- AddModuleScore(obj, features = list(FIP37_related = selected_genes), name = "FIP37")
  p <- VlnPlot(obj, features = "FIP371", split.by = "group", group.by = "CellAnno_Major", pt.size = 0.1)
  ggsave("04_violin_FIP37_related_score.pdf", p, width = 18, height = 5)
}

# -----------------------------
# 5. Mean-expression heatmap
# -----------------------------
if (length(selected_genes) > 0) {
  gene_data <- FetchData(obj, vars = selected_genes)
  gene_data$group <- obj$group
  gene_data$CellAnno_Major <- obj$CellAnno_Major

  gene_avg <- gene_data %>%
    group_by(group, CellAnno_Major) %>%
    summarise(across(all_of(selected_genes), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  gene_avg_wide <- gene_avg %>%
    pivot_wider(names_from = group, values_from = all_of(selected_genes))

  gene_mat <- as.matrix(gene_avg_wide[, -1, drop = FALSE])
  rownames(gene_mat) <- gene_avg_wide$CellAnno_Major

  pdf("04_selected_gene_mean_expression_heatmap.pdf", width = 10, height = 8)
  print(pheatmap(gene_mat, cluster_rows = TRUE, cluster_cols = FALSE, scale = "row",
                 display_numbers = TRUE, fontsize = 10,
                 main = "Selected gene expression"))
  dev.off()
}

saveRDS(obj, "combined_with_marker_scores.rds")
sessionInfo()
