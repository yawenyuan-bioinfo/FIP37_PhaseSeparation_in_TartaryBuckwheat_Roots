#!/usr/bin/env Rscript

# Supplementary Code 8
# Aquaporin/candidate gene heatmaps and FIP37-related expression summaries.
# This script converts the ad hoc heatmap code into a metadata-driven workflow.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
})

source("config/config.R")
outdir <- file.path(result_dir, "candidate_heatmaps")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

obj <- readRDS(config$seurat_rds)
DefaultAssay(obj) <- "RNA"

# FIP37-related candidate genes retained from the original code notes.
fip37_genes <- c("FT01Gene07011", "FT01Gene08370")
fip37_genes <- intersect(fip37_genes, rownames(obj))

if (length(fip37_genes) > 0) {
  pdf("FIP37_candidate_featureplots_by_condition.pdf", width = 12, height = 4)
  print(FeaturePlot(obj, features = fip37_genes, split.by = "condition", ncol = length(unique(obj$condition)),
                    pt.size = 0.05, order = TRUE, reduction = "umap"))
  dev.off()

  obj <- AddModuleScore(obj, features = list(FIP37_signature = fip37_genes), name = "FIP37")
  pdf("FIP37_signature_violin_by_celltype_condition.pdf", width = 12, height = 5)
  print(VlnPlot(obj, features = "FIP371", group.by = "CellAnno_Major", split.by = "condition", pt.size = 0.05) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
}

# Aquaporin Arabidopsis IDs from the original notes. If a feature map is supplied,
# AT IDs are mapped to tartary buckwheat gene IDs before plotting.
aquaporin_at <- c(
  "AT1G31880", "AT1G80760", "AT1G73190", "AT2G45960", "AT3G06100", "AT5G47450",
  "AT3G53420", "AT2G36830", "AT2G37170", "AT2G37180", "AT4G35100", "AT2G29870",
  "AT1G01620", "AT3G61430", "AT3G54820", "AT1G17810", "AT3G47440", "AT2G16850",
  "AT2G39010", "AT3G16240", "AT1G52180", "AT4G23400", "AT2G25810", "AT4G00430",
  "AT5G37810", "AT5G37820", "AT4G17340", "AT4G10380", "AT2G34390", "AT4G01470",
  "AT3G26520", "AT4G18910", "AT5G60660", "AT4G19030", "AT2G21020"
)

map_at_to_ft <- function(at_genes, feature_map_file) {
  if (!file.exists(feature_map_file)) return(character(0))
  fmap <- read.csv(feature_map_file, header = FALSE, stringsAsFactors = FALSE)
  # Expected columns: FT gene ID, AT gene ID, optional symbols/annotations.
  mapped <- fmap$V1[fmap$V2 %in% at_genes]
  unique(mapped[mapped %in% rownames(obj)])
}

aqp_ft <- map_at_to_ft(aquaporin_at, config$feature_map)
if (length(aqp_ft) > 0) {
  expr_df <- FetchData(obj, vars = aqp_ft)
  expr_df$condition <- obj$condition
  expr_df$CellAnno_Major <- obj$CellAnno_Major

  avg_df <- expr_df %>%
    group_by(condition, CellAnno_Major) %>%
    summarise(across(all_of(aqp_ft), mean, na.rm = TRUE), .groups = "drop")

  mat <- avg_df %>%
    unite("group", condition, CellAnno_Major, sep = "_") %>%
    tibble::column_to_rownames("group") %>%
    as.matrix()
  mat <- t(mat)

  pdf("Aquaporin_average_expression_heatmap.pdf", width = 8, height = max(4, length(aqp_ft) * 0.18))
  print(pheatmap(mat, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE,
                 color = colorRampPalette(c("#009DFF", "white", "#D73A34"))(50),
                 fontsize_row = 6, fontsize_col = 6))
  dev.off()
}

writeLines(capture.output(sessionInfo()), "sessionInfo_08_candidate_heatmaps.txt")
