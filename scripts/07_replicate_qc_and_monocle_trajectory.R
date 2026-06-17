#!/usr/bin/env Rscript

# Supplementary Code 7
# Replicate-level correlation QC and Monocle trajectory analysis.
# This script integrates additional analyses from the later code notes, including
# CK replicate correlation checks and cell-type-wise pseudotime inference.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(monocle)
  library(RColorBrewer)
  library(grid)
  library(ggsci)
})

source("config/config.R")
set.seed(135)
outdir <- file.path(result_dir, "monocle")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

obj <- readRDS(config$seurat_rds)
DefaultAssay(obj) <- "RNA"

# -----------------------------------------------------------------------------
# 1. Replicate correlation QC
# -----------------------------------------------------------------------------
# The original code compared average log-normalized expression between CK
# replicates. Here the sample names are read from metadata rather than hard-coded.
meta <- read.csv(config$sample_metadata, stringsAsFactors = FALSE)
ck_samples <- meta$sample_id[meta$condition == "CK"]

if (length(ck_samples) >= 2) {
  expr_list <- lapply(ck_samples[1:2], function(s) {
    rowMeans(GetAssayData(obj, slot = "data")[, obj$sample_id == s, drop = FALSE])
  })
  cor_df <- data.frame(rep1 = expr_list[[1]], rep2 = expr_list[[2]])
  cor_test <- cor.test(cor_df$rep1, cor_df$rep2, method = "pearson")
  write.csv(data.frame(
    sample_1 = ck_samples[1], sample_2 = ck_samples[2],
    pearson_r = unname(cor_test$estimate), p_value = cor_test$p.value
  ), "CK_replicate_expression_correlation.csv", row.names = FALSE)

  pdf("CK_replicate_expression_correlation.pdf", width = 5, height = 5)
  print(ggplot(cor_df, aes(rep1, rep2)) +
          geom_point(size = 0.4, alpha = 0.5) +
          theme_bw() +
          labs(x = ck_samples[1], y = ck_samples[2],
               title = paste0("Pearson R = ", round(cor_test$estimate, 3))))
  dev.off()
}

# -----------------------------------------------------------------------------
# 2. Helper function: Seurat subset to Monocle CellDataSet
# -----------------------------------------------------------------------------
make_monocle_cds <- function(seurat_subset) {
  counts <- GetAssayData(seurat_subset, assay = "RNA", slot = "counts")
  feature_ann <- data.frame(gene_id = rownames(counts), gene_short_name = rownames(counts), row.names = rownames(counts))
  sample_ann <- seurat_subset@meta.data
  sample_ann <- sample_ann[colnames(counts), , drop = FALSE]

  cds <- newCellDataSet(
    as.matrix(counts),
    phenoData = new("AnnotatedDataFrame", data = sample_ann),
    featureData = new("AnnotatedDataFrame", data = feature_ann),
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size()
  )
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  cds
}

run_monocle_for_celltype <- function(celltype, obj, top_n_ordering_genes = 500) {
  message("Running Monocle for cell type: ", celltype)
  sub <- subset(obj, subset = CellAnno_Major == celltype)
  if (ncol(sub) < 100 || length(unique(sub$seurat_clusters)) < 2) {
    message("Skipped ", celltype, ": not enough cells or cluster diversity.")
    return(NULL)
  }

  cds <- make_monocle_cds(sub)
  diff <- differentialGeneTest(cds, fullModelFormulaStr = "~seurat_clusters", cores = 1)
  deg <- diff %>% filter(qval < 0.01) %>% arrange(qval)
  if (nrow(deg) < 10) return(NULL)

  ordering_genes <- rownames(deg)[seq_len(min(top_n_ordering_genes, nrow(deg)))]
  cds <- setOrderingFilter(cds, ordering_genes)
  cds <- reduceDimension(cds, max_components = 2, reduction_method = "DDRTree", num_dim = 15)
  cds <- orderCells(cds)

  safe_name <- gsub("[^A-Za-z0-9_]+", "_", celltype)

  pdf(paste0(safe_name, "_ordering_genes.pdf"), width = 5, height = 4)
  print(plot_ordering_genes(cds))
  dev.off()

  p1 <- plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.3)
  p2 <- plot_cell_trajectory(cds, color_by = "CellAnno_Major", cell_size = 0.3)
  p3 <- plot_cell_trajectory(cds, color_by = "seurat_clusters", cell_size = 0.3)
  p4 <- plot_cell_trajectory(cds, color_by = "condition", cell_size = 0.3)

  pdf(paste0(safe_name, "_trajectory_summary.pdf"), width = 8, height = 8)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2)))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  print(p1, vp = vplayout(1, 1)); print(p2, vp = vplayout(1, 2))
  print(p3, vp = vplayout(2, 1)); print(p4, vp = vplayout(2, 2))
  dev.off()

  saveRDS(cds, paste0("monocle_cds_", safe_name, ".rds"))
  invisible(cds)
}

# Run trajectory analysis separately for each annotated major cell type.
celltypes <- sort(unique(obj$CellAnno_Major))
celltypes <- setdiff(celltypes, c(NA, "Unknown"))
cds_list <- lapply(celltypes, run_monocle_for_celltype, obj = obj)
names(cds_list) <- celltypes

# -----------------------------------------------------------------------------
# 3. Candidate gene visualization on a selected trajectory object
# -----------------------------------------------------------------------------
# These candidate genes are retained from the original notes. The block runs only
# if a valid CellDataSet is available.
candidate_genes <- data.frame(
  symbol = c("ECT2", "EIF3A", "CHS", "UGE3", "FIB1", "TT7"),
  gene_id = c("FT01Gene02068", "FT01Gene21048", "FT01Gene11369", "FT01Gene04660", "FT01Gene18447", "FT01Gene29976")
)

first_cds <- cds_list[[which(!sapply(cds_list, is.null))[1]]]
if (!is.null(first_cds)) {
  for (i in seq_len(nrow(candidate_genes))) {
    gene_id <- candidate_genes$gene_id[i]
    if (!gene_id %in% rownames(first_cds)) next
    pData(first_cds)[[gene_id]] <- log2(exprs(first_cds)[gene_id, ] + 1)
    pdf(paste0(candidate_genes$symbol[i], "_", gene_id, "_monocle_expression.pdf"), width = 7, height = 7)
    print(plot_cell_trajectory(first_cds, color_by = gene_id) + scale_color_gsea() +
            ggtitle(paste(candidate_genes$symbol[i], gene_id)))
    dev.off()
  }
}

writeLines(capture.output(sessionInfo()), "sessionInfo_07_monocle.txt")
