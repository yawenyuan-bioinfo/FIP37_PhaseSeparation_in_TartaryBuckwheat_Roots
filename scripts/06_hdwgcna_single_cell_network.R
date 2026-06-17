#!/usr/bin/env Rscript

# Supplementary Code 6: Cell-type-specific hdWGCNA analysis and Cytoscape export

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(WGCNA)
  library(hdWGCNA)
  library(cowplot)
  library(patchwork)
})

theme_set(theme_cowplot())
set.seed(12345)

# -----------------------------
# 1. User configuration
# -----------------------------
source("config/config.R")
hdwgcna_config <- list(
  seurat_rds = config$seurat_rds,
  outdir = file.path(result_dir, "hdWGCNA"),
  target_celltype = "Xylem",
  target_gene = "FT01Gene07011",
  soft_power = 6,
  min_module_size = 30,
  merge_cut_height = 0.25,
  cytoscape_threshold = 0.2
)

dir.create(hdwgcna_config$outdir, showWarnings = FALSE, recursive = TRUE)
setwd(hdwgcna_config$outdir)

obj <- readRDS(hdwgcna_config$seurat_rds)
DefaultAssay(obj) <- "RNA"

# -----------------------------
# 2. Expression check for target gene
# -----------------------------
if (hdwgcna_config$target_gene %in% rownames(obj)) {
  target_subset <- subset(obj, subset = CellAnno_Major == hdwgcna_config$target_celltype)
  expr_mat <- GetAssayData(target_subset, slot = "data")
  n_expressed <- sum(expr_mat[hdwgcna_config$target_gene, ] > 0)
  n_total <- ncol(expr_mat)
  write.csv(
    data.frame(gene = hdwgcna_config$target_gene, celltype = hdwgcna_config$target_celltype,
               expressed_cells = n_expressed, total_cells = n_total),
    "target_gene_expression_summary.csv",
    row.names = FALSE
  )
}

# -----------------------------
# 3. hdWGCNA preprocessing
# -----------------------------
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)

obj <- SetupForWGCNA(
  obj,
  gene_select = "fraction",
  fraction = 0.1,
  wgcna_name = "root_hdwgcna"
)

obj <- MetacellsByGroups(
  seurat_obj = obj,
  group.by = c("CellAnno_Major"),
  k = 20,
  reduction = "pca",
  slot = "counts",
  max_shared = 10,
  ident.group = "CellAnno_Major"
)

obj <- NormalizeMetacells(obj)

obj <- SetDatExpr(
  obj,
  group_name = hdwgcna_config$target_celltype,
  group.by = "CellAnno_Major",
  assay = "RNA",
  use_metacells = TRUE,
  slot = "data"
)

obj <- TestSoftPowers(obj, powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2)))

obj <- ConstructNetwork(
  obj,
  soft_power = hdwgcna_config$soft_power,
  setDatExpr = FALSE,
  corType = "pearson",
  networkType = "unsigned",
  TOMType = "unsigned",
  minModuleSize = hdwgcna_config$min_module_size,
  mergeCutHeight = hdwgcna_config$merge_cut_height,
  overwrite_tom = TRUE
)

# -----------------------------
# 4. Module eigengenes and connectivity
# -----------------------------
obj <- ScaleData(obj, features = VariableFeatures(obj))
obj <- ModuleEigengenes(obj)
obj <- ModuleConnectivity(obj)

modules <- GetModules(obj)
write.csv(modules, "hdWGCNA_modules.csv", row.names = FALSE)

if (hdwgcna_config$target_gene %in% modules$gene_name) {
  write.csv(subset(modules, gene_name == hdwgcna_config$target_gene), "target_gene_module.csv", row.names = FALSE)
}

# -----------------------------
# 5. Export network to Cytoscape
# -----------------------------
modules_use <- modules %>%
  filter(module != "grey") %>%
  mutate(module = droplevels(module))

genes_use <- unique(modules_use$gene_name)
TOM <- GetTOM(obj)
cur_TOM <- TOM[genes_use, genes_use]
cur_TOM[upper.tri(cur_TOM)] <- t(cur_TOM)[upper.tri(cur_TOM)]

node_info <- data.frame(
  gene = genes_use,
  moduleColor = modules_use$color[match(genes_use, modules_use$gene_name)],
  stringsAsFactors = FALSE
)

exportNetworkToCytoscape(
  cur_TOM,
  nodeNames = genes_use,
  nodeAttr = node_info,
  edgeFile = "CytoscapeEdges.txt",
  nodeFile = "CytoscapeNodes.txt",
  threshold = hdwgcna_config$cytoscape_threshold
)

saveRDS(obj, "hdWGCNA_processed_object.rds")
sessionInfo()
