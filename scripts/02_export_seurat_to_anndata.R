#!/usr/bin/env Rscript

# Supplementary Code 2
# Export the annotated Seurat object to AnnData format for RNA velocity analysis.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(sceasy)
})

source("config/config.R")
dir.create(config$velocity_outdir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(config$seurat_rds)
obj <- UpdateSeuratObject(obj)
DefaultAssay(obj) <- "RNA"

out_h5ad <- file.path(config$velocity_outdir, "combined_integrated_annotated.h5ad")
sceasy::convertFormat(obj, from = "seurat", to = "anndata", outFile = out_h5ad)

writeLines(capture.output(sessionInfo()), file.path(config$velocity_outdir, "sessionInfo_02_export.txt"))
