#!/usr/bin/env Rscript

# Supplementary Code 5: Bulk time-course expression summary, WGCNA, and TO-GCN input preparation

suppressPackageStartupMessages({
  library(WGCNA)
  library(dplyr)
  library(pheatmap)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()
set.seed(123)

# -----------------------------
# 1. User configuration
# -----------------------------
source("config/config.R")
bulk_meta <- read.csv(config$bulk_metadata, stringsAsFactors = FALSE)
wgcna_config <- list(
  workdir = config$wgcna_outdir,
  tpm_csv = config$bulk_tpm,
  count_csv = config$bulk_count,
  annotation_csv = file.path(metadata_dir, "go_merged.csv"),
  sample_time = bulk_meta$time,
  soft_power_default = 7,
  min_module_size = 30,
  merge_cut_height = 0.25
)

dir.create(wgcna_config$workdir, showWarnings = FALSE, recursive = TRUE)
setwd(wgcna_config$workdir)

# -----------------------------
# 2. Helper functions
# -----------------------------
compute_power <- function(sft, default_power = 7) {
  if (!is.na(sft$powerEstimate)) return(sft$powerEstimate)
  if (max(sft$fitIndices$SFT.R.sq, na.rm = TRUE) > 0.8) {
    return(sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)])
  }
  default_power
}

plot_soft_threshold <- function(sft, powers, outfile) {
  pdf(outfile, width = 9, height = 5)
  par(mfrow = c(1, 2))
  cex1 <- 0.9
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft threshold power", ylab = "Scale-free topology model fit, signed R^2",
       type = "n", main = "Scale independence")
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, cex = cex1, col = "red")
  abline(h = 0.90, col = "red")

  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab = "Soft threshold power", ylab = "Mean connectivity",
       type = "n", main = "Mean connectivity")
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
  dev.off()
}

# -----------------------------
# 3. Load expression data
# -----------------------------
tpm <- read.csv(wgcna_config$tpm_csv, row.names = 1, check.names = FALSE)
datExpr0 <- as.data.frame(t(tpm))

sample_trait <- data.frame(
  sample = rownames(datExpr0),
  time = wgcna_config$sample_time[seq_len(nrow(datExpr0))]
)
rownames(sample_trait) <- sample_trait$sample

# Remove genes or samples with excessive missing values
qc <- goodSamplesGenes(datExpr0, verbose = 3)
if (!qc$allOK) {
  datExpr0 <- datExpr0[qc$goodSamples, qc$goodGenes]
  sample_trait <- sample_trait[rownames(datExpr0), , drop = FALSE]
}

datExpr <- datExpr0

# -----------------------------
# 4. Sample clustering and soft-power selection
# -----------------------------
sample_tree <- hclust(dist(datExpr), method = "average")
sample_colors <- numbers2colors(as.numeric(factor(sample_trait$time)),
                                colors = rainbow(length(unique(sample_trait$time))),
                                signed = FALSE)

pdf("step1_sample_dendrogram_and_trait.pdf", width = 8, height = 6)
plotDendroAndColors(sample_tree, sample_colors, groupLabels = "time",
                    cex.dendroLabels = 0.8, marAll = c(1, 4, 3, 1),
                    addGuide = TRUE, main = "Sample dendrogram and trait")
dev.off()

powers <- c(1:10, seq(12, 30, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
plot_soft_threshold(sft, powers, "step2_soft_threshold_selection.pdf")
soft_power <- compute_power(sft, default_power = wgcna_config$soft_power_default)
write.csv(data.frame(soft_power = soft_power), "selected_soft_power.csv", row.names = FALSE)

# -----------------------------
# 5. WGCNA network construction
# -----------------------------
net <- blockwiseModules(
  datExpr,
  power = soft_power,
  TOMType = "signed",
  minModuleSize = wgcna_config$min_module_size,
  reassignThreshold = 0,
  mergeCutHeight = wgcna_config$merge_cut_height,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 1
)

module_colors <- labels2colors(net$colors)
gene_module <- data.frame(gene = colnames(datExpr), module = module_colors)
write.csv(gene_module, "WGCNA_gene_module_membership.csv", row.names = FALSE)

MEs0 <- moduleEigengenes(datExpr, module_colors)$eigengenes
MEs <- orderMEs(MEs0)
design <- model.matrix(~ sample_trait$time)
colnames(design) <- c("Intercept", "Time")

module_trait_cor <- cor(MEs, design[, "Time", drop = FALSE], use = "p")
module_trait_pvalue <- corPvalueStudent(module_trait_cor, nrow(datExpr))

text_matrix <- paste0("r=", signif(module_trait_cor, 2), ", p=", signif(module_trait_pvalue, 1))
dim(text_matrix) <- dim(module_trait_cor)

pdf("step4_module_trait_relationship_heatmap.pdf", width = 6, height = 12)
par(mar = c(5, 9, 3, 3))
labeledHeatmap(
  Matrix = module_trait_cor,
  xLabels = "Time",
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = text_matrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)
dev.off()

selected_modules <- rownames(module_trait_pvalue)[module_trait_pvalue[, 1] < 0.05]
write.csv(data.frame(module = selected_modules), "time_associated_modules_p_lt_0.05.csv", row.names = FALSE)

# -----------------------------
# 6. Prepare TO-GCN input files
# -----------------------------
counts <- read.csv(wgcna_config$count_csv, row.names = 1, check.names = FALSE)
tpm_common <- tpm[rownames(tpm) %in% rownames(counts), , drop = FALSE]
counts_common <- counts[rownames(tpm_common), colnames(tpm_common), drop = FALSE]

if (!all(rownames(counts_common) == rownames(tpm_common))) {
  stop("Counts and TPM matrices must contain the same genes in the same order.")
}

total_counts <- colSums(counts_common)
fpkm <- tpm_common
for (sample in colnames(tpm_common)) {
  fpkm[, sample] <- tpm_common[, sample] * total_counts[sample] / 1e6
}
write.table(fpkm, file = "fpkm_all_genes.tsv", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

annotation <- read.csv(wgcna_config$annotation_csv, check.names = FALSE)
tf_genes <- annotation$GeneID[annotation$type %in% c("Phytohormone;TFgenes", "metabolic_paths;TFgenes", "Redox;TFgenes", "TFgenes")]
tf_genes <- intersect(tf_genes, rownames(fpkm))

fpkm_only_tf <- fpkm[tf_genes, , drop = FALSE]
fpkm_not_tf <- fpkm[setdiff(rownames(fpkm), tf_genes), , drop = FALSE]

write.table(fpkm_only_tf, file = "fpkm_only_tf.tsv", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(fpkm_not_tf, file = "fpkm_not_tf.tsv", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# Shell commands used for TO-GCN, run separately in the TO-GCN installation directory:
# ./Cutoff 5 example_data/fpkm_only_tf.tsv
# ./TO-GCN 5 example_data/fpkm_only_tf.tsv example_data/initial_seed1.txt 0.86
# ./GeneLevel 5 example_data/fpkm_only_tf.tsv example_data/fpkm_not_tf.tsv Node_level.tsv 0.86

sessionInfo()
