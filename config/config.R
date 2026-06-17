# Project-wide configuration for FtFIP37-Nanoplastics-RootAtlas.
# Edit paths here before running the scripts.

project_dir <- normalizePath(getwd(), mustWork = FALSE)

data_dir <- file.path(project_dir, "data")
result_dir <- file.path(project_dir, "results")
metadata_dir <- file.path(project_dir, "metadata")

config <- list(
  # Cell Ranger output directory. Each sample should contain filtered_feature_bc_matrix/.
  cellranger_dir = file.path(data_dir, "cellranger"),
  sample_metadata = file.path(metadata_dir, "sample_metadata.csv"),
  cluster_annotation = file.path(metadata_dir, "cluster_annotation.csv"),
  feature_map = file.path(metadata_dir, "feature_map.csv"),

  # Seurat outputs.
  seurat_outdir = file.path(result_dir, "seurat"),
  seurat_rds = file.path(result_dir, "seurat", "combined_integrated_annotated.rds"),

  # RNA velocity inputs/outputs.
  velocity_outdir = file.path(result_dir, "rna_velocity"),
  loom_metadata = file.path(metadata_dir, "loom_metadata.csv"),

  # Bulk time-course expression and WGCNA inputs/outputs.
  bulk_tpm = file.path(data_dir, "bulk", "transcript_tpm_matrix.csv"),
  bulk_count = file.path(data_dir, "bulk", "gene_count_matrix.csv"),
  bulk_metadata = file.path(metadata_dir, "bulk_sample_metadata.csv"),
  wgcna_outdir = file.path(result_dir, "wgcna"),

  # Main analysis parameters retained from the original workflow.
  min_features = 100,
  max_features = 7500,
  min_counts = 100,
  max_counts = 25000,
  n_variable_features = 2000,
  sample_npcs = 80,
  sample_dims = 50,
  sample_resolution = 1.5,
  integrated_npcs = 50,
  integrated_dims = 25,
  integrated_resolution = 1.0,
  random_seed = 42
)

# Create output folders if they do not exist.
invisible(lapply(c(config$seurat_outdir, config$velocity_outdir, config$wgcna_outdir),
                 dir.create, recursive = TRUE, showWarnings = FALSE))
