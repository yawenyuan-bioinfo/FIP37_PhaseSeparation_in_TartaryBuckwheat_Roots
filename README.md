# FtFIP37-Nanoplastics-RootAtlas

Code for: **FIP37 phase separation facilitates nanoplastics sensing and adaptive responses in tartary buckwheat root systems**.

This repository contains cleaned and annotated analysis scripts for tartary buckwheat root single-nucleus RNA-seq, RNA velocity, pseudotime inference, candidate gene visualization, and co-expression network analyses. The code was reorganized from project analysis notes into a GitHub-ready structure with centralized configuration and metadata templates.

## Repository structure

```text
FtFIP37-Nanoplastics-RootAtlas/
├── config/
│   └── config.R                         # Central path and parameter settings
├── metadata/
│   ├── sample_metadata.csv              # snRNA-seq sample information template
│   ├── cluster_annotation.csv           # Manual cluster-to-cell-type annotation table
│   ├── bulk_sample_metadata.csv         # Bulk time-course metadata template
│   └── loom_metadata.csv                # Velocyto loom file metadata template
├── env/
│   ├── R_packages.txt                   # Required R packages
│   └── requirements.txt                 # Required Python packages
├── scripts/
│   ├── 01_seurat_qc_integration_annotation.R
│   ├── 02_export_seurat_to_anndata.R
│   ├── 03_scvelo_cellrank_velocity.py
│   ├── 04_marker_module_and_gene_plots.R
│   ├── 05_bulk_timecourse_wgcna_togcn.R
│   ├── 06_hdwgcna_single_cell_network.R
│   ├── 07_replicate_qc_and_monocle_trajectory.R
│   └── 08_aquaporin_and_candidate_heatmaps.R
└── results/                             # Generated analysis outputs
```

## Analysis workflow

1. **snRNA-seq preprocessing and annotation**  
   `scripts/01_seurat_qc_integration_annotation.R` reads Cell Ranger matrices, performs QC filtering, normalization, PCA, clustering, DoubletFinder-based doublet removal, Seurat integration, marker gene detection, and marker-based root cell annotation.

2. **Export to AnnData**  
   `scripts/02_export_seurat_to_anndata.R` exports the annotated Seurat object to `.h5ad` for Python-based RNA velocity analysis.

3. **RNA velocity and CellRank**  
   `scripts/03_scvelo_cellrank_velocity.py` merges Seurat-derived AnnData with velocyto loom files and performs scVelo RNA velocity, CellRank transition projection, and diffusion pseudotime analysis.

4. **Marker and candidate gene visualization**  
   `scripts/04_marker_module_and_gene_plots.R` generates marker module scores, dot plots, violin plots, heatmaps, and FIP37-related gene-expression plots.

5. **Bulk time-course WGCNA and TO-GCN preparation**  
   `scripts/05_bulk_timecourse_wgcna_togcn.R` constructs co-expression modules from bulk time-course expression matrices and exports files for TO-GCN analysis.

6. **Cell-type-specific hdWGCNA**  
   `scripts/06_hdwgcna_single_cell_network.R` builds metacell-based co-expression networks for selected root cell types and exports Cytoscape-compatible network files.

7. **Replicate QC and Monocle pseudotime**  
   `scripts/07_replicate_qc_and_monocle_trajectory.R` performs replicate expression-correlation checks and cell-type-wise Monocle trajectory inference.

8. **Aquaporin and candidate heatmaps**  
   `scripts/08_aquaporin_and_candidate_heatmaps.R` summarizes FIP37-related genes and aquaporin-family candidates across conditions and annotated cell types.

## How to run

Edit `config/config.R` and the metadata templates in `metadata/` before running the scripts. From the repository root, run scripts in numerical order:

```bash
Rscript scripts/01_seurat_qc_integration_annotation.R
Rscript scripts/02_export_seurat_to_anndata.R
python scripts/03_scvelo_cellrank_velocity.py
Rscript scripts/04_marker_module_and_gene_plots.R
Rscript scripts/05_bulk_timecourse_wgcna_togcn.R
Rscript scripts/06_hdwgcna_single_cell_network.R
Rscript scripts/07_replicate_qc_and_monocle_trajectory.R
Rscript scripts/08_aquaporin_and_candidate_heatmaps.R
```

## Notes for reproducibility

- Absolute local paths from the original analysis were removed and replaced with `config/config.R`.
- Sample grouping is controlled by `metadata/sample_metadata.csv` rather than inferred from barcode suffixes.
- Manual cell-type annotation is stored in `metadata/cluster_annotation.csv` and can be updated after marker inspection.
- The scripts retain the main thresholds and parameters used in the original analysis, including QC filters, PCA dimensions, clustering resolution, DoubletFinder settings, WGCNA parameters, and candidate gene lists.
- The repository does not include raw sequencing data. Users should place Cell Ranger outputs, loom files, and bulk expression matrices under `data/` or update the paths in `config/config.R`.

## Short description

Code for snRNA-seq QC, integration, clustering, marker-based root cell annotation, FIP37-related gene visualization, pseudotime/RNA velocity analysis, and WGCNA/hdWGCNA/TO-GCN co-expression analysis in tartary buckwheat roots under nanoplastic treatment.
