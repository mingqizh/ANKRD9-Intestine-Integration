# ANKRD9 Expression Correlation and Pathway Analysis

This repository contains the R scripts used for the analysis of ANKRD9 expression patterns and gene correlation networks in human intestinal tissues, as featured in our manuscript for *Nature Communications*.

## Overview
This pipeline performs the following tasks:
1. **Data Integration**: Converts and loads Tabula Sapiens `.h5ad` objects into Seurat.
2. **Correlation Analysis**: Calculates bi-weight midcorrelation (`WGCNA::bicor`) between ANKRD9 and all other detected genes in specific cell subsets (e.g., mature and immature enterocytes).
3. **GSEA**: Conducts Gene Set Enrichment Analysis using `clusterProfiler::gseGO` to identify functional pathways associated with ANKRD9 expression.

## Data Availability
The human single-nuclear sequencing data (Tabula Sapiens v2) was acquired from FigShare:
[https://figshare.com/articles/dataset/Tabula_Sapiens_v2/27921984](https://figshare.com/articles/dataset/Tabula_Sapiens_v2/27921984)

Bulk expression data was analyzed using GTEx v8 (specifically the small intestine dataset).

## Requirements
To run this script, you will need R (>= 4.0.0) and the following libraries:
- `Seurat` & `SeuratDisk`
- `WGCNA`
- `clusterProfiler`
- `org.Hs.eg.db`
- `ggplot2`
- `enrichplot`

## Usage
The main analysis script is `ANKRD9_analysis.R`. 
- Update the `setwd()` path to your local directory.
- Ensure the `.h5ad` file is present in your data path.
- The function `run_paths_incell(cell_subset)` can be called for any specific cell type annotation found in the 'free_annotation' metadata slot.

## Citation
If you use this code, please cite:
[Insert Full Manuscript Citation once published]
