# Pseudotime Trajectory Analysis of Spatial Transcriptomics Data (GSE283269)

This repository contains the R code, figures, and results for **Task 2: Cell-cell communication/interaction analysis** using spatial transcriptomics data from [GSE283269](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283269). The analysis focuses on identifying dynamic changes in gene expression across pseudotime using the `monocle3` package.

---

## 📁 Repository Structure

```
📦 Task2_results/
 ┣ 📂 plots/                   # All figures generated during the analysis
 ┣ 📂 processed_data/          # CSV outputs of differential genes and GO results
 ┣ 📄 test_run_task2.r         # Full R code used for the analysis
 ┣ 📄 question_2_report.docx   # Final report (ready for manuscript submission)
 ┗ 📄 README.md                # This file
```

---

## 🔬 Dataset

- **Source**: [GSE283269 - NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283269)
- **Platform**: Spatial Transcriptomics
- **Organism**: Mouse
- **Purpose**: Characterize spatial and temporal patterns of gene expression and interactions among different spatial clusters.

---

## 🔧 Analysis Pipeline Overview

1. **Download and preprocess raw count data** (`fasterq-dump`, `hisat2`, `samtools`)
2. **Count matrix generation** via `featureCounts`
3. **Create `cell_data_set` object** with `monocle3`
4. **Preprocessing and dimension reduction** (PCA + UMAP)
5. **Cell clustering** and **trajectory graph construction**
6. **Root node selection** for pseudotime ordering
7. **Differential expression analysis** along pseudotime
8. **Top gene visualization** (line plots, heatmaps)
9. **Gene Ontology enrichment analysis** using `clusterProfiler`
10. **Final reporting** and summary generation

---

## 📊 Key Outputs

- `top_genes_in_pseudotime_filtered.png` – Visualization of top dynamic genes along pseudotime
- `go_enrichment_barplot.png` – Top enriched GO biological processes
- `pseudotime_heatmap_top50.png` – Heatmap of top 50 genes clustered across pseudotime
- `deg_pseudotime_filtered.csv` – Full list of differentially expressed genes (q < 0.01)
- `go_results_filtered.csv` – Filtered GO enrichment results

---

## 📝 How to Reproduce

Open RStudio and run the following:

```r
source("test_run_task2.r")
```

This script will:

- Read and clean count data
- Process metadata
- Run trajectory analysis with Monocle3
- Perform differential expression and GO enrichment
- Generate figures and CSV exports

---

## 📦 Requirements

Install the following R packages:

```r
install.packages("tidyverse")
install.packages("patchwork")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("monocle3", "clusterProfiler", "org.Mm.eg.db"))
```

---

## 💡 Suggested Improvements

- Apply additional dimensionality reduction (e.g., Harmony or batch correction)
- Perform sub-clustering and cell-cell interaction modeling using tools like **CellChat** or **NicheNet**
- Annotate gene clusters using public marker gene databases

---

## 📄 License

This project is provided for academic and educational use only. Please cite the original data source and associated tools used in the analysis.
