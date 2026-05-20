

https://github.com/user-attachments/assets/ccc35eed-0bba-4316-9554-69622a24b9b1

---

# GSEA CrossView

`GSEACrossView` is a `{golem}`-based Shiny application for exploring and comparing Gene Set Enrichment Analysis (GSEA) results across tissues, cell types, and experimental conditions.

The app is designed for single-cell and pseudobulk transcriptomic workflows where multiple GSEA outputs need to be integrated into a unified exploratory interface.

Typical use cases include:

* PBMC vs ileum comparisons
* disease vs healthy pathway analysis
* trained immunity studies
* cross-condition pathway concordance
* identifying recurring pathways across immune cell populations

---

# Features

## Cross-Tissue Pathway Scatterplots

Compare pathway NES values between tissues:

* concordant upregulated pathways
* concordant downregulated pathways
* discordant tissue-specific pathways
* significance-aware coloring
* faceted by cell type
* interactive hover tooltips

The scatterplots support:

* adjustable plot dimensions
* custom facet layouts
* top pathway annotations
* filtering non-significant pathways

---

## Shared Pathway Heatmaps

Visualize pathways shared across tissues for a selected cell type.

Features:

* significant-only filtering
* top-N pathway selection
* NES-based ranking
* publication-friendly visualization

Rows represent pathways.
Columns represent tissues.
Values represent normalized enrichment scores (NES).

---

## Recurring Pathway Discovery

Identify pathways that recur across multiple cell types within the same directional/significance quadrant.

Useful for detecting:

* conserved inflammatory programs
* interferon signatures
* tissue-specific immune modules
* broadly shared disease pathways

The recurring pathway table summarizes:

* recurrence frequency across cell types
* average NES values
* minimum adjusted p-values
* pathway directionality groups

---

# Input Format

Upload a single combined `.csv` or `.csv.gz` file.

Required columns:

```text
ID
Description
NES
pvalue
p.adjust or p_adjust
cluster
geneset
dataset
```

The app internally standardizes these columns:

| Input column | Standardized column |
| ------------ | ------------------- |
| ID           | pathway_id          |
| Description  | pathway_name        |
| cluster      | cell_type           |
| dataset      | tissue              |
| geneset      | database            |
| NES          | NES                 |
| p_adjust     | p_adj               |
| pvalue       | p_val               |

---

# Example Workflow

## Step 1 — Generate GSEA results

Run differential expression and GSEA separately for each:

* tissue
* cell type
* comparison

Example:

```text
PBMC/CD4_Tcells/GSEA_results.csv
ILEUM/CD4_Tcells/GSEA_results.csv
```

---

## Step 2 — Merge results using YAML-driven pipeline

The repository includes preprocessing scripts:

```text
scripts/
├── merge_gsea_config.yml
├── merge_gsea_from_yaml.R
```

The YAML configuration controls:

* dataset paths
* regex patterns
* tissue names
* geneset mapping
* output file naming

Example:

```yaml
datasets:
  - name: PBMC
    root: data/pbmc/

  - name: ILEUM
    root: data/ileum/
```

Run:

```r
source("scripts/merge_gsea_from_yaml.R")
```

This generates:

```text
combined_gsea.csv
```

---

## Step 3 — Launch the app

```r
devtools::load_all()
run_app()
```

Upload the combined CSV into the app interface.

---

# Repository Structure

```text
GSEACrossView/
├── R/
│   ├── app_ui.R
│   ├── app_server.R
│   ├── utils_gsea.R
│
├── scripts/
│   ├── merge_gsea_config.yml
│   ├── merge_gsea_from_yaml.R
│
├── inst/
├── DESCRIPTION
├── NAMESPACE
├── README.md
└── app.R
```

---

# Installation

## Local installation

```r
install.packages(c(
  "devtools",
  "golem"
))

devtools::install()
```

Run:

```r
GSEACrossView::run_app()
```

---

# Run From Source

```r
pkgload::load_all()
GSEACrossView::run_app()
```

---

# Deployment

## Shiny Server

Create:

```r
app.R
```

containing:

```r
GSEACrossView::run_app()
```

Then deploy normally using Shiny Server.

---

## Posit Connect / shinyapps.io

Deploy directly from the package root using `{rsconnect}` workflows.

---

# Suggested Future Extensions

Potential future modules:

* pathway clustering
* cross-condition delta-NES analysis
* pathway network visualization
* Reactome hierarchy browsing
* longitudinal trajectory GSEA
* pathway overlap UpSet plots
* multi-contrast comparison mode
* ATAC/RNA integrated pathway analysis

---

# Citation

If you use this application in published work, please cite the associated repository and relevant GSEA methodology papers.

---
