# GSEA CrossView

`GSEACrossView` is a `{golem}` Shiny app for exploring cross-tissue and cross-cell-type GSEA results.

## Input file

Upload one combined `.csv` or `.csv.gz` GSEA table with these columns:

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

The app standardizes these internally:

| Input column | App column |
|---|---|
| ID | pathway_id |
| Description | pathway_name |
| cluster | cell_type |
| dataset | tissue |
| geneset | database |
| NES | NES |
| p_adjust | p_adj |
| pvalue | p_val |

## Main tabs

- **Table**: searchable standardized GSEA results.
- **Shared Pathway Heatmap**: one selected cell type at a time; rows are shared pathways, columns are tissues, values are NES.
- **Tissue Comparison**: scatterplot comparing NES between the first two tissues, faceted by selected cell types.
- **Recurring Pathways**: top recurring pathways per quadrant/significance group across selected cell types.
- **Summary**: basic data summary.

## Development

```r
install.packages(c('devtools', 'golem'))
devtools::load_all()
run_app()
```

## Install locally

```r
devtools::install()
GSEACrossView::run_app()
```

## Run from source

```r
pkgload::load_all()
GSEACrossView::run_app()
```

## Deploy

For Shiny Server, install the package on the server and create an `app.R` containing:

```r
GSEACrossView::run_app()
```

For Posit Connect or shinyapps.io, deploy from the package root using standard `{rsconnect}` workflows.
