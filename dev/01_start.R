# Development setup for GSEACrossView
# Run this once from the package root when starting development.

install.packages(c(
  'golem', 'shiny', 'config', 'dplyr', 'tidyr', 'ggplot2',
  'tibble', 'DT', 'plotly', 'readr', 'devtools', 'roxygen2', 'testthat'
))

golem::fill_desc(
  pkg_name = 'GSEACrossView',
  pkg_title = 'Cross-Tissue GSEA Explorer',
  pkg_description = 'A golem-based Shiny app for exploring cross-tissue and cross-cell-type GSEA results.',
  author_first_name = 'Foo',
  author_last_name = 'Cheung',
  author_email = 'foo@example.com'
)

devtools::document()
devtools::load_all()
run_app()
