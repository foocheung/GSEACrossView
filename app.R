# Convenience app launcher for local testing or Shiny Server.
# In package development, prefer: devtools::load_all(); run_app()

if (requireNamespace('GSEACrossView', quietly = TRUE)) {
  GSEACrossView::run_app()
} else {
  pkgload::load_all()
  run_app()
}
