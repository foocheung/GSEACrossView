suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(fs)
  library(yaml)
})

config <- yaml::read_yaml("merge_gsea.yaml")

clean_cluster_name <- function(x) {
  x %>%
    str_trim() %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("\\+", "") %>%
    str_replace_all("-", "_") %>%
    str_replace_all("__+", "_")
}

extract_geneset <- function(fname, patterns) {
  hits <- names(patterns)[map_lgl(names(patterns), ~ str_detect(fname, .x))]
  if (length(hits) == 0) return("UNKNOWN")
  patterns[[hits[1]]]
}

extract_celltype <- function(f, ds) {
  fname <- basename(f)
  
  if (ds$celltype_from == "parent_after") {
    parts <- str_split(path_norm(f), "/")[[1]]
    idx <- which(parts == ds$parent_anchor)
    
    if (length(idx) > 0 && idx[length(idx)] < length(parts)) {
      return(clean_cluster_name(parts[idx[length(idx)] + 1]))
    }
  }
  
  if (ds$celltype_from == "regex") {
    out <- str_match(fname, ds$celltype_regex)[, 2]
    return(clean_cluster_name(out))
  }
  
  clean_cluster_name(tools::file_path_sans_ext(fname))
}

load_dataset <- function(ds, common_cols, annotation) {
  files <- dir_ls(ds$root, recurse = TRUE, regexp = ds$file_regex)
  
  if (!is.null(ds$exclude_regex)) {
    files <- files[!str_detect(files, ds$exclude_regex)]
  }
  
  map_dfr(files, function(f) {
    df <- fread(f)
    fname <- basename(f)
    
    df$cluster <- extract_celltype(f, ds)
    df$geneset <- extract_geneset(fname, ds$geneset_patterns)
    df$annotation <- annotation
    df$dataset <- ds$name
    
    for (col in common_cols) {
      if (!col %in% colnames(df)) df[[col]] <- NA
    }
    
    df[, common_cols, with = FALSE]
  })
}

common_cols <- unlist(config$columns$common)
annotation <- config$settings$annotation

combined <- map_dfr(
  config$datasets,
  load_dataset,
  common_cols = common_cols,
  annotation = annotation
) %>%
  mutate(
    cluster = clean_cluster_name(cluster),
    pathway_id = paste(ID, Description, sep = "|"),
    celltype_dataset = paste(dataset, cluster, sep = "_")
  )

fwrite(combined, config$project$output_file)

cat("Saved:", config$project$output_file, "\n")
cat("Combined rows:", nrow(combined), "\n")
cat("Datasets:", paste(unique(combined$dataset), collapse = ", "), "\n")