#' Standardize uploaded GSEA table
#'
#' @param df A data frame from a combined GSEA CSV.
#'
#' @return A standardized data frame with pathway, cell type, tissue, database,
#' NES, adjusted p-value, and direction columns.
#' @noRd
standardize_gsea <- function(df) {
  names(df) <- trimws(names(df))
  names(df) <- gsub('\\s+', '_', names(df))
  names(df) <- gsub('\\.+', '_', names(df))

  required_cols <- c(
    'ID', 'Description', 'NES', 'pvalue', 'p_adjust',
    'cluster', 'geneset', 'dataset'
  )

  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    stop(
      paste0(
        'Missing columns: ', paste(missing_cols, collapse = ', '),
        '\nAvailable columns are: ', paste(names(df), collapse = ', ')
      ),
      call. = FALSE
    )
  }

  df |>
    dplyr::rename(
      pathway_name = Description,
      cell_type = cluster,
      tissue = dataset,
      database = geneset,
      p_adj = p_adjust,
      p_val = pvalue
    ) |>
    dplyr::mutate(
      pathway_id = ID,
      NES = as.numeric(NES),
      p_adj = as.numeric(p_adj),
      p_val = as.numeric(p_val),
      sig = p_adj < 0.05,
      direction = dplyr::case_when(
        NES > 0 ~ 'Up',
        NES < 0 ~ 'Down',
        TRUE ~ 'Neutral'
      )
    )
}

#' Build cross-tissue comparison table
#'
#' @noRd
make_cross_tissue_table <- function(df0, database_filter, scatter_celltypes,
                                    padj_cutoff, remove_neither_sig) {
  if (!is.null(database_filter) && database_filter != 'All') {
    df0 <- df0 |> dplyr::filter(database == database_filter)
  }

  if (!is.null(scatter_celltypes) && length(scatter_celltypes) > 0) {
    df0 <- df0 |> dplyr::filter(cell_type %in% scatter_celltypes)
  }

  tissue_names <- sort(unique(df0$tissue))

  if (length(tissue_names) < 2) {
    return(NULL)
  }

  df <- df0 |>
    dplyr::select(pathway_id, pathway_name, cell_type, tissue, NES, p_adj) |>
    dplyr::distinct() |>
    tidyr::pivot_wider(
      names_from = tissue,
      values_from = c(NES, p_adj)
    )

  x_col <- paste0('NES_', tissue_names[1])
  y_col <- paste0('NES_', tissue_names[2])
  x_padj_col <- paste0('p_adj_', tissue_names[1])
  y_padj_col <- paste0('p_adj_', tissue_names[2])

  df <- df |>
    dplyr::filter(
      !is.na(.data[[x_col]]),
      !is.na(.data[[y_col]])
    ) |>
    dplyr::mutate(
      direction_class = dplyr::case_when(
        .data[[x_col]] > 0 & .data[[y_col]] > 0 ~ 'Concordant Up',
        .data[[x_col]] < 0 & .data[[y_col]] < 0 ~ 'Concordant Down',
        TRUE ~ 'Discordant'
      ),
      sig_status = dplyr::case_when(
        .data[[x_padj_col]] <= padj_cutoff &
          .data[[y_padj_col]] <= padj_cutoff ~ 'Both significant',
        .data[[x_padj_col]] <= padj_cutoff &
          .data[[y_padj_col]] > padj_cutoff ~ paste0(tissue_names[1], ' significant only'),
        .data[[x_padj_col]] > padj_cutoff &
          .data[[y_padj_col]] <= padj_cutoff ~ paste0(tissue_names[2], ' significant only'),
        TRUE ~ 'Neither significant'
      ),
      quadrant_group = paste(direction_class, sig_status, sep = ' | '),
      label_group = quadrant_group,
      label_score = abs(.data[[x_col]]) + abs(.data[[y_col]]),
      mean_abs_NES = label_score / 2
    )

  if (isTRUE(remove_neither_sig)) {
    df <- df |> dplyr::filter(sig_status != 'Neither significant')
  }

  attr(df, 'tissue_names') <- tissue_names
  attr(df, 'x_col') <- x_col
  attr(df, 'y_col') <- y_col
  attr(df, 'x_padj_col') <- x_padj_col
  attr(df, 'y_padj_col') <- y_padj_col

  df
}
