#' App server
#'
#' @noRd
app_server <- function(input, output, session) {

  gsea_raw <- shiny::reactive({
    shiny::req(input$gsea_file)

    utils::read.csv(
      input$gsea_file$datapath,
      check.names = FALSE
    )
  })

  gsea_clean <- shiny::reactive({
    standardize_gsea(gsea_raw())
  })

  shiny::observeEvent(gsea_clean(), {
    df <- gsea_clean()

    shiny::updateSelectInput(
      session,
      'cell_type',
      choices = c('All', sort(unique(df$cell_type))),
      selected = 'All'
    )

    shiny::updateSelectInput(
      session,
      'heatmap_cell_type',
      choices = sort(unique(df$cell_type)),
      selected = sort(unique(df$cell_type))[1]
    )

    shiny::updateSelectInput(
      session,
      'database',
      choices = c('All', sort(unique(df$database))),
      selected = 'All'
    )

    shiny::updateSelectizeInput(
      session,
      'scatter_celltypes',
      choices = sort(unique(df$cell_type)),
      selected = head(sort(unique(df$cell_type)), 4),
      server = TRUE
    )
  })

  filtered_gsea <- shiny::reactive({
    df <- gsea_clean()

    if (!is.null(input$cell_type) && input$cell_type != 'All') {
      df <- df |> dplyr::filter(cell_type == input$cell_type)
    }

    if (!is.null(input$database) && input$database != 'All') {
      df <- df |> dplyr::filter(database == input$database)
    }

    df
  })

  cross_tissue_df <- shiny::reactive({
    make_cross_tissue_table(
      df0 = gsea_clean(),
      database_filter = input$database,
      scatter_celltypes = input$scatter_celltypes,
      padj_cutoff = input$padj_cutoff,
      remove_neither_sig = input$remove_neither_sig
    )
  })

  output$gsea_table <- DT::renderDT({
    filtered_gsea() |>
      dplyr::select(
        tissue, cell_type, database,
        pathway_id, pathway_name,
        NES, p_adj, p_val,
        sig, direction,
        dplyr::everything()
      ) |>
      DT::datatable(
        filter = 'top',
        options = list(pageLength = 20, scrollX = TRUE)
      )
  })

  output$heatmap <- plotly::renderPlotly({
    df0 <- gsea_clean()

    if (!is.null(input$database) && input$database != 'All') {
      df0 <- df0 |> dplyr::filter(database == input$database)
    }

    shiny::req(input$heatmap_cell_type)

    df0 <- df0 |>
      dplyr::filter(cell_type == input$heatmap_cell_type)

    tissue_names <- sort(unique(df0$tissue))

    shiny::validate(
      shiny::need(length(tissue_names) >= 2, 'Need at least two tissues for shared pathway heatmap.')
    )

    wide <- df0 |>
      dplyr::select(pathway_id, pathway_name, tissue, NES, p_adj) |>
      dplyr::distinct() |>
      tidyr::pivot_wider(
        names_from = tissue,
        values_from = c(NES, p_adj)
      )

    nes_cols <- paste0('NES_', tissue_names)
    padj_cols <- paste0('p_adj_', tissue_names)

    wide <- wide |>
      dplyr::filter(dplyr::if_all(dplyr::all_of(nes_cols), ~ !is.na(.x)))

    if (isTRUE(input$heatmap_shared_only)) {
      wide <- wide |>
        dplyr::filter(dplyr::if_all(dplyr::all_of(padj_cols), ~ .x <= input$padj_cutoff))
    } else {
      wide <- wide |>
        dplyr::filter(dplyr::if_any(dplyr::all_of(padj_cols), ~ .x <= input$padj_cutoff))
    }

    wide <- wide |>
      dplyr::mutate(
        mean_abs_NES = rowMeans(abs(dplyr::across(dplyr::all_of(nes_cols))), na.rm = TRUE)
      ) |>
      dplyr::arrange(dplyr::desc(mean_abs_NES)) |>
      dplyr::slice_head(n = input$heatmap_top_n)

    shiny::validate(
      shiny::need(nrow(wide) > 0, 'No shared pathways found for this cell type and cutoff.')
    )

    mat <- wide |>
      dplyr::select(pathway_name, dplyr::all_of(nes_cols)) |>
      tibble::column_to_rownames('pathway_name') |>
      as.matrix()

    colnames(mat) <- gsub('^NES_', '', colnames(mat))

    plotly::plot_ly(
      x = colnames(mat),
      y = rownames(mat),
      z = mat,
      type = 'heatmap',
      colorscale = 'RdBu',
      reversescale = TRUE
    ) |>
      plotly::layout(
        width = input$heatmap_width,
        height = input$heatmap_height,
        title = paste0('Shared Pathway NES Heatmap: ', input$heatmap_cell_type),
        xaxis = list(title = 'Tissue'),
        yaxis = list(title = 'Pathway')
      )
  })

  output$scatter <- plotly::renderPlotly({
    df <- cross_tissue_df()

    shiny::validate(shiny::need(!is.null(df), 'Need at least two tissues.'))
    shiny::validate(shiny::need(nrow(df) > 0, 'No scatter points remain after filtering.'))

    tissue_names <- attr(df, 'tissue_names')
    x_col <- attr(df, 'x_col')
    y_col <- attr(df, 'y_col')
    x_padj_col <- attr(df, 'x_padj_col')
    y_padj_col <- attr(df, 'y_padj_col')

    label_df <- df |>
      dplyr::group_by(cell_type, label_group) |>
      dplyr::slice_max(
        order_by = label_score,
        n = input$top_labels,
        with_ties = FALSE
      ) |>
      dplyr::ungroup()

    facet_nrow <- input$facet_rows
    if (is.null(facet_nrow) || is.na(facet_nrow)) {
      facet_nrow <- NULL
    }

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data[[x_col]],
        y = .data[[y_col]],
        color = sig_status,
        shape = direction_class,
        text = paste0(
          pathway_name,
          '<br>Cell Type: ', cell_type,
          '<br>Direction: ', direction_class,
          '<br>Significance: ', sig_status,
          '<br>', tissue_names[1], ' NES: ', round(.data[[x_col]], 2),
          '<br>', tissue_names[1], ' p.adj: ', signif(.data[[x_padj_col]], 3),
          '<br>', tissue_names[2], ' NES: ', round(.data[[y_col]], 2),
          '<br>', tissue_names[2], ' p.adj: ', signif(.data[[y_padj_col]], 3)
        )
      )
    ) +
      ggplot2::geom_point(alpha = 0.75, size = input$scatter_point_size) +
      ggplot2::geom_vline(xintercept = 0, linetype = 'dashed') +
      ggplot2::geom_hline(yintercept = 0, linetype = 'dashed') +
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(label = pathway_name),
        size = input$scatter_label_size,
        show.legend = FALSE,
        vjust = -0.6
      ) +
      ggplot2::facet_wrap(
        ~cell_type,
        scales = 'free',
        ncol = input$facet_cols,
        nrow = facet_nrow
      ) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x = paste0(tissue_names[1], ' NES'),
        y = paste0(tissue_names[2], ' NES'),
        color = 'Significance',
        shape = 'Direction',
        title = 'Cross-Tissue GSEA Comparison'
      )

    plotly::ggplotly(
      p,
      tooltip = 'text',
      width = input$scatter_width,
      height = input$scatter_height
    )
  })

  output$recurring_table <- DT::renderDT({
    df <- cross_tissue_df()

    shiny::validate(shiny::need(!is.null(df), 'Need at least two tissues.'))
    shiny::validate(shiny::need(nrow(df) > 0, 'No recurring pathways remain after filtering.'))

    tissue_names <- attr(df, 'tissue_names')
    x_col <- attr(df, 'x_col')
    y_col <- attr(df, 'y_col')
    x_padj_col <- attr(df, 'x_padj_col')
    y_padj_col <- attr(df, 'y_padj_col')

    recurring <- df |>
      dplyr::group_by(quadrant_group, pathway_id, pathway_name) |>
      dplyr::summarise(
        recurring_cell_types = dplyr::n_distinct(cell_type),
        cell_types = paste(sort(unique(cell_type)), collapse = ', '),
        mean_abs_NES = mean(mean_abs_NES, na.rm = TRUE),
        mean_NES_tissue_1 = mean(.data[[x_col]], na.rm = TRUE),
        mean_NES_tissue_2 = mean(.data[[y_col]], na.rm = TRUE),
        min_padj_tissue_1 = min(.data[[x_padj_col]], na.rm = TRUE),
        min_padj_tissue_2 = min(.data[[y_padj_col]], na.rm = TRUE),
        .groups = 'drop'
      ) |>
      dplyr::mutate(
        tissue_1 = tissue_names[1],
        tissue_2 = tissue_names[2]
      ) |>
      dplyr::group_by(quadrant_group) |>
      dplyr::arrange(dplyr::desc(recurring_cell_types), dplyr::desc(mean_abs_NES), .by_group = TRUE) |>
      dplyr::slice_head(n = input$recurring_top_n) |>
      dplyr::ungroup() |>
      dplyr::select(
        quadrant_group,
        pathway_name,
        pathway_id,
        recurring_cell_types,
        cell_types,
        tissue_1,
        mean_NES_tissue_1,
        min_padj_tissue_1,
        tissue_2,
        mean_NES_tissue_2,
        min_padj_tissue_2,
        mean_abs_NES
      )

    DT::datatable(
      recurring,
      filter = 'top',
      options = list(
        pageLength = 20,
        scrollX = TRUE
      )
    )
  })

  output$summary_text <- shiny::renderPrint({
    df <- filtered_gsea()

    cat('Total pathway results:', nrow(df), '\n')
    cat('Tissues:', paste(sort(unique(df$tissue)), collapse = ', '), '\n')
    cat('Cell types:', paste(sort(unique(df$cell_type)), collapse = ', '), '\n')
    cat('Databases:', paste(sort(unique(df$database)), collapse = ', '), '\n')
    cat('Significant pathways:', sum(df$p_adj <= input$padj_cutoff, na.rm = TRUE), '\n')
  })

  output$download_clean <- shiny::downloadHandler(
    filename = function() {
      'standardized_gsea_results.csv'
    },
    content = function(file) {
      utils::write.csv(gsea_clean(), file, row.names = FALSE)
    }
  )
}
