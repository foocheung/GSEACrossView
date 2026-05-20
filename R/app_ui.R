#' App UI
#'
#' @noRd
app_ui <- function(request) {
  shiny::fluidPage(
    shiny::titlePanel('GSEA CrossView'),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::fileInput(
          'gsea_file',
          'Upload combined GSEA CSV',
          accept = c('.csv', '.csv.gz')
        ),

        shiny::hr(),
        shiny::h4('Global Filters'),

        shiny::selectInput('cell_type', 'Table Cell Type', choices = NULL),
        shiny::selectInput('database', 'Gene Set Database', choices = NULL),

        shiny::numericInput(
          'padj_cutoff',
          'Adjusted p-value cutoff',
          value = 0.05,
          min = 0,
          max = 1,
          step = 0.01
        ),

        shiny::hr(),
        shiny::h4('Shared Pathway Heatmap Options'),

        shiny::selectInput(
          'heatmap_cell_type',
          'Heatmap Cell Type',
          choices = NULL
        ),

        shiny::checkboxInput(
          'heatmap_shared_only',
          'Only pathways significant in both tissues',
          value = TRUE
        ),

        shiny::numericInput(
          'heatmap_top_n',
          'Top shared pathways',
          value = 50,
          min = 5,
          max = 500,
          step = 5
        ),

        shiny::numericInput('heatmap_height', 'Heatmap Height', value = 800, min = 300, max = 3000, step = 50),
        shiny::numericInput('heatmap_width', 'Heatmap Width', value = 900, min = 300, max = 3000, step = 50),

        shiny::hr(),
        shiny::h4('Scatter Plot Options'),

        shiny::selectizeInput(
          'scatter_celltypes',
          'Scatter Cell Types',
          choices = NULL,
          multiple = TRUE
        ),

        shiny::numericInput(
          'top_labels',
          'Top labels per group',
          value = 5,
          min = 0,
          max = 100,
          step = 1
        ),

        shiny::checkboxInput(
          'remove_neither_sig',
          'Remove Neither Significant',
          value = TRUE
        ),

        shiny::numericInput(
          'facet_cols',
          'Facet columns',
          value = 2,
          min = 1,
          max = 10,
          step = 1
        ),

        shiny::numericInput(
          'facet_rows',
          'Facet rows',
          value = NA,
          min = 1,
          max = 10,
          step = 1
        ),

        shiny::numericInput(
          'scatter_height',
          'Scatter Plot Height',
          value = 900,
          min = 300,
          max = 3000,
          step = 50
        ),

        shiny::numericInput(
          'scatter_width',
          'Scatter Plot Width',
          value = 1400,
          min = 300,
          max = 3000,
          step = 50
        ),

        shiny::numericInput(
          'scatter_point_size',
          'Scatter Point Size',
          value = 2,
          min = 0.5,
          max = 10,
          step = 0.5
        ),

        shiny::numericInput(
          'scatter_label_size',
          'Scatter Label Size',
          value = 3,
          min = 1,
          max = 10,
          step = 0.5
        ),

        shiny::hr(),
        shiny::h4('Recurring Pathway Options'),

        shiny::numericInput(
          'recurring_top_n',
          'Top recurring pathways per quadrant',
          value = 10,
          min = 1,
          max = 100,
          step = 1
        ),

        shiny::hr(),

        shiny::downloadButton(
          'download_clean',
          'Download Standardized GSEA'
        )
      ),

      shiny::mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel('Table', DT::DTOutput('gsea_table')),
          shiny::tabPanel('Shared Pathway Heatmap', plotly::plotlyOutput('heatmap', height = '100%', width = '100%')),
          shiny::tabPanel('Tissue Comparison', plotly::plotlyOutput('scatter', height = '100%', width = '100%')),
          shiny::tabPanel('Recurring Pathways', DT::DTOutput('recurring_table')),
          shiny::tabPanel('Summary', shiny::verbatimTextOutput('summary_text'))
        )
      )
    )
  )
}
