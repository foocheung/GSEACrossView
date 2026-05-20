#' Run the GSEA CrossView Shiny app
#'
#' @param ... Arguments passed to [shiny::shinyApp()].
#'
#' @export
run_app <- function(...) {
  options(shiny.maxRequestSize = 500 * 1024^2)

  shiny::shinyApp(
    ui = app_ui(),
    server = app_server,
    ...
  )
}
