#' Access files in the current app
#'
#' @param ... Character vectors, specifying subdirectory and file(s) within
#'   `inst/app/www`.
#'
#' @noRd
app_sys <- function(...) {
  system.file('app/www', ..., package = 'GSEACrossView')
}
