#' Runs drugseqr shiny app
#'
#' @return
#' @export
#'
#' @examples
run_drugseqr <- function() {
  appDir <- system.file("shiny-examples", "drugseqr", package = "drugseqr")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `drugseqr`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
