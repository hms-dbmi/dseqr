#' Run drugseqr app
#'
#' Used for local development
#'
#' @param data_dir Directory containing folders \code{'bulk'} and \code{'single-cell'}.
#'  Ignored if \code{test_data} is \code{TRUE}.
#' @param test Boolean indicating if \code{shinytest} should be run (default is \code{FALSE}).
#'  If \code{TRUE} test data will be used.
#' @param test_data Boolean indicating if test data should be used. Default is \code{TRUE}
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' data_dir <- '~/Documents/Batcave/zaklab/drugseqr/data-raw/patient_data/IBD_20190712'
#' data_dir <- '~/Documents/Batcave/zaklab/drugseqr/data-raw/patient_data/IBD'
#' data_dir <- '~/Documents/Batcave/zaklab/drugseqr/data-raw/patient_data/amnon'
#' data_dir <- '~/Documents/Batcave/zaklab/drugseqr/data-raw/patient_data/sjia'
#'
#' run_drugseqr(data_dir, test_data = TRUE)
#' run_drugseqr(data_dir, test_data = FALSE)
#'
run_drugseqr <- function(data_dir, test = FALSE, test_data = TRUE) {

  app_dir <- 'inst/app'

  # pass arguments to app through options then run
  shiny::shinyOptions(data_dir = data_dir)


  if (test) {
    # run test and return
    shinytest::recordTest(app_dir, seed = 0)
    return(NULL)

  } else if (test_data) {
    # use test data (faster)
    options(shiny.testmode = TRUE)

  } else {
    options(shiny.testmode = FALSE)
  }

  # auto-reload if update app files
  options(shiny.autoreload = TRUE)
  shiny::runApp(app_dir, launch.browser = TRUE)

}


#' Initialize drugseqr folders/files for a new app
#'
#' Creates necessary folders/files for a new drugseqr app inside of /srv/shiny-server/drugseqr.
#'
#' @param app_name Name for new drugseqr app.
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' init_drugseqr('example')
#'
init_drugseqr <- function(app_name) {
  app_dir <- file.path('/srv/shiny-server/drugseqr', app)
  data_dir <- file.path(app_dir, 'data_dir')

  dir.create(file.path(data_dir, 'bulk'), recursive = TRUE)
  dir.create(file.path(data_dir, 'single-cell'))
  dir.create(file.path(data_dir, 'custom_queries'))

  anals <- data.frame(matrix(ncol = 3, nrow = 0), stringsAsFactors = FALSE)
  colnames(anals) <- c("dataset_name", "dataset_dir", "anal_name")
  saveRDS(anals, file.path(data_dir, 'bulk', 'anals.rds'))

  saveRDS(NULL, file.path(data_dir, 'bulk', 'datasets.rds'))
  saveRDS(NULL, file.path(data_dir, 'single-cell', 'integrated.rds'))
}
