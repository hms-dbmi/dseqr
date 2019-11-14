#' Run drugseqr app
#'
#' Used for local development
#'
#' @param data_dir Directory containing folders \code{'bulk'}, \code{'single-cell'}, and \code{'custom_queries'}.
#'  Ignored if \code{test_data} is \code{TRUE}.
#' @param app_dir Directory containing drugseqr shiny app files.
#' @param pert_query_dir Path to directory where pert query results (using CMAP02/L1000 as query signature) will be downloaded as requested.
#' @param test Boolean indicating if \code{shinytest} should be run (default is \code{FALSE}).
#'  If \code{TRUE} test data will be used.
#' @param test_data Boolean indicating if test data should be used. Default is \code{TRUE}
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' library(drugseqr)
#'
#' # override default app_dir etc for development
#' app_dir <- 'inst/app'
#' data_dir <- 'data-raw/patient_data/sjia'
#' pert_query_dir <- 'data-raw/drug_gene_queries/data'
#' pert_signature_dir <- 'data-raw/drug_es/signatures'
#'
#' run_drugseqr(data_dir, app_dir, pert_query_dir, pert_signature_dir, test_data = TRUE)
#' run_drugseqr(data_dir, app_dir, pert_query_dir, pert_signature_dir, test_data = FALSE)
#' run_drugseqr(data_dir, app_dir, pert_query_dir, pert_signature_dir, test = TRUE)
#'
run_drugseqr <- function(data_dir,
                         app_dir = system.file('app', package = 'drugseqr', mustWork = TRUE),
                         pert_query_dir = '/srv/drugseqr/pert_query_dir',
                         pert_signature_dir = '/srv/drugseqr/pert_signature_dir',
                         indices_dir = '/srv/drugseqr/indices',
                         test = FALSE,
                         test_data = FALSE,
                         host = '0.0.0.0',
                         port = 3838) {


  # pass arguments to app through options then run
  shiny::shinyOptions(data_dir = normalizePath(data_dir),
                      pert_query_dir = normalizePath(pert_query_dir),
                      pert_signature_dir = normalizePath(pert_signature_dir),
                      indices_dir = normalizePath(indices_dir))


  if (test) {
    # reset data for testing
    data_dir <- 'inst/app/tests/data/test/example'
    static_dir <- 'inst/app/tests/data/static/example'
    unlink(data_dir, recursive = TRUE)
    dir.create(data_dir)
    file.copy(list.files(static_dir, full.names = TRUE), data_dir, recursive = TRUE)

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
  shiny::runApp(app_dir, launch.browser = TRUE, host = host, port = port)
}


#' Initialize drugseqr folders/files for a new app
#'
#' Creates necessary folders/files for a new drugseqr app inside of /srv/shiny-server/drugseqr.
#'
#' @param app_name Name for new drugseqr app.
#' @param app_dir Path to put \code{app_name} directory where app will be initialized. Default is \code{'/srv/drugseqr'} for shiny.
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' # app_dir for local development
#' init_drugseqr('example', app_dir = 'data-raw/patient_data')
#'
init_drugseqr <- function(app_name, app_dir = '/srv/drugseqr') {

  data_dir <- file.path(app_dir, app_name)
  dir.create(data_dir, recursive = TRUE)
  dir.create(file.path(data_dir, 'bulk'))
  dir.create(file.path(data_dir, 'single-cell'))
  dir.create(file.path(data_dir, 'custom_queries'))

  anals <- data.frame(matrix(ncol = 3, nrow = 0), stringsAsFactors = FALSE)
  colnames(anals) <- c("dataset_name", "dataset_dir", "anal_name")
  saveRDS(anals, file.path(data_dir, 'bulk', 'anals.rds'))

  datasets <- data.frame(matrix(ncol = 2, nrow = 0), stringsAsFactors = FALSE)
  colnames(datasets) <- c("dataset_name", "dataset_dir")
  saveRDS(datasets, file.path(data_dir, 'bulk', 'datasets.rds'))

  saveRDS(NULL, file.path(data_dir, 'single-cell', 'integrated.rds'))
}
