#' Run drugseqr app
#'
#' Used for local development
#'
#' @param data_dir Directory containing folders \code{'bulk'}, \code{'single-cell'}, and \code{'custom_queries'}.
#'  Ignored if \code{test_data} is \code{TRUE}.
#' @param data_dir Directory containing drugseqr shiny app files.
#' @param pert_query_dir Path to directory where pert query results (using CMAP02/L1000 as query signature) will be downloaded as requested.
#' @param pert_signature_dir Path to directory where pert signatures for CMAP02/L1000 will be downloaded as requested.
#' @param indices_dir Path to directory containing \code{kallisto} indices and whitelists.
#' @param test Boolean indicating if \code{shinytest} should be run (default is \code{FALSE}).
#'  If \code{TRUE} test data will be used.
#' @param test_data Boolean indicating if test data should be used. Default is \code{TRUE}
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' app_name <- 'example'
#' data_dir <- 'data_dir'
#' drugseqr::run_drugseqr(app_name, data_dir)
#'
#'
#' # override default data_dir etc for development
#' app_dir <- 'inst/app'
#' data_dir <- 'data-raw/patient_data'
#' pert_query_dir <- 'data-raw/drug_gene_queries/data'
#' pert_signature_dir <- 'data-raw/drug_es/signatures'
#' indices_dir <- '/srv/drugseqr/indices'
#'
#' drugseqr::run_drugseqr('covid19', data_dir, app_dir, pert_query_dir, pert_signature_dir, indices_dir, port = 3839)
#'
run_drugseqr <- function(app_name,
                         data_dir = '/srv/drugseqr',
                         app_dir = system.file('app', package = 'drugseqr', mustWork = TRUE),
                         pert_query_dir = file.path(data_dir, 'pert_query_dir'),
                         pert_signature_dir = file.path(data_dir, 'pert_signature_dir'),
                         indices_dir = file.path(data_dir, 'indices'),
                         tabs = c('Single Cell', 'Bulk Data', 'Drugs'),
                         test = FALSE,
                         test_data = FALSE,
                         host = '0.0.0.0',
                         port = 3838) {

  data_dir <- file.path(data_dir, app_name)
  if (!dir.exists(data_dir)) dir.create(data_dir)

  # pass arguments to app through options then run
  shiny::shinyOptions(data_dir = normalizePath(data_dir),
                      pert_query_dir = normalizePath(pert_query_dir),
                      pert_signature_dir = normalizePath(pert_signature_dir),
                      indices_dir = normalizePath(indices_dir),
                      tabs = tabs)


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
#' init_drugseqr('example', data_dir = 'data-raw/patient_data')
#'
init_drugseqr <- function(app_name, data_dir = '/srv/drugseqr') {

  data_dir <- file.path(data_dir, app_name)
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
