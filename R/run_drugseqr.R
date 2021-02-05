#' Run drugseqr app
#'
#' Run drugseqr application to explore single-cell and bulk RNA-seq datasets.
#'
#' @inheritParams init_drugseqr
#' @inheritParams shiny::runApp
#' @param data_dir Directory containing folders \code{'bulk'}, \code{'single-cell'}, and \code{'custom_queries'}.
#'  Ignored if \code{test_data} is \code{TRUE}.
#' @param app_dir Directory containing shiny app files. Can be 'inst/app' if working from source code.
#' @param data_dir Directory containing folder \code{app_name} with saved analyses.
#' @param pert_query_dir Path to directory where pert query results (using CMAP02/L1000 as query signature) will be downloaded as requested.
#' @param pert_signature_dir Path to directory where pert signatures for CMAP02/L1000 will be downloaded as requested.
#' @param indices_dir Path to directory containing \code{kallisto} indices and whitelists.
#' @param tabs Character vector of tabs to include in order desired. Must be subset of 'Single Cell', 'Bulk Data', and 'Drugs'.
#' @param test Boolean indicating if \code{shinytest} should be run (default is \code{FALSE}).
#'  If \code{TRUE} test data will be used.
#' @param test_data Boolean indicating if test data should be used. Default is \code{TRUE}
#'
#' @import rintrojs
#' @import shiny
#' @importFrom shinyjs toggle toggleClass toggleState html addClass removeClass hidden runjs
#' @import org.Hs.eg.db org.Mm.eg.db
#'
#' @return Runs drugseqr app
#' @export
#'
#' @examples
#'
#' # create directory structure for new datasets
#' data_dir <- tempdir()
#' app_name <- 'example'
#' init_drugseqr(app_name, data_dir)
#'
#'
#' # run app
#' # run_drugseqr(app_name, data_dir)
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
                         port = 3838,
                         logout_url = NULL) {

  data_dir <- file.path(data_dir, app_name)
  if (!dir.exists(data_dir)) dir.create(data_dir)

  # pass arguments to app through options then run
<<<<<<< HEAD
  shiny::shinyOptions(data_dir = normalizePath(data_dir),
                      pert_query_dir = normalizePath(pert_query_dir),
                      pert_signature_dir = normalizePath(pert_signature_dir),
                      indices_dir = normalizePath(indices_dir),
                      tabs = tabs)
=======
  shinyOptions(data_dir = normalizePath(data_dir),
               pert_query_dir = normalizePath(pert_query_dir),
               pert_signature_dir = normalizePath(pert_signature_dir),
               indices_dir = normalizePath(indices_dir),
               tabs = tabs,
               logout_url = logout_url)
>>>>>>> master


  if (test) {
    # run test and return
    # shinytest::recordTest(app_dir, seed = 0)
    return(NULL)

  } else if (test_data) {
    # use test data (faster)
    options(shiny.testmode = TRUE)

  } else {
    options(shiny.testmode = FALSE)
  }

  # auto-reload if update app files
  if (grepl('/srv/drugseqr', data_dir)) options(shiny.autoreload = TRUE)

  shiny::runApp(app_dir, launch.browser = TRUE, host = host, port = port)
}


#' Initialize drugseqr folders/files for a new app
#'
#' Creates necessary folders/files for a new drugseqr app inside of /srv/shiny-server/drugseqr.
#'
#' @param app_name Name for new drugseqr app.
#' @param data_dir Path to put \code{app_name} directory where app will be
#' initialized. Default is \code{'/srv/drugseqr'} (for hosting app on server).
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' data_dir <- tempdir()
#' init_drugseqr('example', data_dir)
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
