#' Run dseqr app
#'
#' Run dseqr application to explore single-cell and bulk RNA-seq datasets.
#'
#' @inheritParams init_dseqr
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
#' @return Runs dseqr app
#' @export
#'
#' @examples
#'
#' # create directory structure for new datasets
#' data_dir <- tempdir()
#' app_name <- 'example'
#' init_dseqr(app_name, data_dir)
#'
#'
#' # run app
#' # run_dseqr(app_name, data_dir)
#'
run_dseqr <- function(app_name,
                      data_dir = '/srv/dseqr',
                      app_dir = system.file('app', package = 'dseqr', mustWork = TRUE),
                      pert_query_dir = file.path(data_dir, 'pert_query_dir'),
                      pert_signature_dir = file.path(data_dir, 'pert_signature_dir'),
                      gs_dir = file.path(data_dir, 'gs_dir'),
                      indices_dir = file.path(data_dir, 'indices'),
                      tx2gene_dir = file.path(data_dir, 'tx2gene'),
                      tabs = c('Single Cell', 'Bulk Data', 'Drugs'),
                      test = FALSE,
                      test_data = FALSE,
                      host = '0.0.0.0',
                      port = 3838,
                      logout_url = NULL,
                      is_example = FALSE) {

  user_dir <- file.path(data_dir, app_name)
  if (!dir.exists(user_dir)) init_dseqr(app_name, data_dir)

  # pass arguments to app through options then run
  shinyOptions(data_dir = normalizePath(user_dir),
               pert_query_dir = normalizePath(pert_query_dir),
               pert_signature_dir = normalizePath(pert_signature_dir),
               gs_dir = normalizePath(gs_dir),
               indices_dir = normalizePath(indices_dir),
               tx2gene_dir = normalizePath(tx2gene_dir),
               tabs = tabs,
               logout_url = logout_url,
               is_example = is_example)


  if (test) {
    # TODO: run test and return
    shinytest::recordTest(app_dir, seed = 0, loadTimeout = 100000)
    options(shiny.testmode = TRUE)
    return(NULL)

  } else {
    options(shiny.testmode = FALSE)
  }

  # allow up to 30GB uploads
  options(shiny.maxRequestSize=30*1024*1024^2)

  # for moving EFS files out of IA
  is_aws <- !is.null(logout_url)
  if (is_aws) Sys.setenv('EFS_LIFECYCLE'=7)

  # auto-reload if update app files
  if (!is_aws) options(shiny.autoreload = TRUE)
  runApp(app_dir, launch.browser = TRUE, host = host, port = port)
}


#' Initialize dseqr folders/files for a new app
#'
#' Creates necessary folders/files for a new dseqr app inside of /srv/shiny-server/dseqr.
#'
#' @param app_name Name for new dseqr app.
#' @param data_dir Path to put \code{app_name} directory where app will be
#' initialized. Default is \code{'/srv/dseqr'} (for hosting app on server).
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' data_dir <- tempdir()
#' init_dseqr('example', data_dir)
#'
init_dseqr <- function(app_name, data_dir = '/srv/dseqr') {

  user_dir <- file.path(data_dir, app_name)
  dir.create(user_dir, recursive = TRUE)
  dir.create(file.path(user_dir, 'bulk'))
  dir.create(file.path(user_dir, 'single-cell'))
  dir.create(file.path(user_dir, 'custom_queries'))

  anals <- data.frame(matrix(ncol = 3, nrow = 0), stringsAsFactors = FALSE)
  colnames(anals) <- c("dataset_name", "dataset_dir", "anal_name")
  qs::qsave(anals, file.path(user_dir, 'bulk', 'anals.qs'))

  datasets <- data.frame(matrix(ncol = 2, nrow = 0), stringsAsFactors = FALSE)
  colnames(datasets) <- c("dataset_name", "dataset_dir")
  qs::qsave(datasets, file.path(user_dir, 'bulk', 'datasets.qs'))
}
