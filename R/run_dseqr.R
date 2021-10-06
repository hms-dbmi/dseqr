#' Run dseqr app
#'
#' Run dseqr application to explore single-cell and bulk RNA-seq datasets.
#'
#' @param project_name Name of project folder in \code{data_dir}. Will be created if doesn't exist.
#' @param data_dir Directory containing project sub-folders. By default also will contain sub-folders
#'  \code{pert_query_dir}, \code{pert_signature_dir}, and \code{indices_dir}.
#' @param tabs Character vector of tabs to include in order desired. Must be subset of 'Single Cell', 'Bulk Data', and 'Drugs'.
#' @param pert_query_dir Path to directory where pert query results (using CMAP02/L1000 as query signature) will be downloaded as requested.
#' @param pert_signature_dir Path to directory where pert signatures for CMAP02/L1000 will be downloaded as requested.
#' @param indices_dir Path to directory containing \code{kallisto} indices and whitelists.
#' @param app_dir Directory containing shiny app files. Default is to use 'app' directory of
#'  installed dseqr package. Can be 'inst/app' if working from source code.
#' @param test Boolean indicating if \code{shinytest} should be run (default is \code{FALSE}).
#'  If \code{TRUE} test data will be used.
#' @param logout_url URL used to log users out. Default is \code{NULL} for local use.
#' @param is_local Is dseqr running locally? If \code{FALSE}, uses CDNs for
#'  some dependencies and sets up other dseqr.com specific tags. Default is \code{TRUE}
#'  if \code{logout_url} is \code{NULL}, otherwise \code{FALSE}.
#' @inheritParams init_dseqr
#' @inheritParams shiny::runApp
#'
#' @import rintrojs
#' @import shiny
#' @importFrom shinyjs toggle toggleClass toggleState html addClass removeClass hidden runjs
#'
#' @return Runs dseqr app
#' @export
#'
#' @examples
#'
#' if (interactive()) {
#'
#'   data_dir <- tempdir()
#'   project_name <- 'example'
#'   run_dseqr(project_name, data_dir)
#' }
#'
run_dseqr <- function(project_name,
                      data_dir,
                      tabs = c('Single Cell', 'Bulk Data', 'Drugs'),
                      pert_query_dir = file.path(data_dir, '.pert_query_dir'),
                      pert_signature_dir = file.path(data_dir, '.pert_signature_dir'),
                      gs_dir = file.path(data_dir, '.gs_dir'),
                      indices_dir = file.path(data_dir, '.indices_dir'),
                      tx2gene_dir = file.path(data_dir, '.tx2gene_dir'),
                      app_dir = system.file('app', package = 'dseqr', mustWork = TRUE),
                      test = FALSE,
                      host = '0.0.0.0',
                      port = 3838,
                      logout_url = NULL,
                      is_local = is.null(logout_url),
                      is_example = FALSE) {

  if (missing(data_dir) & !is_local) {
    message('Setting data_dir to /srv/dseqr for hosted application.')
    data_dir <- '/srv/dseqr'
  }

  if (missing(data_dir)) stop('data_dir not specified.')

  user_dir <- file.path(data_dir, project_name)
  if (!dir.exists(user_dir)) init_dseqr(project_name, data_dir)

  # pass arguments to app through options then run
  shiny::shinyOptions(
    data_dir = normalizePath(user_dir),
    pert_query_dir = normalizePath(pert_query_dir),
    pert_signature_dir = normalizePath(pert_signature_dir),
    gs_dir = normalizePath(gs_dir),
    indices_dir = normalizePath(indices_dir),
    tx2gene_dir = normalizePath(tx2gene_dir),
    tabs = tabs,
    logout_url = logout_url,
    is_local = is_local,
    is_example = is_example)


  if (test) {
    # TODO: run test and return
    shinytest::recordTest(app_dir, seed = 0, loadTimeout = 100000)
    options(shiny.testmode = TRUE)
    return(NULL)

  } else {
    options(shiny.testmode = FALSE)
  }

  # on remote: send logins errors to slack
  if (!is_local) {
    options(shiny.error = send_slack_error)
    user <- Sys.getenv('SHINYPROXY_USERNAME', 'localhost')
    slack <- readRDS(system.file('extdata/slack.rds', package = 'dseqr'))

    httr::POST(
      url = slack$logins,
      httr::add_headers('Content-Type' = 'application/json'),
      body = sprintf('{"text": "⭐⭐⭐ \n\n *user*: %s 🧑"}', user)
    )
  }

  # allow up to 30GB uploads
  options(shiny.maxRequestSize=30*1024*1024^2)

  # for moving EFS files out of IA
  is_aws <- !is.null(logout_url)
  if (is_aws) Sys.setenv('EFS_LIFECYCLE'=7)

  # auto-reload if update app files
  if (!is_aws) options(shiny.autoreload = TRUE)
  shiny::runApp(app_dir, launch.browser = TRUE, host = host, port = port)
}


#' Initialize dseqr folders/files for a new app
#'
#' Creates necessary folders/files for a new dseqr app inside of /srv/shiny-server/dseqr.
#'
#' @param anal_name Name for new dseqr app.
#' @param data_dir Path to put \code{anal_name} directory where app will be
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
init_dseqr <- function(project_name, data_dir = '/srv/dseqr') {

  user_dir <- file.path(data_dir, project_name)
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
