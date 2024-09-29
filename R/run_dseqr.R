#' Run dseqr app
#'
#' Run dseqr application to explore single-cell and bulk RNA-seq datasets.
#'
#' @param app_name Name of user folder in \code{data_dir}. Will be created if doesn't exist.
#' @param data_dir Directory containing app folders. By default also will contain folders
#'  \code{.pert_query_dir}, \code{.pert_signature_dir}, and \code{.indices_dir}.
#' @param tabs Character vector of tabs to include in order desired. Must be subset of 'Single Cell', 'Bulk Data', and 'Drugs'.
#' @param pert_query_dir Path to directory where pert query results (using CMAP02/L1000 as query signature) will be downloaded as requested.
#' @param pert_signature_dir Path to directory where pert signatures for CMAP02/L1000 will be downloaded as requested.
#' @param gs_dir Path to directory where gene to Gene Ontology maps produced by \link{get_genego} are saved.
#' @param indices_dir Path to directory containing \code{kallisto} indices and whitelists.
#' @param tx2gene_dir Path to directory containing transcript to gene maps
#'  produced by \link[dseqr.data]{load_tx2gene}.
#' @param app_dir Directory containing shiny app files. Default is to use 'app' directory of
#'  installed dseqr package. Can be 'inst/app' if working from source code.
#' @param logout_url URL used to log users out. Default is \code{NULL} for local use.
#' @param is_local Is dseqr running locally? If \code{FALSE}, uses CDNs for
#'  some dependencies and sets up other dseqr.com specific tags. Default is \code{TRUE}
#'  if \code{logout_url} is \code{NULL}, otherwise \code{FALSE}.
#' @param test Create a test with shinytest2? Default is \code{FALSE}.
#' @param is_example Is the app for demonstration purposes? If \code{TRUE}, various features
#'  are disabled to prevent dataset changes.
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
#'   app_name <- 'example'
#'   run_dseqr(app_name, data_dir)
#' }
#'
run_dseqr <- function(app_name,
                      data_dir,
                      tabs = c('Single Cell', 'Bulk Data', 'Drugs'),
                      pert_query_dir = file.path(data_dir, '.pert_query_dir'),
                      pert_signature_dir = file.path(data_dir, '.pert_signature_dir'),
                      gs_dir = file.path(data_dir, '.gs_dir'),
                      indices_dir = file.path(data_dir, '.indices_dir'),
                      tx2gene_dir = file.path(data_dir, '.tx2gene_dir'),
                      app_dir = system.file('app', package = 'dseqr', mustWork = TRUE),
                      host = '0.0.0.0',
                      port = 3838,
                      logout_url = NULL,
                      is_local = is.null(logout_url),
                      test = FALSE,
                      is_example = FALSE) {

  if (missing(data_dir) & !is_local) {
    message('Setting data_dir to /srv/dseqr for hosted application.')
    data_dir <- '/srv/dseqr'
  }

  # gather options
  opts <- list()

  # allow up to 30GB uploads
  opts$shiny.maxRequestSize <- 30*1024*1024^2

  # for violin plot dot size
  opts$shiny.usecairo <- FALSE

  # record shinytest2
  if (test) {
    if (!requireNamespace("shinytest2", quietly = TRUE))
      stop('install shinytest2')

    opts$shiny.testmode <- TRUE

    app <- shinytest2::AppDriver$new(app_dir, name = 'blah', options = opts)
    shinytest2::record_test(app, allow_no_input_binding = TRUE)
    return(NULL)

  } else {
    opts$shiny.testmode <- FALSE
  }

  if (missing(data_dir)) stop('data_dir not specified.')

  # pass arguments to app through options then run
  shiny::shinyOptions(
    data_dir = normalizePath(data_dir, mustWork = FALSE),
    app_name = app_name,
    pert_query_dir = normalizePath(pert_query_dir, mustWork = FALSE),
    pert_signature_dir = normalizePath(pert_signature_dir, mustWork = FALSE),
    gs_dir = normalizePath(gs_dir, mustWork = FALSE),
    indices_dir = normalizePath(indices_dir, mustWork = FALSE),
    tx2gene_dir = normalizePath(tx2gene_dir, mustWork = FALSE),
    tabs = tabs,
    logout_url = logout_url,
    is_local = is_local,
    is_example = is_example)

  # if developing
  # options(shiny.error = browser)

  # partial stack sometimes obscures errors
  opts$shiny.fullstacktrace <- FALSE

  # auto-reload if update app files
  is_aws <- !is.null(logout_url)
  if (!is_aws) opts$shiny.autoreload <- TRUE

  # set options and run app
  options(opts)
  shiny::runApp(app_dir, launch.browser = TRUE, host = host, port = port)
}


#' Initialize dseqr folders/files for a new app
#'
#' Creates necessary folders/files for a new dseqr app inside of /srv/shiny-server/dseqr.
#'
#' @param user_name User name for new dseqr app.
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
init_dseqr <- function(user_name, data_dir = '/srv/dseqr') {

  default_dir <- file.path(data_dir, user_name, 'default')
  dir.create(default_dir, recursive = TRUE)
  dir.create(file.path(default_dir, 'bulk'))
  dir.create(file.path(default_dir, 'single-cell'))
  dir.create(file.path(default_dir, 'custom_queries'))

}

run_dseqr_shinyproxy <- function(project_name, host_url) {

  is_example <- project_name == 'example'

  logout_url <- sprintf('https://%s/logout', host_url)

  message('project_name: ', project_name)
  message('is_example: ', is_example)
  message('logout_url: ', logout_url)

  # where to download/load drug and reference data from dseqr.data
  Sys.setenv('DSEQR_DATA_PATH' = '/srv/dseqr/.data')

  run_dseqr(project_name,
            logout_url = logout_url,
            is_example = is_example)

}
