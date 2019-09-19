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
#' data_dir <- '~/Documents/Batcave/zaklab/drugseqr/data-raw/patient_data/example'
#' pert_query_dir <- '~/Documents/Batcave/zaklab/drugseqr/data-raw/drug_gene_queries/data'
#'
#' run_drugseqr(data_dir, pert_query_dir, test_data = TRUE)
#' run_drugseqr(data_dir, pert_query_dir, test_data = FALSE)
#'
run_drugseqr <- function(data_dir, pert_query_dir = NULL, test = FALSE, test_data = FALSE) {

  app_dir <- 'inst/app'

  # pass arguments to app through options then run
  shiny::shinyOptions(data_dir = data_dir, pert_query_dir = pert_query_dir)


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
#' init_drugseqr('example', local_dir = 'data-raw/patient_data')
#'
init_drugseqr <- function(app_name, local_dir = NULL) {

  if (is.null(local_dir)) {
    # sync the drugseqr app components
    drugseqr_dir <- system.file('app', package = 'drugseqr', mustWork = TRUE)
    app_dir <- file.path('/srv/shiny-server/drugseqr', app_name)
    dir.create(app_dir, recursive = TRUE)

    system2('rsync', args = c('-av', paste0(drugseqr_dir, '/'), paste0(app_dir, '/')))

    # create necessary folders/blank files to initialize new app
    data_dir <- file.path(app_dir, 'data_dir')
  } else {
    data_dir <- file.path(local_dir, app_name)
  }

  dir.create(data_dir)
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
