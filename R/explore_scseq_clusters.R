#' Explore Single Cell Clusters
#'
#' @param scseq \code{SingleCellExperiment} or \code{Seurat} object
#' @param markers Named list of \code{data.frame}s where \code{row.names} are the marker genes. One list per cluster in \code{scseq}
#'  with the same name as the cluster.
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' data_dir <- '~/Documents/Batcave/zaklab/drugseqr/data-raw/single-cell/example-anals/sjia'
#'
#' explore_scseq_clusters(data_dir)
#'

explore_scseq_clusters <- function(data_dir, test = FALSE, test_data = TRUE) {

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

  }

  # auto-reload if update app files
  options(shiny.autoreload = TRUE)
  shiny::runApp(app_dir, launch.browser = TRUE)

}

#' Integrate previously saved scseqs
#'
#' Performs integration and saves as a new analysis.
#' Used by \code{explore_scseq_clusters} shiny app.
#'
#' @param data_dir Directory with saved analyses.
#' @param test Character vector of test analysis names.
#' @param ctrl Character vector of control analysis names.
#' @param anal_name Name for new integrated analysis.
#' @param progress optional Shiny \code{Progress} object.
#'
#' @return NULL
#' @export
#' @keywords internal
#'
#' @examples
integrate_saved_scseqs <- function(data_dir, test, ctrl, anal_name, updateProgress = NULL) {

  # save dummy data if testing shiny
  if (isTRUE(getOption('shiny.testmode'))) {
    save_combined(combined = NULL, markers = NULL, data_dir = data_dir, anal_name = anal_name)
    return(NULL)
  }

  # default updateProgress and number of steps
  if (is.null(updateProgress)) updateProgress <- function(...) {NULL}
  n = 6

  # get paths for saved scseqs
  test_paths <- scseq_part_path(data_dir, test, 'scseq')
  ctrl_paths <- scseq_part_path(data_dir, ctrl, 'scseq')

  updateProgress(1/n, 'loading')
  test_scseqs <- lapply(test_paths, readRDS)
  ctrl_scseqs <- lapply(ctrl_paths, readRDS)

  # set orig.ident to ctrl/test and integrate
  test_scseqs <- lapply(test_scseqs, function(x) {x$orig.ident <- factor('test'); x})
  ctrl_scseqs <- lapply(ctrl_scseqs, function(x) {x$orig.ident <- factor('ctrl'); x})

  updateProgress(2/n, 'integrating')
  combined <- integrate_scseqs(c(test_scseqs, ctrl_scseqs))

  updateProgress(3/n, 'clustering')
  combined <- add_scseq_clusters(combined)

  updateProgress(4/n, 'reducing')
  combined <- run_umap(combined)

  updateProgress(5/n, 'getting markers')
  markers <- get_scseq_markers(combined)

  updateProgress(6/n, 'saving')
  save_combined(combined, markers, data_dir, anal_name)
}

save_combined <- function(combined, markers, data_dir, anal_name) {
  int_path <- file.path(data_dir, 'integrated.rds')
  int_options <- readRDS(int_path)
  saveRDS(c(int_options, anal_name), int_path)

  dir.create(file.path(data_dir, anal_name))
  saveRDS(combined, scseq_part_path(data_dir, anal_name, 'scseq'))
  saveRDS(markers, scseq_part_path(data_dir, anal_name, 'markers'))
  saveRDS(names(markers), scseq_part_path(data_dir, anal_name, 'annot'))

  return(NULL)
}

#' Validate dataset selection for integration
#'
#' @param test Character vector of test dataset names
#' @param ctrl Character vector of control dataset names
#'
#' @return \code{NULL} is valid, otherwise an error message
#' @export
#' @keywords internal
#'
#' @examples
validate_integration <- function(test, ctrl, anal_name, anal_options) {
  msg <- NULL
  # make sure both control and test analyses provided
  if (is.null(anal_name) || anal_name == '') {
    msg <- 'Provide a name for integrated analysis'

  } else if (anal_name %in% unlist(anal_options)) {
    msg <- 'Analysis name already exists'

  } else if (is.null(test) || is.null(ctrl)) {
    msg <- 'Need control and test datasets'
  }

  return(msg)
}


#' Get path to saved scseq part
#'
#' @param data_dir Path to directory with analyses.
#' @param anal_name Name of analysis.
#' @param part either \code{'annot'}, \code{'scseq'}, or \code{'markers'}.
#'
#' @return Path to analysis \code{part}.
#' @export
#' @keywords internal
#'
#' @examples
scseq_part_path <- function(data_dir, anal_name, part) {
  fname <- paste0(anal_name, '_', part, '.rds')
  file.path(data_dir, anal_name, fname)
}


