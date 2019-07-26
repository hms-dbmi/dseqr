
#' Get percentage of cells expressing each gene
#'
#' @param scseq \code{Seurat} object
#' @param ident.1 Test group level in \code{Idents(scseq)}.
#' @param ident.2 Control group level in \code{Idents(scseq)}.
#'
#' @return data.frame with columns pct.1 and pct.2 indicating fraction of cells that express each gene.
#' @export
#' @keywords internal
get_cell_pcts <- function(scseq, ident.1, ident.2) {
  data <- scseq[['SCT']]@data

  cells <- Seurat::Idents(scseq)
  cells.1 <- cells == ident.1
  cells.2 <- cells == ident.2

  pct.1 <- round(
    x = Matrix::rowSums(x = data[, cells.1, drop = FALSE] > 0) / sum(cells.1),
    digits = 3
  )

  pct.2 <- round(
    x = Matrix::rowSums(x = data[, cells.2, drop = FALSE] > 0) / sum(cells.2),
    digits = 3
  )

  return(cbind(pct.1, pct.2))
}



#' Utility to generate filename for single cell download csv
#'
#' @param cluster Character vector of cluster names
#' @param anal result of diff_expr_scseq
#' @param comparison_type either 'samples' or 'clusters'
#'
#' @return File name string
#' @export
#' @keywords internal
sc_dl_filename <- function(cluster, anal, comparison_type) {

  # remove underscores and spaces for cluster(s) and analysis name
  cluster <- gsub('[_ ]', '-', cluster)
  anal <- gsub('[_ ]', '-', anal)

  if (comparison_type == 'samples')
    cluster <- paste0('test-vs-ctrl_', paste(cluster, collapse = '-'))


  paste('single-cell', anal, cluster, paste0(Sys.Date(), '.csv'), sep='_')
}



#' Get cluster choices data.frame for selectize dropdown
#'
#' @param clusters Character vector of cluster names
#' @param scseq \code{Seurat} object
#' @param value Character vector for value column which is returned from \code{selectizeInput}. Default is \code{clusters}.
#'
#' @return data.frame with columns for rendering selectizeInput cluster choices
#' @export
#' @keywords internal
get_cluster_choices <- function(clusters, scseq, value = clusters) {

  # show the cell numbers/percentages
  ncells <- tabulate(scseq$seurat_clusters)
  pcells <- round(ncells / sum(ncells) * 100)
  pspace <- strrep('&nbsp;&nbsp;', 2 - nchar(pcells))

  # cluster choices are the clusters themselves
  testColor <- get_palette(clusters)
  cluster_choices <- data.frame(name = stringr::str_trunc(clusters, 27),
                                value = value,
                                label = clusters,
                                testColor,
                                ncells, pcells, pspace, row.names = NULL)

  return(cluster_choices)
}

#' Get contrast choices data.frame for selectize dropdown
#'
#' @param clusters Character vector of cluster names
#' @param test Name of test contrast
#'
#' @return data.frame with columns for rendering selectizeInput contrast choices
#' @export
#' @keywords internal
get_contrast_choices <- function(clusters, test) {

  # group choices are as compared to other clusters
  ctrls <- clusters[clusters != test]

  colours <- get_palette(clusters)
  names(colours) <- clusters

  contrast_choices <- data.frame(test = stringr::str_trunc(test, 11, ellipsis = '..'),
                                 ctrl = stringr::str_trunc(c('all', ctrls), 11, ellipsis = '..'),
                                 value = c(test, paste0(test, ' vs ', ctrls)),
                                 testColor = colours[test],
                                 ctrlColor = c('white', colours[ctrls]), row.names = NULL)

  return(contrast_choices)

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
integrate_saved_scseqs <- function(data_dir, test, ctrl, anal_name, updateProgress = NULL) {

  # save dummy data if testing shiny
  if (isTRUE(getOption('shiny.testmode'))) {
    scseq_data <- list(scseq = NULL, markers = NULL, annot = NULL)
    save_scseq_data(scseq_data, anal_name, data_dir, integrated = TRUE)
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
  scseq_data <- list(scseq = combined, markers = markers, annot = names(markers))
  save_scseq_data(scseq_data, anal_name, data_dir, integrated = TRUE)
}

#' Save Single Cell RNA-seq data for app
#'
#' @param scseq_data Named list with \code{scseq}, \code{markers}, and/or \code{annot}
#' @param anal_name The analysis name.
#' @param data_dir Path to directory to save in
#' @param integrated is the analysis integration. Default is \code{FALSE}
#'
#' @return NULL
#' @export
save_scseq_data <- function(scseq_data, anal_name, data_dir, integrated = FALSE, reduce_size = FALSE) {
  if (integrated) {
    int_path <- file.path(data_dir, 'integrated.rds')
    int_options <- readRDS(int_path)
    saveRDS(c(int_options, anal_name), int_path)
  }


  if (reduce_size) {
    # optionally seperate larger parts off
    scseq <- scseq_data$scseq
    sct <- scseq[['SCT']]

    # keep @data which has corrected log counts for visualization
    scseq[['SCT']]@scale.data <- matrix(nrow = 0, ncol = 0)
    scseq[['SCT']]@counts <- matrix(nrow = 0, ncol = 0)
    scseq[['SCT']]@misc <- NULL

    # redundant
    stopifnot(all.equal(scseq[['RNA']]@counts, scseq[['RNA']]@data))
    scseq[['RNA']]@counts <- matrix(nrow = 0, ncol = 0)

    scseq_data$scseq <- scseq
    scseq_data$sct <- sct
  }


  dir.create(file.path(data_dir, anal_name))
  for (type in names(scseq_data)) {
    saveRDS(scseq_data[[type]], scseq_part_path(data_dir, anal_name, type))
  }

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
scseq_part_path <- function(data_dir, anal_name, part) {
  fname <- paste0(anal_name, '_', part, '.rds')
  file.path(data_dir, anal_name, fname)
}
