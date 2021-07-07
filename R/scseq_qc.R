
#' Determine non-empty droplets
#'
#' @param counts dgCMatrix of counts
#'
#' @return Indices of columns of counts corresponding to non-empty droplets
#' @keywords internal
#'
detect_cells <- function(counts, qcgenes) {

  # omit mito/ribo genes for cell calling (see MarioniLab/DropletUtils#33)
  keep.genes <- !row.names(counts) %in% unlist(qcgenes)

  set.seed(100)
  e.out <- DropletUtils::emptyDrops(counts[keep.genes, ], retain = Inf)
  keep.cells <- which(e.out$FDR <= 0.001)

  return(keep.cells)
}

#' Run quality control for single-cell dataset
#'
#' @param sce \code{SingleCellExperiment}.
#' @param metrics Character vector of metrics to remove outliers for.
#'
#' @return \code{sce} with outliers removed.
#' @keywords internal
run_scseq_qc <- function(sce, metrics = c('low_lib_size',
                                          'low_n_features',
                                          'high_subsets_mito_percent',
                                          'low_subsets_ribo_percent',
                                          'high_doublet_score')) {

  # remove low lib size/high mito
  df <- sce@colData
  reasons <- scater::quickPerCellQC(df, percent_subsets=c("subsets_mito_percent"))

  # remove low ribo
  reasons$low_subsets_ribo_percent <- scater::isOutlier(df$subsets_ribo_percent, type = 'lower')

  # remove high doublet score
  reasons$high_doublet_score <- df$doublet_class == 'doublet'

  # allow to see outliers for non-selected metrics/combination of all
  not_selected <- setdiff(colnames(reasons), c('discard', metrics))

  if (length(not_selected)) {
    reasons$outlier_any <- apply(reasons, 1, any)
    sce@colData <- cbind(df, reasons[, c(not_selected, 'outlier_any')])
  }

  if (!is.null(metrics)) {
    # subset to specified metrics
    reasons <- reasons[, colnames(reasons) %in% metrics, drop = FALSE]

    # update discard reasons
    reasons$discard <- apply(reasons, 1, any)
    message('keeping ', sum(!reasons$discard), '/', ncol(sce), ' non-empty droplets.')

    sce <- sce[ ,!reasons$discard]
  }

  return(sce)
}
