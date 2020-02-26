
#' Determine non-empty droplets
#'
#' @param counts dgCMatrix of counts
#' @param species Character vector with names of species. Supports Homo sapiens or Mus musculus.
#'
#' @return Indices of columns of counts corresponding to non-empty droplets
#' @export
#'
detect_cells <- function(counts, species = 'Homo sapiens') {

  # omit mito/ribo genes for cell calling (see MarioniLab/DropletUtils#33)
  qc <- load_scseq_qcgenes(species)
  keep.genes <- !row.names(counts) %in% unlist(qc)

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
#' @export
run_scseq_qc <- function(sce, metrics = c('low_lib_size',
                                          'low_n_features',
                                          'high_subsets_mito_percent',
                                          'low_subsets_ribo_percent',
                                          'high_doublet_score',
                                          'high_outlyingness')) {

  # remove low lib size/high mito
  df <- sce@colData
  reasons <- scater::quickPerCellQC(df, percent_subsets=c("subsets_mito_percent"))

  # remove low ribo
  reasons$low_subsets_ribo_percent <- scater::isOutlier(df$subsets_ribo_percent, type = 'lower')

  # remove high doublet score
  reasons$high_doublet_score <- scater::isOutlier(df$doublet_score, type = 'higher')

  # remove high outlyingness for QC metrics
  stats <- cbind(log10(df$sum), log10(df$detected), df$subsets_mito_percent, df$subsets_ribo_percent)
  outlying <- robustbase::adjOutlyingness(stats, only.outlyingness = TRUE)
  reasons$high_outlyingness <- scater::isOutlier(outlying, type = 'higher')

  # allow to see outliers for non-selected metrics/combination of all
  not_selected <- setdiff(colnames(reasons), c('discard', metrics))

  if (length(not_selected)) {
    reasons$discard <- apply(reasons, 1, any)
    sce@colData <- cbind(df, reasons[, c(not_selected, 'discard')])
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
