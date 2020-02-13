#' Get whitelist for high quality cells
#'
#' @param counts dgTMatrix of counts for all barcodes.
#' @inheritParams load_scseq
#'
#' @return Character vector of barcodes called as high quality cells.
#' @export
#' @keywords internal
get_scseq_whitelist <- function(counts, data_dir, overwrite = TRUE, species = 'Homo sapiens') {

  # check for previous whitelist
  whitelist_path <- file.path(data_dir, 'whitelist.txt')
  if (file.exists(whitelist_path) & !overwrite) {
    whitelist <- readLines(whitelist_path)
    return(whitelist)
  }

  ncount <- Matrix::colSums(counts)

  # if already filtered cellranger matrix, can't call if empty
  if (min(ncount) < 10) {

    # omit mito/ribo genes for cell calling (see MarioniLab/DropletUtils#33)
    qc <- load_scseq_qcgenes(species)
    keep <- !row.names(counts) %in% unlist(qc)

    set.seed(100)
    e.out <- DropletUtils::emptyDrops(counts[keep, ], retain = Inf)
    counts <- counts[, which(e.out$FDR <= 0.001)]

  } else {
    message('Skipping emptyDrops as counts previously filtered.')
  }

  # add qc metrics
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  sce <- add_scseq_qc_metrics(sce, species = species)

  # setup stats for outlier detection
  df <- sce@colData
  stats <- cbind(log10(df$sum), log10(df$detected), df$subsets_mito_percent, df$subsets_ribo_percent)

  outlying <- robustbase::adjOutlyingness(stats, only.outlyingness = TRUE)
  outl.drop <- scater::isOutlier(outlying, type = "higher")

  kneelist  <- colnames(sce)
  whitelist <- colnames(sce)[!outl.drop]
  message('keeping ', length(whitelist), '/', length(kneelist), ' non-empty droplets.')

  writeLines(whitelist, whitelist_path)
  writeLines(kneelist, file.path(data_dir, 'kneelist.txt'))
  return(whitelist)
}
