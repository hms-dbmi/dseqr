#' Construct Pseudobulk QC expression sets
#'
#' Used for testing if significant changes in e.g. sum of mtRNA or ribosomal counts
#'
#' @param pbulk_esets List of pseudobulk expression sets
#'
#' @return List of pseodobulk qc expression sets
#' @keywords internal
#'
#'
construct_pbulk_qcsets <- function(pbulk_esets) {
  qc <- load_scseq_qcgenes()

  qcsets <- list()
  for (i in seq_along(pbulk_esets)) {
    eset <- pbulk_esets[[i]]
    mset <- eset[row.names(eset) %in% qc$mrna, ]
    rset <- eset[row.names(eset) %in% qc$rrna, ]

    Biobase::exprs(eset)[1, ] <- colSums(Biobase::exprs(mset))
    Biobase::exprs(eset)[2, ] <- colSums(Biobase::exprs(rset))

    eset <- eset[1:2, ]
    fnames <- c('mito_sum', 'ribo_sum')
    row.names(eset) <- fnames
    Biobase::fData(eset) <- data.frame(gene_name = fnames, ENTREZID = fnames, row.names = fnames)
    qcsets[[as.character(i)]] <- eset
  }

  return(qcsets)
}

