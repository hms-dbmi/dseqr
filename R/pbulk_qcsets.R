#' Construct Pseudobulk QC expression sets
#'
#' Used for testing if significant changes in e.g. sum of mtRNA or ribosomal counts
#'
#' @param pbulk_esets List of pseudobulk expression sets
#'
#' @return List of pseodobulk qc expression sets
#' @export
#'
#' @examples
#' pbulk_esets <- readRDS("~/Documents/Batcave/zaklab/drugseqr/data-raw/patient_data/sjia/single-cell/sjia_pbmcs_lungdisease_vs_healthy_harmony/pbulk_esets.rds")
#' qcsets <- construct_pbulk_qcsets(pbulk_esets)
#' lm_fit <- run_limma_scseq(qcsets)
#' get_top_table(lm_fit[[1]], robust = FALSE)
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

  return(qcset)
}
