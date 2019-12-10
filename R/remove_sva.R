
#' Run surrogate variable analysis and factor out surrogate variables and any pairs
#'
#' Returned result useful for heatmaps, MDS plots, machine learning, etc
#'
#' @param eset ExpressionSet with group column in \code{pData(eset)} and optional pair column.
#' @param run_sva If \code{FALSE} (default \code{TRUE}) then only regresses out pair column.
#'
#' @return matrix of vsd (if RNA-seq) RMA (microarray) normalized expression with
#'  surrogate variables and optional pairs factored out
#' @export
remove_sva <- function(eset, run_sva = TRUE) {

  # determine if this is rna seq data
  rna_seq <- 'norm.factors' %in% colnames(Biobase::pData(eset))

  # make full and null model matrix
  group_levels = unique(Biobase::pData(eset)$group)
  group <- factor(Biobase::pData(eset)$group, levels = group_levels)
  pair <- factor(as.numeric(Biobase::pData(eset)$pair))

  if (length(pair)) {
    mod <- stats::model.matrix(~0 + group + pair)
    mod0 <- stats::model.matrix(~1 + pair, data = group)
    colnames(mod)[1:length(group_levels)] <- group_levels

    # remove non matched pairs
    has.pair <- colSums(mod0[, 2:ncol(mod0)]) >= 2
    has.pair <- names(has.pair[has.pair])

    mod <- mod[, c(group_levels, has.pair)]
    mod0 <- mod0[, c("(Intercept)", has.pair)]

  } else {
    mod <- stats::model.matrix(~0 + group)
    mod0 <- stats::model.matrix(~1, data = group)
    colnames(mod)[1:length(group_levels)] <- group_levels
  }

  # remove duplicated rows (from 1:many PROBE:SYMBOL) as affect sva
  PROBE <- Biobase::fData(eset)[,1]
  expr <- unique(data.table::data.table(Biobase::exprs(eset), PROBE))[, PROBE := NULL]
  expr <- as.matrix(expr)

  # sva or svaseq
  sva_fun <-ifelse(rna_seq, sva::svaseq, sva::sva)
  if (run_sva) {
    svobj <- sva_fun(expr, mod, mod0)
  } else {
    svobj <- list(sv = NULL)
  }

  # get DESeq::vsd transformed counts for RNA-Seq
  el <- ifelse(rna_seq, 'vsd', 'exprs')
  vsd <- Biobase::assayDataElement(eset, el)

  # clean svs and pairs
  if (length(pair)) {
    to_clean <- cbind(svobj$sv, mod0[, has.pair])
  } else {
    to_clean <- svobj$sv
  }

  vsd_clean <- clean_y(vsd, mod[, group_levels], to_clean)
  return(vsd_clean)
}
