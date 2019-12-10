
#' Run surrogate variable analysis and factor out surrogate variables
#'
#' Returned result usefull for heatmaps, MDS plots, machine learning, etc
#'
#' @param eset ExpressionSet with group column in \code{pData(eset)}
#'
#' @return matrix of vsd (if RNA-seq) RMA (microarray) normalized expression with
#'  surrogate variables factored out
#' @export
remove_sva <- function(eset) {

  # determine if this is rna seq data
  rna_seq <- 'norm.factors' %in% colnames(Biobase::pData(eset))

  # make full and null model matrix
  group_levels = unique(Biobase::pData(eset)$group)
  group <- factor(Biobase::pData(eset)$group, levels = group_levels)
  pair <- Biobase::pData(eset)$pair


  mod <- stats::model.matrix(~0 + group)
  mod0 <- stats::model.matrix(~1, data = group)

  colnames(mod)[1:length(group_levels)] <- group_levels


  # remove duplicated rows (from 1:many PROBE:SYMBOL) as affect sva
  PROBE <- Biobase::fData(eset)[,1]
  expr <- unique(data.table::data.table(Biobase::exprs(eset), PROBE))[, PROBE := NULL]
  expr <- as.matrix(expr)

  # sva or svaseq
  sva_fun <-ifelse(rna_seq, sva::svaseq, sva::sva)
  svobj <- sva_fun(expr, mod, mod0)

  # get DESeq::vsd transformed counts for RNA-Seq
  el <- ifelse(rna_seq, 'vsd', 'exprs')
  vsd <- Biobase::assayDataElement(eset, el)

  # return cleaned
  vsd_clean <- clean_y(vsd, mod, svobj$sv)
  return(vsd_clean)

}
