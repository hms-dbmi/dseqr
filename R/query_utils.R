#' Get dprime effect size values.
#'
#' These are used to query against drug effect size matrices.
#'
#' @param top_table data.frame of differential expression results.
#' @return Named numeric vector. Names are gene names, values are effect size values.
#' @keywords internal
#'
#' @examples
#'
#' # generate fake previous result
#' top_table <- data.frame(dprimes=rnorm(5), row.names = paste0('gene', 1:5))
#'
#' dprimes <- dseqr:::get_dprimes(top_table)
#'
get_dprimes <- function(top_table) {
  dprimes <- top_table$dprime
  names(dprimes) <- row.names(top_table)
  return(dprimes)
}


#' Get correlation between query and drug signatures.
#'
#' Determines the pearson correlation between the query and each drug signature.
#'
#' Drugs with the largest positive and negative pearson correlation are predicted to,
#' respectively, mimic and reverse the query signature. Values range from +1 to -1.
#'
#'
#' @param query_genes Named numeric vector of differential expression values for
#'   query genes. Names are HGNC symbols.
#' @param drug_es A matrix of differential expression values for drugs.
#' @param ngenes The number of top differentially-regulated (up and down) query genes to use.
#'
#' @return Named vector of pearson correlations between query and drug combination signatures.
#'
#' @keywords internal
#'
#' @examples
#'
#' # load CMAP02 data
#' cmap_es <- dseqr.data::load_data('cmap_es_ind.qs')
#'
#' # use first cmap_es signature as query
#' dprimes <- cmap_es[1:100, 1]
#' names(dprimes) <- row.names(cmap_es)[1:100]
#'
#' # get correlations between query and drug signatures
#' res <- dseqr:::query_drugs(dprimes, cmap_es)
#'
query_drugs <- function(query_genes, drug_es, ngenes = 200) {

  # use only common genes
  query_genes <- query_genes[names(query_genes) %in% row.names(drug_es)]

  # top up/down ngenes
  top_ngenes  <- utils::head(names(sort(abs(query_genes), TRUE)), ngenes)
  query_genes <- query_genes[top_ngenes]
  drug_es <- drug_es[names(query_genes), ,drop = FALSE]

  if (length(query_genes) == 1) {
    # single gene: sort by opposing sign
    sim <- sort(drug_es[1,], decreasing = query_genes < 0)
    sim <- scales::rescale(sim, to = c(-1,1))
  } else {
    # pearson correlation
    sim <- stats::cor(query_genes, drug_es, method="pearson")
    sim <- structure(c(sim), names=colnames(sim))
  }

  return(sim)
}

#' Find drugs that maximally affect a selection of genes
#'
#' Results are based on the average effect on the query genes. Genes to upregulate are multiplied by -1 so that
#' strong positive results for these genes contribute towards a more negative average effect. All results are divided
#' by the absolute of the average minimum effect (after -1 multiplication of genes to upregulate) of the query genes.
#' This ensures that results will be between -1 and 1 for consistency with correlation values for \code{\link{query_drugs}}.
#'
#' @param query_genes Named list of character vectors with \code{'dn'} indicating genes that want to down-regulated and
#'   \code{'up'} indicating genes that want to up-regulate.
#' @inheritParams query_drugs
#'
#' @return Named numeric vector where most negative results are predicted to have the strongest desired effect as
#'  indicated by \code{query_genes}.
#'
#' @keywords internal
query_budger <- function(query_genes, drug_es) {

  # will sort by increasing so multiply upregulated genes by -1
  is.up <- row.names(drug_es) %in% query_genes$up
  drug_es[is.up, ] <- drug_es[is.up, ] * -1

  # use only common genes
  query_genes <- unlist(query_genes, use.names = FALSE)
  query_genes <- query_genes[query_genes %in% row.names(drug_es)]

  # put results between -1 and 1 with most negative result having most desired effect
  # so that consistent with query_drugs
  norm <- apply(drug_es[query_genes,, drop = FALSE], 1, min)
  sim <- colMeans(drug_es[query_genes,, drop = FALSE]) / abs(mean(norm))

  return(sim)
}


