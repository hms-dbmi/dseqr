#' Add metaMA effectsize values to top table.
#'
#' Adds moderated unbiased standardised effect sizes (dprimes) to top table
#' from differential expression analysis.
#'
#' @param diff_exprs Result from call.
#' @param cols Columns from \code{\link[metaMA]{effectsize}} result to add to
#'    top table.
#'
#' @export
#'
#' @return diff_exprs with specified columns added to top_table.
#'
#' @examples
#' # location of previously saved analysis
#' data_dir <- system.file('extdata', 'IBD', package = 'drugseqr')
#'
#' # load previous analysis for eset
#' anal <- readRDS(file.path(data_dir, 'diff_exprs_symbol.rds'))
#'
#' # add dprime and vardprime to top tables
#' anal <- add_es(anal)
add_es <- function(tt, ebfit, cols = c("dprime", "vardprime"), groups = c('test', 'ctrl')) {


  # get study degrees of freedom and group classes
  df <- ebfit$df.residual + ebfit$df.prior

  # get sample sizes for groups
  ni <- sum(ebfit$design[, groups[2]])
  nj <- sum(ebfit$design[, groups[1]])

  # bind effect size values with top table
  es <- metaMA::effectsize(tt$t, ((ni * nj)/(ni + nj)), df)[, cols, drop = FALSE]
  tt <- cbind(tt, es)

  return(tt)
}
