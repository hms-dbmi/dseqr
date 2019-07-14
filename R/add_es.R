#' Add metaMA effectsize values to top table.
#'
#' Adds moderated unbiased standardised effect sizes (dprimes) to top table
#' from differential expression analysis.
#'
#' @param diff_exprs Result from call to \code{\link{diff_expr}}.
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

add_es <- function(diff_exprs, cols = c("dprime", "vardprime")) {


  # get study degrees of freedom and group classes

  df <- diff_exprs$ebayes_sv$df.residual + diff_exprs$ebayes_sv$df.prior
  classes <- diff_exprs$pdata$group

  # group names for contrast
  groups <- c('test', 'ctrl')

  # get sample sizes for groups
  ni <- sum(classes == groups[2])
  nj <- sum(classes == groups[1])

  # bind effect size values with top table
  tt <- diff_exprs$top_table
  es <- metaMA::effectsize(tt$t, ((ni * nj)/(ni + nj)), df)[, cols, drop = FALSE]
  tt <- cbind(tt, es)

  # store results
  diff_exprs$top_table <- tt

  return(diff_exprs)
}
