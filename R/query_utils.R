#' Get dprime effect size values.
#'
#' These are used to query against drug effect size matrices.
#'
#' @param diff_exprs List returned by \code{\link{dif_expr}}
#'
#' @return
#' @export
#'
#' @examples
#'
#' # load result of previous differential expression analysis
#' data_dir <- file.path('data-raw/example-data')
#' anal <- readRDS(file.path(data_dir, 'diff_expr_symbol.rds'))
#'
#' dprimes <- get_dprimes(anal)
#'
get_dprimes <- function(diff_exprs) {
  diff_exprs <- add_es(diff_exprs)

  dprimes <- diff_exprs$top_table$dprime
  names(dprimes) <- row.names(diff_exprs$top_table)

  return(dprimes)
}
