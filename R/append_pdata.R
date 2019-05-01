#' Join CMAP02/L1000 pdata and query results.
#'
#' @param query_res Named numeric vector. Result of previous call to \code{\link{query_drugs}}
#' @param study Character vector. Either \code{'CMAP02'} or \code{'L1000'}.
#'
#' @return \code{data.frame} with correlation column showing results from \code{query_res} and other columns giving information about each treatment.
#' @export
#'
#' @examples
#'
#' # load CMAP02 data
#' cmap_path <- system.file('extdata', 'cmap_es_ind.rds', package = 'drugseqr')
#' cmap_es <- readRDS(cmap_path)
#'
#' # load previous analysis
#' anal_path <- file.path('data-raw/example-data', 'diff_expr_symbol.rds')
#' anal <- readRDS(anal_path)
#'
#' # get effect size values
#' dprimes <- get_dprimes(anal)
#'
#' # run query
#' res <- query_drugs(dprimes, cmap_es)
#'
#' # append drug pdata
#' res <- append_pdata(res, 'CMAP02')
#'
append_pdata <- function(query_res, study) {

  # load pdata for study
  pdata_fname <- paste0(study, '_pdata.rds')
  pdata_path <- system.file('extdata', pdata_fname, package = 'drugseqr')
  pdata <- readRDS(pdata_path)

  if (!all(pdata$title == names(query_res)))
    stop("Query results and pdata titles don't match.")

  # add result
  query_res <- tibble::add_column(pdata, Correlation = query_res, .before=0)

  # destructure title and return
  query_res <- destructure_title(query_res, .after = 'Correlation')
  return(query_res)
}

#' Split CMAP02/L1000 title into multiple columns.
#'
#' Function used by \code{\link{append_pdata}}.
#'
#' @param pdata \code{data.frame} with and \code{'title'} columns. \code{'title'} column
#'   should contain 4 underscore seperated parts (e.g. \code{'10-DEBC_A375_20um_24h'}).
#' @param drop Should the \code{'title} column be droped? Default is \code{TRUE}.
#' @param ... Additional arguments to \code{\link[tibble]{add_column}} (e.g. \code{.before} or \code{.after}).
#'
#' @return \code{pdata} without \code{'title'} column and 4 new columns: \code{'Compound'}, \code{'Cell Line'}, \code{'Dose'},
#'   and \code{'Duration'}.
#' @export
#'
#' @examples
destructure_title <- function(pdata, drop = TRUE, ...) {
  title_split <- strsplit(pdata$title, '_')

  pdata <- tibble::add_column(pdata,
                              'Compound'  = sapply(title_split, `[`, 1),
                              'Cell Line' = sapply(title_split, `[`, 2),
                              'Dose'      = sapply(title_split, `[`, 3),
                              'Duration'  = sapply(title_split, `[`, 4), ...)

  if (drop) pdata$title <- NULL
  return(pdata)
}
