#' Covert from crossmeta to dseqr formats
#'
#' Convert a microarray ExpressionSet saved with 'crossmeta' into a format
#' for direct usage in 'dseqr'.
#'
#' @param data_dir directory that contains \code{gse_name} folder
#' @param gse_name Name of GSE accession dataset, must be folder in data_dir
#'
#' @export
#' @examples
#' library(Biobase)
#'
#' # generate example eset
#' sd <- 0.3*sqrt(4/rchisq(1000,df=4))
#' y <- matrix(rnorm(1000*6,sd=sd),1000,6)
#' rownames(y) <- paste("Gene",1:1000)
#' y[1:2,4:6] <- y[1:2,4:6] + 2
#'
#' y[100:200, c(1,2,6)] <- y[100:200, c(1,2,6)] + 0.5
#'
#' pdata <- data.frame(group = c(rep('healthy', 3), rep('disease', 3)))
#' fdata <- data.frame(SYMBOL = row.names(y),
#'                    PROBE = row.names(y),
#'                    row.names = row.names(y))
#'
#' eset <- ExpressionSet(y,
#'                      phenoData = as(pdata, 'AnnotatedDataFrame'),
#'                      featureData = as(fdata, 'AnnotatedDataFrame'))
#'
#'
#' # save
#' eset <- list(GSE1=eset)
#' data_dir <- tempdir()
#' gse_dir <- file.path(data_dir, 'GSE1')
#' dir.create(gse_dir)
#' saveRDS(eset, file.path(gse_dir, 'GSE1_eset.rds'))
#'
#' from_crossmeta('GSE1', data_dir)
#'
from_crossmeta <- function(gse_name, data_dir) {

  # change saved eset name to eset.rds
  dataset_dir <- file.path(data_dir, gse_name)
  eset_path <- list.files(dataset_dir, '^.+?_eset.rds', full.names = TRUE)
  eset <- readRDS(eset_path)[[1]]
  saveRDS(eset, file.path(dataset_dir, 'eset.rds'))
}
