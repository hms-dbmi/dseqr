#' Covert from crossmeta to drugseqr formats
#'
#' @param data_dir directory that contains \code{gse_name} folder
#' @param gse_name Name of GSE accession dataset, must be folder in data_dir
#'
#' @return
#' @export
#'
#' @examples
#' data_dir <- 'data-raw/patient_data/covid19/bulk'
#' gse_name <- 'GSE17400'
from_crossmeta <- function(gse_name, data_dir) {

  # change saved eset name to eset.rds
  dataset_dir <- file.path(data_dir, gse_name)
  eset_path <- list.files(dataset_dir, '^.+?_eset.rds', full.names = TRUE)
  eset <- readRDS(eset_path)[[1]]
  saveRDS(eset, file.path(dataset_dir, 'eset.rds'))
}
