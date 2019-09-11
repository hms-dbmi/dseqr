#' Validate pdata before differential expression analysis
#'
#' @param pdata data.frame with column \code{'Group'}
#'
#' @return NULL if valid, otherwise a character vector indicating what's wrong
#' @export
#' @keywords internal
validate_pdata <- function(pdata) {
  group <- pdata$Group
  group <- group[!is.na(group)]

  if (length(unique(group)) != 2) {
    msg <- 'Analysis requires test and control groups'

  } else if (length(group) < 3) {
    msg <- 'At least three samples are required for analysis'

  } else {
    msg <- NULL
  }
  return(msg)
}


#' Load previous bulk anals dataframe
#'
#' @param data_dir Path to folder container \code{'bulk'} and \code{'single-cell'} directories
#' @return data.frame with columns "dataset_name" "dataset_dir" and "anal_name".
#' @export
#' @keywords internal
load_bulk_anals <- function(data_dir, with_type = FALSE) {
  anals_path <- file.path(data_dir, 'bulk', 'anals.rds')

  if (file.exists(anals_path)) {
    anals <- readRDS(anals_path)

  } else {
    anals <- data.frame(matrix(ncol = 3, nrow = 0), stringsAsFactors = FALSE)
    colnames(anals) <- c("dataset_name", "dataset_dir", "anal_name")
    saveRDS(anals, anals_path)
  }

  anals$label <- anals$anal_name
  anals$value <- seq_len(nrow(anals))

  if (with_type) {
    anals$type <- anals$dataset_name
    if (nrow(anals)) anals$type <- paste0('Bulk - ', anals$type)

  }

  return(anals)
}

#' Save new analysis info to anals dataframe
#'
#' @param dataset_name Name of dataset
#' @param dataset_dir Folder name for dataset
#' @param anal_name Name of new analysis
#' @param data_dir Path to folder container \code{'bulk'} and \code{'single-cell'} directories
#' @return NULL
#' @export
#' @keywords internal
save_bulk_anals <- function(dataset_name, dataset_dir, anal_name, data_dir) {
  anals_path <- file.path(data_dir, 'bulk', 'anals.rds')
  anals <- readRDS(anals_path)

  anals[nrow(anals)+1, ] <- c(dataset_name, dataset_dir, anal_name)
  saveRDS(anals, anals_path)
}

#' Load previous bulk datasets dataframe
#'
#' @param data_dir Path to folder container \code{'bulk'} and \code{'single-cell'} directories
#' @return data.frame with columns "dataset_name" and "dataset_dir"
#' @export
#' @keywords internal
load_bulk_datasets <-function(data_dir) {
  datasets_path <- file.path(data_dir, 'bulk', 'datasets.rds')

  if (file.exists(datasets_path)) {
    datasets <- readRDS(datasets_path)

  } else {
    datasets <- data.frame(matrix(ncol = 2, nrow = 0), stringsAsFactors = FALSE)
    colnames(datasets) <- c("dataset_name", "dataset_dir")
    saveRDS(datasets, datasets_path)
  }

  datasets$value <-  datasets$label <- datasets$dataset_name

  return(datasets)
}

#' Save new dataset info to datasets dataframe
#'
#' @param dataset_name Name of dataset
#' @param dataset_dir Folder name for dataset
#' @param data_dir Path to folder container \code{'bulk'} and \code{'single-cell'} directories
#' @return NULL
#' @export
#' @keywords internal
save_bulk_dataset <- function(dataset_name, dataset_dir, data_dir) {
  datasets_path <- file.path(data_dir, 'bulk', 'datasets.rds')
  datasets <- readRDS(datasets_path)

  datasets[nrow(datasets)+1, ] <- c(dataset_name, dataset_dir)
  saveRDS(datasets, datasets_path)
}

