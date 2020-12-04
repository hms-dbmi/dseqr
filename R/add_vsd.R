
#' Add VST normalized assay data element to expression set
#'
#' For microarray datasets duplicates exprs slot into vsd slot.
#'
#' @param eset ExpressionSet with group column in \code{pData(eset)}
#' @param vsd_path Path to save result to. Allows skipping running transform
#'   on each load.
#'
#' @return \code{eset} with \code{'vsd'} \code{assayDataElement} added.
#' @export
add_vsd <- function(eset, rna_seq = TRUE, pbulk = FALSE, vsd_path = NULL) {

  # for cases where manually added (e.g. nanostring dataset)
  els <- Biobase::assayDataElementNames(eset)
  if ('vsd' %in% els) return(eset)

  if (!rna_seq) {
    # for microarray use exprs
    vsd <- Biobase::assayDataElement(eset, 'exprs')

  } else if (!is.null(vsd_path) && file.exists(vsd_path)) {
    vsd <- readRDS(vsd_path)

  } else if (pbulk) {
    pdata <- Biobase::pData(eset)
    dds <- DESeq2::DESeqDataSetFromMatrix(Biobase::exprs(eset), pdata, design = ~group)
    dds <- DESeq2::estimateSizeFactors(dds)
    vsd <- tryCatch(DESeq2::rlog(dds, blind = FALSE),
                    warning = function(e) {
                      if (grepl('varianceStabilizingTransformation', e$message))
                        DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
                    })

    vsd <- SummarizedExperiment::assay(vsd)
    if (!is.null(vsd_path)) saveRDS(vsd, vsd_path)

  } else {
    # rlog fast enough for most rna-seq
    vsd <- GEOkallisto::get_vsd(eset)
    vsd <- SummarizedExperiment::assay(vsd)
    if (!is.null(vsd_path)) saveRDS(vsd, vsd_path)
  }

  Biobase::assayDataElement(eset, 'vsd') <- vsd
  return(eset)
}


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

is.installed <- function(packages, level = c('quiet', 'error', 'message')) {
  level <- level[1]

  for (package in packages) {

    if (!requireNamespace(package, quietly = TRUE)) {
      msg <- paste0("Package \"", package, "\" needed for this function to work. Please install it.")
      if (level == 'error') stop(msg, call. = FALSE)
      if (level == 'message') message(msg)
      return(FALSE)
    }
  }

  return(TRUE)
}
