#' Load single-cell RNA-Seq data into a SingleCellExperiment.
#'
#' @param  Directory with raw and quantified single-cell RNA-Seq files.
#'
#' @return \code{SingleCellExperiment} with \code{colData} attribute including alevin whitelist.
#' @export
#'
#' @examples
#'
#' data_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'
#' load_scseq(data_dir)
#'
load_scseq <- function(data_dir) {

  # import alevin quants
  alevin_dir <- file.path(data_dir, 'alevin_output', 'alevin')
  counts <- tximport::tximport(file.path(alevin_dir, 'quants_mat.gz'), type = 'alevin')$counts

  # final alevin whitelist
  whitelist <- read.delim1(file.path(alevin_dir, 'whitelist.txt'))

  # load mito and ribo genes
  rrna <- read.delim1(system.file('extdata', 'rrna.csv', package = 'drugseqr'))
  mrna <- read.delim1(system.file('extdata', 'mrna.csv', package = 'drugseqr'))

  # Convert to SingleCellExperiment
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts),
                                                    colData = data.frame(whitelist = colnames(counts) %in% whitelist),
                                                    metadata = list(rrna = rrna, mrna = mrna))

  return(sce)
}

#' Helper function to read single column text files
#'
#' @param file File to read from
#'
#' @return Character vector from \code{file}.
#' @export
#'
#' @examples
read.delim1 <- function(file) {
  return(read.delim(file, header = FALSE, as.is = TRUE)$V1)
}
