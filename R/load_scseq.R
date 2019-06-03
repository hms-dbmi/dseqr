#' Load alevin quantification into a Seurat or SingleCellExperiment object.
#'
#' @param  Directory with raw and alevin-quantified single-cell RNA-Seq files.
#' @param  type Object type to return. Either \code{'Seurat'} or \code{'SingleCellExperiment'}.
#'
#' @return \code{Seurat} (default) or \code{SingleCellExperiment} with alevin whitelist meta data.
#' @export
#'
#' @examples
#'
#' data_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'
#' load_scseq(data_dir)
#'
load_scseq <- function(data_dir, type = 'Seurat') {

  # import alevin quants
  alevin_dir <- file.path(data_dir, 'alevin_output', 'alevin')
  counts <- tximport::tximport(file.path(alevin_dir, 'quants_mat.gz'), type = 'alevin')$counts

  # final alevin whitelist
  whitelist <- read.delim1(file.path(alevin_dir, 'whitelist.txt'))
  whitelist <- data.frame(whitelist = colnames(counts) %in% whitelist, row.names = colnames(counts))

  # covert to Seurat object
  srt <- Seurat::CreateSeuratObject(counts, meta.data = whitelist)
  if (type == 'Seurat') return(srt)

  # convert to SingleCellExperiment
  return(srt_to_sce(srt))
}

#' Convert Seurat object to SingleCellExperiment
#'
#' Also adds mrna and rrna
#'
#' @param srt
#'
#' @return
#' @export
#'
#' @examples
srt_to_sce <- function(srt) {
  sce <- Seurat::as.SingleCellExperiment(srt)

  qcgenes <- load_scseq_qcgenes()
  sce@metadata$mrna <- qcgenes$mrna
  sce@metadata$rrna <- qcgenes$rrna
  return(sce)
}

#' Load mitochondrial and ribsomal gene names
#'
#' This are the genes used by alevin for whitelisting.
#'
#' @return Named list with \code{rrna} and \code{mrna} character vectors.
#' @export
#'
#' @examples
load_scseq_qcgenes <- function() {

  # load mito and ribo genes
  rrna <- read.delim1(system.file('extdata', 'rrna.csv', package = 'drugseqr'))
  mrna <- read.delim1(system.file('extdata', 'mrna.csv', package = 'drugseqr'))

  return(list(rrna=rrna, mrna=mrna))
}

#' Normalize single-cell libraries for cell-specific biases
#'
#' @param sce \code{SingleCellExperiment} returned from \code{\link{load_scseq}}.
#'
#' @return Size-factor normalized \code{SingleCellExperiment}.
#' @export
#'
#' @examples
norm_scseq <- function(sce) {
  # paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0947-7
  # example: https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/tenx.html#6_normalizing_for_cell-specific_biases

  set.seed(1000)
  clusters <- scran::quickCluster(sce, use.ranks=FALSE, BSPARAM=BiocSingular::IrlbaParam())
  sce <- scran::computeSumFactors(sce, cluster=clusters, min.mean=0.1)
  sce <- scater::normalize(sce)
  return(sce)
}


#' Model and account for mean-variance relationship
#'
#' @param sce
#'
#' @return
#' @export
#'
#' @examples
stabilize_scseq <- function(sce) {
  # model mean-variance relationship
  sctech <- scran::makeTechTrend(x=sce)

  # remove sctech portion (keeps only biological component)
  set.seed(1000)
  sce <- scran::denoisePCA(sce, technical=sctech, BSPARAM=BiocSingular::IrlbaParam())
  return(sce)
}

#' Calculate QC metrics for SingleCellExperiment
#'
#' @param sce
#'
#' @return
#' @export
#'
#' @examples
qc_scseq <- function(sce) {
  # calculate qc metrics if haven't previous
  if (is.null(sce$total_counts))
    sce <- scater::calculateQCMetrics(sce,
                                      feature_controls=list(mito = which(row.names(sce) %in% sce@metadata$mrna),
                                                            ribo = which(row.names(sce) %in% sce@metadata$rrna)))
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


#' Add SNN Graph based clusters to SingleCellExperiment
#'
#' @param sce \code{SingleCellExperiment} returned by \code{\link{stabilize_scseq}} (contains \code{sce@reducedDims$PCA})
#'
#' @return \code{sce} with column \code{cluster} in \code{colData(sce)}
#' @export
#'
#' @examples
add_scseq_clusters <- function(sce, use.dimred = 'PCA') {
  snn.gr <- scran::buildSNNGraph(sce, use.dimred=use.dimred)
  clusters <- igraph::cluster_walktrap(snn.gr)
  sce$cluster <- factor(clusters$membership)

  return(sce)
}
