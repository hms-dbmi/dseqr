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
load_scseq <- function(data_dir, type = 'Seurat', project = 'SeuratProject') {

  # import alevin quants
  alevin_dir <- file.path(data_dir, 'alevin_output', 'alevin')
  counts <- tximport::tximport(file.path(alevin_dir, 'quants_mat.gz'), type = 'alevin')$counts

  # final alevin whitelist
  whitelist <- read.delim1(file.path(alevin_dir, 'whitelist.txt'))
  whitelist <- data.frame(whitelist = colnames(counts) %in% whitelist, row.names = colnames(counts))

  # covert to Seurat object
  srt <- Seurat::CreateSeuratObject(counts, meta.data = whitelist, project = project)
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
srt_to_sce <- function(srt, assay = NULL) {
  if (class(srt) == 'SingleCellExperiment') return(srt)

  sce <- as.SingleCellExperiment(srt, assay)

  # add qc genes as metadata
  qcgenes <- load_scseq_qcgenes()
  sce@metadata$mrna <- qcgenes$mrna
  sce@metadata$rrna <- qcgenes$rrna

  # Seuray assay used by prevent_integrated
  sce@metadata$seurat_assay <- ifelse(is.null(assay), Seurat::DefaultAssay(srt), assay)

  # for compatibility in explore_scseq_clusters
  sce$cluster <- sce$seurat_clusters

  return(sce)
}

#' Coerce Seurat to SingleCellExperiment
#'
#' This exists because of bug satijalab/seurat#1626
#'
#' @param x
#' @param assay
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
as.SingleCellExperiment <- function(x, assay = NULL, ...) {
  if (!Seurat:::PackageCheck('SingleCellExperiment', error = FALSE)) {
    stop("Please install SingleCellExperiment from Bioconductor before converting to a SingeCellExperiment object")
  }
  assay <- ifelse(is.null(assay), Seurat::DefaultAssay(object = x), assay)

  assays = list(
    counts = Seurat::GetAssayData(object = x, assay = assay, slot = "counts"),
    logcounts = Seurat::GetAssayData(object = x, assay = assay, slot = "data")
  )

  assays <- assays[sapply(assays, nrow) != 0]
  sce <- SingleCellExperiment::SingleCellExperiment(assays = assays)

  metadata <- x[[]]
  metadata$ident <- Seurat::Idents(object = x)
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(metadata)
  SummarizedExperiment::rowData(sce) <- S4Vectors::DataFrame(x[[assay]][[]])
  for (dr in Seurat:::FilterObjects(object = x, classes.keep = "DimReduc")) {
    SingleCellExperiment::reducedDim(sce, toupper(x = dr)) <- Seurat::Embeddings(object = x[[dr]])
  }
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

#' Utility wrapper to run normalization and variance stabilization
#'
#' If \code{scseq} is a \code{SingleCellExperiment} object then runs the simpleSingleCell workflow.
#' If \code{scseq} is a \code{Seurat} object then uses \code{SCTransform}.
#'
#' @param scseq
#'
#' @return
#' @export
#'
#' @examples
preprocess_scseq <- function(scseq) {

  if (class(scseq) == 'SingleCellExperiment') {
    scseq <- norm_scseq(scseq)
    scseq <- stabilize_scseq(scseq)
  }

  if (class(scseq) == 'Seurat') {
    # alevin has non-integer values that sctransform turns to Inf
    scseq <- Seurat::SetAssayData(scseq, 'counts', round(Seurat::GetAssayData(scseq, 'counts')))
    scseq <- Seurat::SCTransform(scseq, verbose = FALSE, return.only.var.genes = FALSE)
  }
  return(scseq)
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
add_scseq_clusters <- function(scseq, use.dimred = 'PCA') {

  if (class(scseq) == 'SingleCellExperiment') {
    snn.gr <- scran::buildSNNGraph(scseq, use.dimred=use.dimred)
    clusters <- igraph::cluster_walktrap(snn.gr)
    scseq$cluster <- factor(clusters$membership - 1)

  } else if (class(scseq) == 'Seurat') {
    suppressWarnings(scseq <- Seurat::RunPCA(scseq, verbose = FALSE))
    scseq <- Seurat::FindNeighbors(scseq, dims=1:30, verbose = FALSE)
    scseq <- Seurat::FindClusters(scseq, verbose = FALSE)
  } else {
    stop('scseq must be either class SingleCellExperiment or Seurat')
  }

  return(scseq)
}

#' Get markers genes for single cell clusters
#'
#' If \code{scseq} is a \code{SingleCellExperiment} object then uses \code{scran::findMarkers}.
#' If \code{scseq} is a \code{Seurat} object then uses \code{Seurat::FindAllMarkers}.
#'
#' @param scseq \code{SingleCellExperiment} or \code{Seurat} object.
#' @param assay.type Used by \code{\link[scran]{findMarkers}}
#'
#' @return List of \code{data.frame}s, one for each cluster.
#' @export
#'
#' @examples
get_scseq_markers <- function(scseq, assay.type = 'logcounts') {

  if (!exist_clusters(scseq)) return(NULL)
  scseq <- prevent_integrated(scseq)

  # only upregulated as more useful for positive id of cell type
  if (class(scseq) == 'SingleCellExperiment') {
    markers <- scran::findMarkers(scseq, clusters=scseq$cluster, direction="up", assay.type = assay.type)

  } else if (class(scseq) == 'Seurat') {
    suppressWarnings(markers <- Seurat::FindAllMarkers(scseq, only.pos = TRUE, verbose = FALSE))
    markers <- split(markers, markers$cluster)
    markers <- lapply(markers, function(df) {row.names(df) <- df$gene; return(df)})
  }

  return(markers)
}

#' Switch away from integrated assay slot
#'
#' The "integrated" assay slot is inappropriate for differential expression analysis. This will switch to \code{SCTransform}'ed
#' slot for \code{Seurat} objects. \code{SingleCellExperiment} objects check to make sure they weren't generated from
#' "integrated" data.
#'
#' @param scseq \code{Seurat} or \code{SingleCellExperiment} object.
#'
#' @return \code{scseq} with "SCT" as the default assay if it was previously "integrated".
#' @export
#'
#' @examples
prevent_integrated <- function(scseq) {

  if (class(scseq) == 'Seurat' &&
      Seurat::DefaultAssay(srt) == 'integrated') {
    Seurat::DefaultAssay(srt) <- 'SCT'

  } else if (class(scseq) == 'SingleCellExperiment' &&
             isTRUE(scseq@meta.data$seurat_assay) == 'integrated') {
    stop("SingleCellExperiment object was generated from integrated Seurat assay.")
  }

  return(scseq)
}

exist_clusters <- function(scseq) {
  if (class(scseq) == 'SingleCellExperiment') {
    exist_clusters <- length(unique(scseq$clusters)) > 1

  } else if (class(scseq) == 'Seurat') {
    exist_clusters <- length(unique(Seurat::Idents(scseq))) > 1
  }
  return(exist_clusters)
}

#' Run TSNE for visualizing single cell data.
#'
#' If \code{scseq} is a \code{SingleCellExperiment} object then uses \code{scater::runTSNE}.
#' If \code{scseq} is a \code{Seurat} object then uses \code{Seurat::RunTSNE}.
#'
#' @param scseq \code{SingleCellExperiment} or \code{Seurat} object.
#'
#' @return \code{scseq} with TSNE results.
#' @export
#'
#' @examples
run_tsne <- function(scseq) {

  set.seed(1000)
  if (class(scseq) == 'SingleCellExperiment') {
    scseq <- scater::runTSNE(scseq, use_dimred="PCA")

  } else if (class(scseq) == 'Seurat') {
    scseq <- Seurat::RunTSNE(scseq, dims = 1:30, verbose = FALSE)

  } else {
    stop('scseq must be either class SingleCellExperiment or Seurat')
  }

  return(scseq)
}



