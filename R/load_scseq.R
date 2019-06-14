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
  if (type == 'Seurat') {
    return(srt)

  } else if (type == 'SingleCellExperiment') {
    # convert to SingleCellExperiment
    return(srt_to_sce(srt))

  } else {
    stop('type must be either Seurat or SingleCellExperiment')
  }

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
    logcounts = Seurat::GetAssayData(object = x, assay = assay, slot = "data"),
    scale.data = Seurat::GetAssayData(object = x, assay = assay, slot = "scale.data")
  )

  assays <- assays[sapply(X = assays, FUN = nrow) != 0]
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

  } else {
    stop('scseq must be either a SingleCellExperiment or Seurat object.')
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


#' Add clusters single cell RNA-seq object
#'
#' Uses either \code{SingleCellExperiment} or \code{Seurat} workflows based on class of \code{scseq}.
#'
#' @param scseq \code{SingleCellExperiment} or \code{Seurat} returned by \code{\link{preprocess_scseq}}.
#' @param use.dimred A string specifying reduced dimension to use e.g. \code{'PCA'} (default).
#'  Only used if \code{scseq} has class \code{SingleCellExperiment}.
#' @param resolution Use larger values to obtain more clusters (default is 0.8).
#'  Only used if \code{scseq} has class \code{Seurat}.
#'
#' @return If \code{scseq} is a \code{SingleCellExperiemnt} object, column \code{cluster} in \code{colData(sce)} is added.
#'  If \code{scseq} is a \code{Seurat} object, the result of \code{\link[Seurat]{FindClusters}} is returned.
#' @export
#'
#' @examples
add_scseq_clusters <- function(scseq, use.dimred = 'PCA', resolution = 0.8) {

  if (class(scseq) == 'SingleCellExperiment') {
    snn.gr <- scran::buildSNNGraph(scseq, use.dimred=use.dimred)
    clusters <- igraph::cluster_walktrap(snn.gr)
    scseq$cluster <- factor(clusters$membership - 1)

  } else if (class(scseq) == 'Seurat') {
    suppressWarnings(scseq <- Seurat::RunPCA(scseq, verbose = FALSE))
    scseq <- Seurat::FindNeighbors(scseq, dims=1:30, verbose = FALSE)
    scseq <- Seurat::FindClusters(scseq, verbose = FALSE, resolution = resolution)

  } else {
    stop('scseq must be either class SingleCellExperiment or Seurat')
  }

  return(scseq)
}

#' Get markers genes for single cell clusters
#'
#' If \code{scseq} is a \code{SingleCellExperiment} object then uses \code{\link[scran]{findMarkers}}.
#' If \code{scseq} is a \code{Seurat} object then uses either \code{\link[Seurat]{FindAllMarkers}} or
#' \code{Seurat::FindMarkers} if both \code{ident.1} and \code{ident.2} are not NULL.
#'
#' @param scseq \code{SingleCellExperiment} or \code{Seurat} object.
#' @param assay.type Only used if \code{scseq} has class \code{SingleCellExperiment}.
#' A string specifying which assay values to use, e.g., \code{"counts"} or \code{"logcounts"} (default).
#' @param ident.1 Identity class to define markers for.
#'  Only used if \code{scseq} has class \code{Seurat} and \code{ident.2} is not \code{NULL}.
#' @param ident.2 A second identity class for comparison.
#'  Only used if \code{scseq} has class \code{Seurat} and \code{ident.1} is not \code{NULL}.
#'
#' @return List of \code{data.frame}s, one for each cluster.
#' @export
#'
#' @examples
get_scseq_markers <- function(scseq, assay.type = 'logcounts', ident.1 = NULL, ident.2 = NULL) {

  # dont get markers if no clusters
  if (!exist_clusters(scseq)) return(NULL)

  # only upregulated as more useful for positive id of cell type
  if (class(scseq) == 'SingleCellExperiment') {
    markers <- scran::findMarkers(scseq, clusters=scseq$cluster, direction="up", assay.type = assay.type)

  } else if (class(scseq) == 'Seurat') {
    if (!is.null(ident.1) & !is.null(ident.2)) {
      markers <- list()
      markers[[ident.1]] <- Seurat::FindMarkers(scseq, assay = 'SCT', only.pos = TRUE, ident.1 = ident.1, ident.2 = ident.2)

    } else {
      # no scale.data for integrated
      slot <- ifelse(nrow(Seurat::GetAssayData(scseq, slot = 'scale.data', assay = 'SCT')),
                     'scale.data', 'data')
      markers <- Seurat::FindAllMarkers(scseq, assay = 'SCT', slot = slot, only.pos = TRUE)
      markers <- split(markers, markers$cluster)
      markers <- lapply(markers, function(df) {row.names(df) <- df$gene; return(df)})
    }
  }

  return(markers)
}

#' Integrate multiple scRNA-seq samples
#'
#' @param scseqs List of \code{Seurat} objects
#'
#' @return Integrated \code{Seurat} object with default assay of \code{"integrated"}
#' @export
#'
#' @examples
integrate_scseq <- function(scseqs) {

  anchors <- Seurat::FindIntegrationAnchors(scseqs)
  combined <- Seurat::IntegrateData(anchors)
  Seurat::DefaultAssay(combined) <- "integrated"
  combined <- Seurat::ScaleData(combined, verbose = FALSE)

  return(combined)
}





#' Test is there is at lest two clusters
#'
#' Used by \code{\link{get_scseq_markers}} to prevent getting markers if there are no clusters to compare
#'
#' @param scseq
#'
#' @return TRUE if more than one cluster exists
#' @export
#'
#' @examples
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



