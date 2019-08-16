#' Load kallisto/bustools quantification into a Seurat or SingleCellExperiment object.
#'
#' @param data_dir Directory with raw and kallisto/bustools or CellRanger quantified single-cell RNA-Seq files.
#' @param project String identifying sample.
#' @param type Quantification file type. One of either \code{'kallisto'} or \code{'cell_ranger'}.
#' @param h5 Boolean indicating, for \code{type = 'cell_ranger'}, if \code{data_dir} container a hdf5 file.
#' @param soupx Boolean indicating if \code{SoupX} should be used to remove backgroud counts (defaul is \code{FALSE}).
#'    Experimental and not currently recommended.
#'
#' @return \code{Seurat} object with whitelist meta data.
#' @export
#'
#' @examples
#'
#' data_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'
#' load_scseq(data_dir)
#'
load_scseq <- function(data_dir, project = 'SeuratProject', type = c('kallisto', 'cell_ranger'), h5 = FALSE, soupx = FALSE) {

  # load counts
  if (type[1] == 'kallisto') {
    data_dir <- file.path(data_dir, 'bus_output')
    counts <- load_kallisto_counts(data_dir)

  } else if (type[1] == 'cell_ranger') {
    counts <- load_cell_ranger_counts(data_dir, h5 = h5)
  }

  # generate/load whitelist
  whitelist <- get_scseq_whitelist(counts, data_dir)
  whitelist <- data.frame(whitelist = colnames(counts) %in% whitelist, row.names = colnames(counts))
  kneelist  <- readLines(file.path(data_dir, 'kneelist.txt'))

  # get ambient expression profile/determine outlier genes
  pct_ambient <- get_pct_ambient(counts)
  out_ambient <- get_outliers(pct_ambient)

  if (soupx) {
    empty <- !whitelist[[1]]
    counts <- strain_scseq(counts, empty)
    kneelist <- kneelist[kneelist %in% colnames(counts)]
  }


  # covert to Seurat object
  srt <- Seurat::CreateSeuratObject(counts[, kneelist], meta.data = whitelist, project = project)

  # add ambient metadata for genes
  srt[['RNA']]@meta.features$pct_ambient <- pct_ambient
  srt[['RNA']]@meta.features$out_ambient <- out_ambient

  return(srt)
}

#' Determine ambient percent for each gene
#'
#' Looks at droplets with counts less than or equal to 10 calculates the total percent for each gene
#'
#' @param counts \code{dgTMatrix} of counts. Rows are genes, columns are droplets.
#'
#' @return Named numeric vector of percentages for each gene.
#' @export
#' @keywords internal
get_pct_ambient <- function(counts) {

  # get drops with less than 10 counts
  ncount <- Matrix::colSums(counts)
  ambient <- counts[, ncount <= 10]

  # percentage of counts per gene
  nambient <- Matrix::rowSums(ambient)
  pct_ambient <- nambient / sum(nambient) * 100
  return(pct_ambient)
}

#' Flag outliers
#'
#'
#' @param x Named numeric vector
#'
#' @return Boolean vector with \code{length(x)} indicating if values of \code{x} are outliers (TRUE) or not (FALSE).
#' @export
#' @keywords internal
get_outliers <- function(x) {
  outliers <- graphics::boxplot(x, plot = FALSE)$out
  is.outlier <- names(x) %in% names(outliers)
  return(is.outlier)
}


#' Read kallisto/bustools market matrix and annotations
#'
#'
#' @inheritParams load_scseq
#'
#' @return sparse dgTMatrix with barcodes in columns and genes in rows.
#' @export
load_kallisto_counts <- function(data_dir) {

  # read sparse matrix
  counts <- Matrix::readMM(file.path(data_dir, 'genecount', 'genes.mtx'))
  counts <- Matrix::t(counts)

  # read annotations
  row.names(counts) <- readLines(file.path(data_dir, 'genecount', 'genes.genes.txt'))
  colnames(counts) <- readLines(file.path(data_dir, 'genecount', 'genes.barcodes.txt'))

  # remove non-expressed genes
  counts <- counts[Matrix::rowSums(counts) > 0, ]

  return(counts)
}

#' Load cell ranger counts
#'
#' Mainly to avoid having to download massive datasets that have already been quantified.
#'))))
#' @inheritParams load_scseq
#' @param h5
#' @importFrom magrittr "%>%"
#'
#' @return dgCMatrix
#' @export
load_cell_ranger_counts <- function(data_dir, h5) {
  # read the data in using ENSG features
  if (h5) {
    h5file <- list.files(data_dir, '.h5$', full.names = TRUE)
    counts <- Seurat::Read10X_h5(h5file, use.names = FALSE)

  } else {
    counts <- Seurat::Read10X(data_dir, gene.column = 1)
  }

  if (class(counts) == 'list') counts <- counts$`Gene Expression`

  counts <- counts[row.names(counts) %in% tx2gene$gene_id, ]

  # name genes by tx2gene so that best match with cmap/l1000 data
  map <- tx2gene %>%
    dplyr::filter(gene_id %in% row.names(counts)) %>%
    dplyr::select(gene_id, gene_name) %>%
    dplyr::distinct() %>%
    dplyr::arrange(match(gene_id, row.names(counts)))

  stopifnot(setequal(row.names(counts), map$gene_id))

  row.names(counts) <- map$gene_name

  # remove non-expressed genes
  counts <- counts[Matrix::rowSums(counts) > 0, ]

  # sum counts in rows with same gene
  counts <- Matrix.utils::aggregate.Matrix(counts, row.names(counts), fun = 'sum')

  return(counts)
}


#' Convert Seurat object to SingleCellExperiment
#'
#' Also adds mrna and rrna
#'
#' @param srt \code{Seurat} object
#'
#' @return \code{SingleCellExperiment} object
#' @export
srt_to_sce <- function(srt, assay = NULL) {
  if (class(srt) == 'SingleCellExperiment') return(srt)

  sce <- as.SingleCellExperiment(srt, assay)

  # for compatibility in explore_scseq_clusters
  sce$cluster <- sce$seurat_clusters

  return(sce)
}


#' Coerce Seurat to SingleCellExperiment
#'
#' This exists because of bug satijalab/seurat#1626. Remove once release accounts for it.
#'
#' @param x \code{Seurat} object
#' @param assay The assay to use. Default (NULL) uses the default assay from \code{x}.
#'
#' @return \code{SingleCellExperiment}
#'
as.SingleCellExperiment <- function(x, assay = NULL) {
  if (!Seurat:::PackageCheck('SingleCellExperiment', error = FALSE)) {
    stop("Please install SingleCellExperiment from Bioconductor before converting to a SingeCellExperiment object")
  }
  assay <- ifelse(is.null(assay), Seurat::DefaultAssay(object = x), assay)

  assays = list(
    counts = Seurat::GetAssayData(object = x, assay = assay, slot = "counts"),
    logcounts = Seurat::GetAssayData(object = x, assay = assay, slot = "data")
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
load_scseq_qcgenes <- function() {

  # load mito and ribo genes
  rrna <- readLines(system.file('extdata', 'rrna.csv', package = 'drugseqr', mustWork = TRUE))
  mrna <- readLines(system.file('extdata', 'mrna.csv', package = 'drugseqr', mustWork = TRUE))

  return(list(rrna=rrna, mrna=mrna))
}

add_qc_genes <- function(sce) {
  # add qc genes as metadata
  qcgenes <- load_scseq_qcgenes()
  sce@metadata$mrna <- qcgenes$mrna
  sce@metadata$rrna <- qcgenes$rrna

  return(sce)
}


#' Utility wrapper to run normalization and variance stabilization
#'
#' If \code{scseq} is a \code{SingleCellExperiment} object then runs the simpleSingleCell workflow.
#' If \code{scseq} is a \code{Seurat} object then uses \code{SCTransform}.
#'
#' @param scseq \code{Seurat} or \code{SingleCellExperiment} object
#'
#' @return Normalized and variance stabilized \code{scseq}.
#' @export
preprocess_scseq <- function(scseq, use_sctransform = TRUE) {

  if (use_sctransform) {
    scseq <- Seurat::SCTransform(scseq, verbose = FALSE)
  }

  else {
    scseq <- Seurat::NormalizeData(scseq)
    scseq <- Seurat::FindVariableFeatures(scseq)
    scseq <- Seurat::ScaleData(scseq)
  }

  return(scseq)
}


#' Normalize single-cell libraries for cell-specific biases
#'
#' @param sce \code{SingleCellExperiment} returned from \code{\link{load_scseq}}.
#'
#' @return Size-factor normalized \code{SingleCellExperiment}.
#' @export
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
#' @param sce \code{SingleCellExperiment}
#'
#' @return Variance stabilized \code{sce}
#' @export
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
#' @param sce \code{SingleCellExperiment}
#'
#' @return \code{sce} with qc metrics added by \code{\link[scater]{calculatedQCMetrics}}
#' @export
add_scseq_qc_metrics <- function(sce) {
  # calculate qc metrics if haven't previous
  if (!is.null(sce$total_counts)) return(sce)

  sce <- add_qc_genes(sce)
  sce <- scater::calculateQCMetrics(sce,
                                    feature_controls=list(mito = which(row.names(sce) %in% sce@metadata$mrna),
                                                          ribo = which(row.names(sce) %in% sce@metadata$rrna)))
  return(sce)
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
add_scseq_clusters <- function(scseq, reduction = 'pca', resolution = 0.8, dims = 1:30) {

  if (reduction == 'pca') suppressWarnings(scseq <- Seurat::RunPCA(scseq, verbose = FALSE))

  scseq <- Seurat::FindNeighbors(scseq, reduction = reduction, dims=dims, verbose = FALSE)
  scseq <- Seurat::FindClusters(scseq, verbose = FALSE, resolution = resolution)


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
get_scseq_markers <- function(scseq, ident.1 = NULL, ident.2 = NULL, min.diff.pct = 0.25, only.pos = TRUE, ...) {

  assay <- ifelse('SCT' %in% names(scseq@assays), 'SCT', 'RNA')

  # dont get markers if no clusters
  if (!exist_clusters(scseq)) return(NULL)

  # only upregulated as more useful for positive id of cell type
  if (!is.null(ident.1) | !is.null(ident.2)) {
    markers <- Seurat::FindMarkers(scseq,
                                   assay = assay,
                                   ident.1 = ident.1,
                                   ident.2 = ident.2,
                                   only.pos = only.pos,
                                   min.diff.pct = min.diff.pct,
                                   ...)

  } else {
    markers <- Seurat::FindAllMarkers(scseq,
                                      assay = assay,
                                      only.pos = only.pos,
                                      min.diff.pct = min.diff.pct,
                                      ...)
    markers <- split(markers, markers$cluster)
    markers <- lapply(markers, function(df) {row.names(df) <- df$gene; return(df)})
  }

  return(markers)
}




#' Integrate multiple scRNA-seq samples
#'
#' @param scseqs List of \code{Seurat} objects
#'
#' @return Integrated \code{Seurat} object with default assay of \code{"integrated"}
#' @export
integrate_scseqs <- function(scseqs, scalign = FALSE) {

  genes  <- Seurat::SelectIntegrationFeatures(object.list = scseqs, nfeatures = 3000)
  scseqs <- Seurat::PrepSCTIntegration(object.list = scseqs, anchor.features = genes)

  ambient <- get_integrated_ambient(scseqs)
  k.filter <- min(200, min(sapply(scseqs, ncol)))

  anchors <- Seurat::FindIntegrationAnchors(scseqs, k.filter = k.filter, normalization.method = "SCT",
                                            anchor.features = anchor.features)

  rm(scseqs); gc()

  combined <- Seurat::IntegrateData(anchors, normalization.method = "SCT")
  combined$orig.ident <- factor(combined$orig.ident)

  # add ambient outlier info
  combined <- add_integrated_ambient(combined, ambient)

  return(combined)
}

scalign_scseqs <- function(scseqs, genes) {

  sce.list <- lapply(scseqs, function(x) {
    SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = x[['SCT']]@data[genes, ],
                                                             scale.data = x[['SCT']]@scale.data[genes, ]))
  })

  sce.object = scAlign::scAlignCreateObject(sce.objects = sce.list, project.name = "sjia")

  sce.object = scAlign::scAlignMulti(sce.object,
                                     options=scAlign::scAlignOptions(steps=5000, log.every=5000, norm=TRUE, early.stop=FALSE, architecture="small"),
                                     encoder.data="scale.data",
                                     decoder.data="logcounts",
                                     supervised='none',
                                     run.encoder=TRUE,
                                     run.decoder=TRUE,
                                     log.results=TRUE,
                                     log.dir=file.path('./tmp'),
                                     device="CPU")

  return(sce.object)
}



get_integrated_ambient <- function(scseqs) {

  # datasets that are test samples
  is.test <- sapply(scseqs, function(x) levels(x$orig.ident) == 'test')

  # genes that are ambient in at least one test sample
  ambient.test <- lapply(scseqs[is.test], function(x) x[['RNA']]@meta.features)
  ambient.test <- lapply(ambient.test, function(x) row.names(x)[x$out_ambient])
  ambient.test <- unique(unlist(ambient.test))

  # genes that are ambient in at least one ctrl sample
  ambient.ctrl <- lapply(scseqs[!is.test], function(x) x[['RNA']]@meta.features)
  ambient.ctrl <- lapply(ambient.ctrl, function(x) row.names(x)[x$out_ambient])
  ambient.ctrl <- unique(unlist(ambient.ctrl))

  return(list(test = ambient.test, ctrl = ambient.ctrl))
}

#' Mark ambient outliers in combined dataset
#'
#' A gene is marked as an ambient outlier if it is an ambient outlier in at least one of the datasets.
#'
#' @param scseqs the original scseqs
#' @param combined the combined scseqs
#'
#' @return \code{combined} with \code{out_ambient} column added to \code{meta.features} slot of \code{SCT} assay.
#' @export
#' @keywords internal
add_integrated_ambient <- function(combined, ambient) {

  # genes to keep in order that appear
  keep <- row.names(combined[['SCT']])

  # set in combined
  combined[['SCT']]@meta.features$test_ambient <- FALSE
  combined[['SCT']]@meta.features$ctrl_ambient <- FALSE
  combined[['SCT']]@meta.features$test_ambient[keep %in% ambient$test] <- TRUE
  combined[['SCT']]@meta.features$ctrl_ambient[keep %in% ambient$ctrl] <- TRUE

  return(combined)
}

#' Get predicted cells types for a query dataset based on a reference dataset
#'
#' @param reference \code{Seurat} object to use as a reference
#' @param query \code{Seurat} object to obtain predictions for
#'
#' @return \code{data.frame} with predictions
#' @export
transfer_labels <- function(reference, query) {

  k.filter <- min(201, ncol(reference), ncol(query))-1

  anchors <- Seurat::FindTransferAnchors(reference, query, k.filter = k.filter)
  predictions <- Seurat::TransferData(anchorset = anchors,
                                      refdata = reference$seurat_clusters,
                                      dims = 1:30)

  return(predictions)
}


#' Test is there is at lest two clusters
#'
#' Used by \code{\link{get_scseq_markers}} to prevent getting markers if there are no clusters to compare
#'
#' @param scseq
#'
#' @return TRUE if more than one cluster exists
#' @export
exist_clusters <- function(scseq) {
  if (class(scseq) == 'SingleCellExperiment') {
    exist_clusters <- length(unique(scseq$clusters)) > 1

  } else if (class(scseq) == 'Seurat') {
    exist_clusters <- length(unique(Seurat::Idents(scseq))) > 1
  }
  return(exist_clusters)
}

#' Run UMAP for visualizing single cell data.
#'
#' If \code{scseq} is a \code{SingleCellExperiment} object then uses \code{scater::runUMAP}.
#' If \code{scseq} is a \code{Seurat} object then uses \code{Seurat::RunUMAP}.
#'
#' @param scseq \code{SingleCellExperiment} or \code{Seurat} object.
#'
#' @return \code{scseq} with TSNE results.
#' @export
run_umap <- function(scseq, dims = 1:30, reduction = 'pca') {

  set.seed(1000)
  if (class(scseq) == 'SingleCellExperiment') {
    scseq <- scater::runUMAP(scseq, use_dimred="PCA")

  } else if (class(scseq) == 'Seurat') {
    scseq <- Seurat::RunUMAP(scseq, dims = dims, reduction = reduction, verbose = FALSE)

  } else {
    stop('scseq must be either class SingleCellExperiment or Seurat')
  }

  return(scseq)
}

#' Add jitter to UMAP embeddings
#'
#' @param scseq \code{Seurat} object with umap reduction
#' @inheritParams base::jitter
#'
#' @return \code{scseq} with jitter added to umap embedding
#' @export
jitter_umap <- function(scseq, factor = 1, amount = 0) {
  scseq[['umap']]@cell.embeddings <-
    apply(scseq[['umap']]@cell.embeddings, 2, jitter, factor, amount)

  return(scseq)
}

