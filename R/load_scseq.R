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
load_scseq <- function(data_dir, project = 'SeuratProject', type = c('kallisto', 'cellranger'), soupx = FALSE) {

  # load counts
  if (type[1] == 'kallisto') {
    data_dir <- file.path(data_dir, 'bus_output')
    counts <- load_kallisto_counts(data_dir)

  } else if (type[1] == 'cellranger') {
    counts <- load_cellranger_counts(data_dir)
  }

  # generate/load whitelist
  whitelist <- get_scseq_whitelist(counts, data_dir)
  whitelist <- data.frame(whitelist = colnames(counts) %in% whitelist, row.names = colnames(counts))
  kneelist  <- readLines(file.path(data_dir, 'kneelist.txt'))

  # get ambient expression profile/determine outlier genes
  # if pre-filtered cellranger, can't determine outliers
  ncount <- Matrix::colSums(counts)
  if (min(ncount) > 10) {
    pct_ambient <- 0
    out_ambient <- FALSE

  } else {
    pct_ambient <- get_pct_ambient(counts)
    out_ambient <- get_outliers(pct_ambient)
  }

  if (soupx) {
    empty <- !whitelist[[1]]
    counts <- strain_scseq(counts, empty)
    kneelist <- kneelist[kneelist %in% colnames(counts)]
  }

  # covert to SingleCellExperiment object
  # add ambient metadata for genes
  rowData <- DataFrame(pct_ambient, out_ambient)
  colData <- DataFrame(project = rep(project, length(kneelist)))

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts[, kneelist]),
    rowData = rowData,
    colData = colData
  )

  return(sce)
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
#' @importFrom magrittr "%>%"
#'
#' @return dgCMatrix
#' @export
load_cellranger_counts <- function(data_dir) {
  # read the data in using ENSG features
  h5file <- list.files(data_dir, '.h5$', full.names = TRUE)
  if (length(h5file)) {
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

#' Determine if the selected folder has CellRanger files
#'
#' @param data_dir path to directory to check.
#'
#' @return \code{TRUE} if CellRanger files detected, otherwise \code{FALSE}.
#' @export
check_is_cellranger <- function(data_dir) {

  # cellranger file names
  files <- list.files(data_dir)
  mtx.file <- grep('.mtx', files, fixed = TRUE, value = TRUE)
  genes.file <- grep('features.tsv|genes.tsv', files, value = TRUE)
  barcodes.file <-  grep('barcodes.tsv', files, fixed = TRUE, value = TRUE)

  if (length(mtx.file) & length(genes.file) & length(barcodes.file)) return(TRUE)
  return(FALSE)
}

#' Rename CellRanger files for loading by Seurat::Read10X
#'
#' @param data_dir Path to folder with CellRanger files.
#'
#' @return NULL
#' @export
standardize_cellranger <- function(data_dir) {

  # cellranger file names
  files <- list.files(data_dir)
  mtx.file <- grep('.mtx', files, fixed = TRUE, value = TRUE)
  genes.file <- grep('features.tsv|genes.tsv', files, value = TRUE)
  barcodes.file <-  grep('barcodes.tsv', files, fixed = TRUE, value = TRUE)

  # rename for ?Seurat::Read10X
  file.rename(file.path(data_dir, c(mtx.file, genes.file, barcodes.file)),
              file.path(data_dir, c('matrix.mtx', 'genes.tsv', 'barcodes.tsv')))
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


#' Utility wrapper to run normalization and log transformation
#'
#'
#' @param scseq \code{SingleCellExperiment} object
#'
#' @return Normalized and log transformed \code{scseq}.
#' @export
normalize_scseq <- function(scseq) {

  ncores <- min(parallel::detectCores(), 7)
  doParallel::registerDoParallel(ncores)

  set.seed(100)
  preclusters <- scran::quickCluster(scseq, BPPARAM = BiocParallel::DoparParam())
  scseq <- scran::computeSumFactors(scseq, cluster=preclusters)
  scseq <- scater::logNormCounts(scseq)

  return(scseq)
}

#' Get top highly variable genes for subsequent dimensionality reduction
#'
#' Runs after \code{preprocess_scseq}
#'
#' @param sce \code{SingleCellExperiment} object
#'
#' @return
#' @export
add_hvgs <- function(sce) {

  dec <- scran::modelGeneVar(sce)
  hvg <- row.names(sce) %in% scran::getTopHVGs(dec, prop=0.1)
  SummarizedExperiment::rowData(sce)$hvg <- hvg

  return(sce)
}

#' Perform PCA and TSNE dimensionality reduction
#'
#' Runs after \code{add_hvgs} so that runs faster and on interesting genes.
#'
#' @param sce \code{SingleCellExperiment} object with \code{'hvg'} in \code{rowData}
#'
#' @return \code{sce} with \code{'PCA'} and \code{'TSNE'} reducedDim slots.
#' @export
reduce_dims <- function(sce) {

  # run on HVGs
  rdata <- SummarizedExperiment::rowData(sce)
  subset_row <- row.names(rdata[rdata$hvg, ])

  # run PCA and pick number of PCs
  set.seed(100)
  sce <- scater::runPCA(sce, subset_row = subset_row)
  sce <- pick_npcs(sce)

  # TSNE on top PCs
  sce <- run_tsne(sce)

  return(sce)
}

run_tsne <- function(sce, dimred = 'PCA') {
  set.seed(1100101001)
  sce <- scater::runTSNE(sce, dimred = dimred, n_dimred = sce@metadata$npcs)
  colnames(SingleCellExperiment::reducedDim(sce, 'TSNE')) <- c('TSNE1', 'TSNE2')
  return(sce)
}


add_scseq_clusters <- function(sce, dimred = 'PCA') {

  FUN <- function(x, ...) {
    g <- scran::buildSNNGraph(x, ..., transposed = TRUE)
    igraph::cluster_louvain(g)$membership
  }

  # pick number of PCs
  pcs <- SingleCellExperiment::reducedDim(sce, type = dimred)
  choices <- scran::getClusteredPCs(pcs, FUN = FUN)
  npcs <- S4Vectors::metadata(choices)$chosen
  sce@metadata$npcs <- npcs


  # add clusters
  sce$cluster <- factor(choices$clusters[[npcs]])
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




#' Run pairwise wilcox tests between single cell clusters
#'
#' @param scseq \code{SingleCellExperiment} object.
#'
#' @return List of \code{data.frame}s, one for each cluster.
#' @export
pairwise_wilcox <- function(scseq, groups = scseq$cluster, direction = 'up', block = NULL, restrict = NULL) {

  # dont get markers if no clusters
  if (!exist_clusters(scseq)) return(NULL)
  groups <- as.character(groups)

  # only upregulated as more useful for positive id of cell type
  wilcox_tests <- scran::pairwiseWilcox(SingleCellExperiment::logcounts(scseq),
                                        groups = groups,
                                        direction = direction,
                                        block = block,
                                        restrict = restrict)
  return(wilcox_tests)
}

#' Combine pairwise wilcox tests between single-cell clusters
#'
#' @param tests Result of \code{pairwise_wilcox}
#' @param pval.type
#'
#' @return List of data.frames
#' @export
get_scseq_markers <- function(tests, pval.type = 'some', effect.field = 'AUC', keep = NULL) {
  if (is.null(keep)) keep <- rep(TRUE, nrow(tests$pairs))

  markers <- scran::combineMarkers(tests$statistics[keep],
                                   tests$pairs[keep, ],
                                   pval.type = pval.type,
                                   effect.field = effect.field)
  lapply(markers, as.data.frame)
}








#' Integrate multiple scRNA-seq samples
#'
#' @param scseqs List of \code{SingleCellExperiment} objects
#' @param type One of \code{'clusterMNN'} (default) or \code{'fastMNN'} specifying cluster function to use.
#'
#' @return Integrated \code{SingleCellExperiment} object.
#' @export
integrate_scseqs <- function(scseqs, type = c('clusterMNN', 'fastMNN')) {

  # all common genes
  universe <- Reduce(intersect, lapply(scseqs, row.names))

  # variance modelling results
  decs <- lapply(scseqs, scran::modelGeneVar)

  # subset scseqs and decs
  decs <- lapply(decs, function(x) x[universe, ])
  scseqs <- lapply(scseqs, `[`, universe)

  # rescale each batch for depths
  scseqs <- do.call('multiBatchNorm', scseqs, envir = loadNamespace('batchelor'))

  # feature selection
  decs <- do.call('combineVar', decs, envir = loadNamespace('scran'))
  hvgs <- decs$bio > 0

  # mnn integration
  # TODO use fastMNN restriction to exclude batch specific cells
  combined <- do.call('noCorrect', scseqs, envir = loadNamespace('batchelor'))


  set.seed(1000101001)
  if (type[1] == 'clusterMNN') {
    mnn.fun <- function(...) batchelor::clusterMNN(
      ...,
      batch = combined$batch,
      clusters = lapply(scseqs, `[[`, 'cluster'),
      subset.row = hvgs,
      auto.merge = TRUE,
      correct.all = TRUE,
      cos.norm = FALSE)

  } else if (type[1] == 'fastMNN') {
    mnn.fun <- function(...) batchelor::fastMNN(
      ...,
      batch = combined$batch,
      subset.row = hvgs,
      auto.merge = TRUE,
      correct.all = TRUE,
      cos.norm = FALSE,
      prop.k = 0.05)
  }

  mnn.out <- do.call(mnn.fun, scseqs)
  mnn.out$orig.ident <- unlist(lapply(scseqs, `[[`, 'orig.ident'))
  mnn.out$orig.cluster <- unlist(lapply(scseqs, `[[`, 'cluster'))

  # store merged (batch normalized) for DE
  SummarizedExperiment::assay(mnn.out, 'logcounts') <- SummarizedExperiment::assay(combined, 'merged')

  # need corrected as.matrix for as.Seurat before plots
  SingleCellExperiment::reducedDim(mnn.out, 'corrected') <- as.matrix(SingleCellExperiment::reducedDim(mnn.out, 'corrected'))

  return(mnn.out)
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


#' Get genes that are ambient in at least one test and control sample
#'
#' @param scseqs List of \code{SingleCellExperiment} objects.
#'
#' @return List with test and control ambient genes
#' @export
#' @keywords internal
get_integrated_ambient <- function(scseqs) {

  # datasets that are test samples
  is.test <- sapply(scseqs, function(x) levels(x$orig.ident) == 'test')

  # genes that are ambient in at least one test sample
  ambient.test <- lapply(scseqs[is.test], function(x) SingleCellExperiment::rowData(x))
  ambient.test <- lapply(ambient.test, function(x) row.names(x)[x$out_ambient])
  ambient.test <- unique(unlist(ambient.test))

  # genes that are ambient in at least one ctrl sample
  ambient.ctrl <- lapply(scseqs[!is.test], function(x) SingleCellExperiment::rowData(x))
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

  genes <- row.names(combined)
  SummarizedExperiment::rowData(combined)$test_ambient <- genes %in% ambient$test
  SummarizedExperiment::rowData(combined)$ctrl_ambient <- genes %in% ambient$ctrl

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
exist_clusters <- function(scseq) {
  if (class(scseq) == 'SingleCellExperiment') {
    exist_clusters <- length(unique(scseq$cluster)) > 1

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

