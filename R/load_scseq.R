#' Load kallisto/bustools quantification into a Seurat or SingleCellExperiment object.
#'
#' @param data_dir Directory with raw and kallisto/bustools or CellRanger quantified single-cell RNA-Seq files.
#' @param project String identifying sample.
#' @param type Quantification file type. One of either \code{'kallisto'} or \code{'cellranger'}.
#' @param knee_type Knee type used for cell quality whitelist modeling by \code{get_scseq_whitelist}.
#'
#' @seealso \link{pick_roryk_cutoff} for details of \code{'roryk'} \code{knee_type}.
#' \link[DropletUtils]{barcodeRanks} for \code{'inflection'} and \code{'knee'} \code{knee_type}.
#'

#' @return \code{Seurat} object with whitelist meta data.
#' @export
#'
#' @examples
#'
#' data_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'
#' load_scseq(data_dir)
#'
load_scseq <- function(data_dir, project, type = c('kallisto', 'cellranger'), knee_type = c('inflection', 'roryk', 'knee')) {

  # load counts
  if (type[1] == 'kallisto') {
    data_dir <- file.path(data_dir, 'bus_output')
    counts <- load_kallisto_counts(data_dir)

  } else if (type[1] == 'cellranger') {
    counts <- load_cellranger_counts(data_dir)
  }

  # generate/load whitelist
  whitelist <- get_scseq_whitelist(counts, data_dir, knee_type = knee_type)
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

  # covert to SingleCellExperiment object
  # add ambient metadata for genes
  rowData <- S4Vectors::DataFrame(pct_ambient, out_ambient)
  colData <- S4Vectors::DataFrame(project = rep(project, length(kneelist)),
                       whitelist = kneelist %in% whitelist)

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
  counts <- as(counts, 'dgCMatrix')

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
  preclusters <- scran::quickCluster(scseq, BPPARAM = BiocParallel::DoparParam(), min.size = 5)
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


#' Run TSNE
#'
#' @param sce \code{SingleCellExperiment}
#' @param dimred reducedDim to run TSNE on
#'
#' @return \code{sce} with \code{'TSNE'} \code{reducedDim}
#' @export
run_tsne <- function(sce, dimred = 'PCA') {
  set.seed(1100101001)
  sce <- scater::runTSNE(sce, dimred = dimred, n_dimred = sce@metadata$npcs)
  colnames(SingleCellExperiment::reducedDim(sce, 'TSNE')) <- c('TSNE1', 'TSNE2')
  return(sce)
}


#' Get number of clusters different number of PCs
#'
#' Used to pick number of PCs to retain
#'
#' @param sce \code{SingleCellExperiement}
#'
#' @return result of \code{scran::getClusteredPCs}
#' @export
#'
get_npc_choices <- function(sce, type = 'PCA') {

  # walktrap very slow if too many cells
  cluster_fun <- ifelse(ncol(sce) > 10000,
                        igraph::cluster_louvain,
                        igraph::cluster_walktrap)

  FUN <- function(x, ...) {
    g <- scran::buildSNNGraph(x, ..., transposed = TRUE)
    cluster_fun(g)$membership
  }

  # at most 50 pcs
  pcs <- SingleCellExperiment::reducedDim(sce, type = type)
  pcs <- pcs[,seq_len(min(50, ncol(pcs)))]

  choices <- scran::getClusteredPCs(pcs, FUN = FUN)
  names(choices$clusters) <- choices$n.pcs
  return(choices)
}

#' Cluster SingleCellExperiment
#'
#' @param sce \code{SingleCellExperiment}
#'
#' @return \code{sce} with column \code{cluster} in colData and \code{'npcs'} in metadata
#' @export
add_scseq_clusters <- function(sce) {

  # run PCA on HVGs
  rdata <- SummarizedExperiment::rowData(sce)
  subset_row <- row.names(rdata[rdata$hvg, ])

  set.seed(100)
  sce <- scater::runPCA(sce, subset_row = subset_row)

  # pick number of PCs
  choices <- get_npc_choices(sce)
  npcs <- S4Vectors::metadata(choices)$chosen
  sce@metadata$npcs <- npcs

  # add clusters
  sce$cluster <- factor(choices$clusters[[as.character(npcs)]])
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
  sce <- scater::addPerCellQC(sce,
                              subsets=list(mito = which(row.names(sce) %in% sce@metadata$mrna),
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

  markers <- lapply(markers, as.data.frame)
  ord <- order(as.numeric(names(markers)))
  markers[ord]
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
  no_correct <- function(assay.type) function(...) batchelor::noCorrect(..., assay.type = assay.type)
  combined <- do.call(no_correct('logcounts'), scseqs)


  set.seed(1000101001)
  if (type[1] == 'clusterMNN') {
    mnn.fun <- function(...) batchelor::clusterMNN(
      ...,
      clusters = lapply(scseqs, `[[`, 'cluster'),
      subset.row = hvgs,
      auto.merge = TRUE,
      correct.all = TRUE,
      cos.norm = FALSE)

  } else if (type[1] == 'fastMNN') {
    mnn.fun <- function(...) batchelor::fastMNN(
      ...,
      subset.row = hvgs,
      auto.merge = TRUE,
      correct.all = TRUE,
      cos.norm = FALSE,
      prop.k = 0.05)
  }

  mnn.out <- do.call(mnn.fun, scseqs) ; gc()
  mnn.out$orig.ident <- unlist(lapply(scseqs, `[[`, 'orig.ident'), use.names = FALSE)
  mnn.out$orig.cluster <- unlist(lapply(scseqs, `[[`, 'cluster'), use.names = FALSE)

  # store merged (batch normalized) for DE
  SummarizedExperiment::assay(mnn.out, 'logcounts') <- SummarizedExperiment::assay(combined, 'merged')
  rm(combined); gc()

  # get counts for pseudobulk
  counts <- do.call(no_correct('counts'), scseqs)
  SummarizedExperiment::assay(mnn.out, 'counts') <- SummarizedExperiment::assay(counts, 'merged')
  rm(counts, scseqs); gc()

  # need corrected as.matrix for as.Seurat before plots
  SingleCellExperiment::reducedDim(mnn.out, 'corrected') <- as.matrix(SingleCellExperiment::reducedDim(mnn.out, 'corrected'))

  return(mnn.out)
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
  length(unique(scseq$cluster)) > 1
}

