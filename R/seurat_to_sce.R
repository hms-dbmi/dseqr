
seurat_to_sce <- function(sdata, dataset_name) {

  # set to RNA to get the counts and logcounts from RNA assay
  Seurat::DefaultAssay(sdata) <- 'RNA'

  # join layers if v5
  if (methods::is(sdata[['RNA']], 'Assay5')) {
    sdata[['RNA']] <- SeuratObject::JoinLayers(sdata[['RNA']])
  }

  # normalize if need to
  if (identical(sdata[['RNA']]$counts,
                sdata[['RNA']]$data)) {

    # NOTE: identical if normalize each sample or all together
    sdata <- Seurat::NormalizeData(sdata, assay = 'RNA')
  }

  # remove layers with fewer ncols (SCT integration bug?)
  rna_assay <- sdata[['RNA']]
  layers <- SeuratObject::Layers(rna_assay)
  ncol <- sapply(layers, function(x) ncol(SeuratObject::LayerData(rna_assay, x)))
  discard <- names(ncol)[ncol < ncol(rna_assay)]

  if (length(discard)) {
    for (layer in discard)
      SeuratObject::LayerData(rna_assay, layer) <- NULL

    # get normalized and scaled data
    sdata[['RNA']] <- rna_assay
    sdata <- Seurat::NormalizeData(sdata)
    sdata <- Seurat::ScaleData(sdata)
  }


  # meta.features need to have same row.names as assay
  # assays need to have same nrows

  for (assay_name in Seurat::Assays(sdata)) {
    assay <- sdata[[assay_name]]

    meta_name <- 'meta.features'
    if(methods::is(assay, 'Assay5')) {
      assay <- SeuratObject::JoinLayers(assay)
      meta_name <- 'meta.data'
    }
    methods::slot(assay, meta_name) <- methods::slot(assay, meta_name)[row.names(assay), ]
    sdata[[assay_name]] <- assay
  }

  # transfer QC features
  sdata@meta.data$mito_percent <- sdata@meta.data$percent.mt
  sdata@meta.data$ribo_percent <- sdata@meta.data$percent.ribo

  meta_cols <- colnames(sdata@meta.data)
  ncount_cols <- c('nCount_RNA', 'nCount_SCT')
  nfeatr_cols <- c('nFeature_RNA', 'nFeature_SCT')
  ncount_col <- ncount_cols[(ncount_cols %in% meta_cols)][1]
  nfeatr_col <- nfeatr_cols[(nfeatr_cols %in% meta_cols)][1]

  if (length(ncount_col))
    sdata@meta.data$log10_sum <- log10(sdata@meta.data[[ncount_col]])

  if (length(nfeatr_col))
    sdata@meta.data$log10_detected <- log10(sdata@meta.data[[nfeatr_col]])

  sce <- Seurat::as.SingleCellExperiment(sdata)

  # transfer reference info
  sce@metadata$ref_name <- sdata@misc$ref_name
  sce@metadata$resoln <- sdata@misc$resoln

  # transfer samples
  samples <- sce$sample
  if (is.null(samples)) samples <- sce$orig.ident
  have.samples <- length(unique(samples)) > 1

  alt.names <- SingleCellExperiment::altExpNames(sce)
  have.integrated <- 'integrated' %in% alt.names

  if (have.samples) sce$batch <- samples
  else sce$batch <- dataset_name

  # transfer clusters
  sce$cluster <- unname(Seurat::Idents(sdata))

  # default HVGs unless integrated assay
  hvgs <- SeuratObject::VariableFeatures(sdata, assay = 'RNA')

  is.integrated <- have.samples | have.integrated
  if (is.integrated) {

    # get HVGs and transfer reductions from integrated assay if present
    if (have.integrated) {
      hvgs <- SeuratObject::VariableFeatures(sdata, assay = 'integrated')
      SingleCellExperiment::reducedDims(sce) <-
        SingleCellExperiment::reducedDims(SingleCellExperiment::altExp(sce, 'integrated'))
    }

    SingleCellExperiment::altExps(sce) <- NULL
    rm(sdata); gc()

    # Seurat integrated 'PCA' or 'HARMONY' equivalent to SingleCellExperiment 'corrected'
    # (used to run UMAP/TSNE)
    red.names <- SingleCellExperiment::reducedDimNames(sce)

    # prefer harmony over PCA as 'corrected'
    cor.ind <- lapply(c('HARMONY', 'PCA'), grep, red.names)
    cor.ind <- cor.ind[sapply(cor.ind, length) > 0]
    cor.ind <- unlist(cor.ind)[1]

    if (length(cor.ind)) {
      red.names[cor.ind] <- 'corrected'
      SingleCellExperiment::reducedDimNames(sce) <- red.names
    }
  }


  # HVG order indicates biological variance
  rdata <- SummarizedExperiment::rowData(sce)

  rdata$bio <- 0
  rdata[hvgs, 'bio'] <- rev(seq_along(hvgs))
  rdata$hvg <- row.names(sce) %in% hvgs
  SummarizedExperiment::rowData(sce) <- rdata

  # remove reductions that won't use
  reds <- SingleCellExperiment::reducedDims(sce)
  reds.keep <- c('corrected', 'PCA', 'UMAP', 'TSNE')

  # no PCA if have corrected
  if ('corrected' %in% names(reds)) reds.keep <- reds.keep[-2]
  reds <- reds[names(reds) %in% reds.keep]

  SingleCellExperiment::reducedDims(sce) <- reds

  # remove cdata that won't use
  keep_cols <- c('cluster', 'batch', const$features$qc)
  is.ref <- !is.null(sce@metadata$ref_name)
  if (is.ref) {
    cluster_cols <- get_ref_cols(colnames(sce@colData), type = 'cluster')
    keep_cols <- c(keep_cols, cluster_cols)
  }

  sce@colData <- sce@colData[, colnames(sce@colData) %in% keep_cols]

  return(sce)
}
