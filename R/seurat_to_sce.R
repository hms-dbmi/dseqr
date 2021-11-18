
seurat_to_sce <- function(sdata, dataset_name) {

    # get the counts and logcounts from RNA assay
    # normalize if need to
    if (identical(sdata[['RNA']]@counts,
                  sdata[['RNA']]@data)) {

        # NOTE: identical if normalize each sample or all together
        sdata <- Seurat::NormalizeData(sdata, assay = 'RNA')
    }

    # transfer QC features
    sdata@meta.data$mito_percent <- sdata@meta.data$percent.mt
    sdata@meta.data$ribo_percent <- sdata@meta.data$percent.ribo
    sdata@meta.data$log10_sum <- log10(sdata@meta.data$nCount_RNA)
    sdata@meta.data$log10_detected <- log10(sdata$nFeature_RNA)

    sce <- Seurat::as.SingleCellExperiment(sdata)

    # transfer samples
    samples <- sce$sample
    if (is.null(samples)) samples <- sce$orig.ident
    is.integrated <- length(unique(samples)) >1

    # transfer clusters
    sce$cluster <- unname(Seurat::Idents(sdata))

    # default HVGs unless integrated assay
    hvgs <- sdata[['RNA']]@var.features

    if (!is.integrated) {
        sce$batch <- dataset_name

    } else {
        sce$batch <- samples

        # get HVGs and transfer reductions from integrated assay if present
        alt.names <- SingleCellExperiment::altExpNames(sce)
        if ('integrated' %in% alt.names) {
            hvgs <- sdata[['integrated']]@var.features
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

    return(sce)
}
