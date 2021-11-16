
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
    is.integrated <- length(unique(sce$orig.ident)) >1

    if (!is.integrated) {
        sce$batch <- dataset_name

    } else if (is.integrated) {
        # treat orig.ident as sample identifier
        sce$batch <- sce$orig.ident
        if (length(unique(sce$batch)) < 2) stop("indicate samples in 'orig.ident'")

        # get HVGs and transfer reductions from integrated assay if present
        alt.names <- SingleCellExperiment::altExpNames(sce)
        if (!'integrated' %in% alt.names) {
            hvgs <- sdata[['RNA']]@var.features

        } else {
            hvgs <- sdata[['integrated']]@var.features

            SingleCellExperiment::reducedDims(sce) <-
            SingleCellExperiment::reducedDims(SingleCellExperiment::altExp(sce, 'integrated'))
        }

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

        rdata <- SummarizedExperiment::rowData(sce)

        rdata$bio <- 0
        rdata[hvgs, 'bio'] <- rev(seq_along(hvgs))
        rdata$hvg <- row.names(sce) %in% hvgs
        SummarizedExperiment::rowData(sce) <- rdata
    }

    # transfer clusters
    # prefer ident as may have annotations
    ident_is_cluster <- identical(as.numeric(sce$ident),
                                  as.numeric(sce$seurat_clusters))

    if (ident_is_cluster) sce$cluster <-sce$ident
    else sce$cluster <- sce$seurat_clusters

    return(sce)
}
