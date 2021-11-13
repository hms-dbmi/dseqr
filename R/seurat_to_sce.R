# library(Seurat)
# library(SeuratData)
# library(patchwork)
#
# InstallData("ifnb")
#
# # load dataset
# LoadData("ifnb")
#
# # split the dataset into a list of two seurat objects (stim and CTRL)
# ifnb.list <- SplitObject(ifnb, split.by = "stim")
#
# # normalize and identify variable features for each dataset independently
# ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
#     x <- NormalizeData(x)
#     x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })
#
# # select features that are repeatedly variable across datasets for integration
# features <- SelectIntegrationFeatures(object.list = ifnb.list)
#
# immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
#
# # this command creates an 'integrated' data assay
# immune.combined <- IntegrateData(anchorset = immune.anchors)
#
# # specify that we will perform downstream analysis on the corrected data note that the
# # original unmodified data still resides in the 'RNA' assay
# DefaultAssay(immune.combined) <- "integrated"
#
# # Run the standard workflow for visualization and clustering
# immune.combined <- ScaleData(immune.combined, verbose = FALSE)
# immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
# immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
# immune.combined <- FindClusters(immune.combined, resolution = 0.5)
#
#
# # Visualization
# p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
# p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
# p1 + p2
#
# DimPlot(immune.combined, reduction = "umap", split.by = "stim")
#
# # For performing differential expression after integration, we switch back to the original
# # data
# DefaultAssay(immune.combined) <- "RNA"
# nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
# head(nk.markers)
#
# FeaturePlot(immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A",
#                                           "CCL2", "PPBP"), min.cutoff = "q9")
#
#
# immune.combined <- RenameIdents(immune.combined, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T",
#                                 `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "NK", `7` = "T activated", `8` = "DC", `9` = "B Activated",
#                                 `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets", `14` = "HSPC")
# DimPlot(immune.combined, label = TRUE)
#
# # SCT integration
# LoadData("ifnb")
# ifnb.list <- SplitObject(ifnb, split.by = "stim")
# ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
# features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
# ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
#
# immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
#                                          anchor.features = features)
# immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
#
# immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
# immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)
#
#
# p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "stim")
# p2 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "seurat_annotations", label = TRUE,
#               repel = TRUE)
# p1 + p2

# immune.combined.sct$seurat_clusters <- factor(immune.combined.sct$seurat_annotations)
# immune.combined.sct <-NormalizeData(immune.combined.sct, assay = 'RNA')

seurat_to_sce <- function(sdata, dataset_name) {

    # get the counts and logcounts from RNA assay
    # normalize if need to
    if (identical(sdata[['RNA']]@counts,
                  sdata[['RNA']]@data)) {

        # NOTE: identical if normalize each sample or all together
        sdata <- Seurat::NormalizeData(sdata, assay = 'RNA')
    }

    sce <- Seurat::as.SingleCellExperiment(sdata, assay = 'RNA')
    is.integrated <- 'integrated' %in% Seurat::Assays(sdata)

    if (!is.integrated) {
        sce$batch <- dataset_name

    } else if (is.integrated) {
        # treat orig.ident as sample identifier
        sce$batch <- sce$orig.ident
        if (length(unique(sce$batch)) < 2) stop("indicate samples in 'orig.ident'")

        # get corrected and UMAP/TSNE reductions from integrated assay
        sce.int <- Seurat::as.SingleCellExperiment(sdata, assay = 'integrated')

        reducedDims(sce) <- reducedDims(sce.int)

        # Seurat integrated 'PCA' or 'HARMONY' equivalent to SingleCellExperiment 'corrected'
        # (used to run UMAP/TSNE)
        red.names <- reducedDimNames(sce)

        # prefer harmony over PCA as 'corrected'
        cor.ind <- lapply(c('HARMONY', 'PCA'), grep, red.names)
        cor.ind <- cor.ind[sapply(cor.ind, length) > 0]
        cor.ind <- unlist(cor.ind)[1]

        if (length(cor.ind)) {
            red.names[cor.ind] <- 'corrected'
            reducedDimNames(sce) <- red.names
        }

        # Seurat HVGs are ordered - use to mock biological variance component
        hvgs <- sdata[['integrated']]@var.features
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
