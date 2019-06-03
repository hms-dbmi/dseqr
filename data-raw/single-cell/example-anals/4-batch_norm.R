# http://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/batch.html

# import alevin quants
ctrl_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'
test_dir <- 'data-raw/single-cell/example-data/Run2643-10X-Lung/10X_FID12518_Diseased_3hg'

ctrl_sce <- load_scseq(ctrl_dir)
test_sce <- load_scseq(test_dir)

# subset by alevin whitelist norm/stabilize using good cell only
ctrl_sce <- ctrl_sce[, ctrl_sce$whitelist]
ctrl_sce <- norm_scseq(ctrl_sce)

test_sce <- test_sce[, test_sce$whitelist]
test_sce <- norm_scseq(test_sce)

rescaled <- batchelor::multiBatchNorm(
  test_sce,
  ctrl_sce
)

sce <- cbind(rescaled[[1]], rescaled[[2]])
sce$group <- c(rep('test', ncol(test_sce)), rep('ctrl', ncol(ctrl_sce)))
sce <- stabilize_scseq(sce)


set.seed(1000)
sce <- scater::runTSNE(sce, use_dimred="PCA")
scater::plotTSNE(sce, colour_by = "group")


set.seed(1000)
unc.ctrl <- logcounts(sce)[, sce$group == 'ctrl']
unc.test <- logcounts(sce)[, sce$group == 'test']

mnn.out <- batchelor::fastMNN(
  ctrl=unc.ctrl, test=unc.test,
  k=20, d=50, BSPARAM=BiocSingular::IrlbaParam(deferred=TRUE)
)
mnn.out

# Corrected.
set.seed(1000)
sce.cor <- runTSNE(mnn.out, use_dimred="corrected")
plotTSNE(sce.cor, colour_by="batch", point_alpha = ifelse(sce.cor$group == 'ctrl', 0.2, 1)) + ggtitle("Corrected")

sce.cor <- add_scseq_clusters(sce.cor, "corrected")
markers <- scran::findMarkers(sce.cor, clusters=sce.cor$cluster, direction="up", assay.type="reconstructed")

plotTSNE(sce.cor, by_exprs_values="reconstructed", colour_by="NOC2L")
explore_scseq_clusters(sce.cor, markers, "reconstructed")
