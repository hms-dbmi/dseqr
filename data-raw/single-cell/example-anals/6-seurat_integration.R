library(Seurat)
library(drugseqr)

# import alevin quants
ctrl_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'
test_dir <- 'data-raw/single-cell/example-data/Run2643-10X-Lung/10X_FID12518_Diseased_3hg'

ctrl_scseq <- load_scseq(ctrl_dir, project = 'ctrl')
test_scseq <- load_scseq(test_dir, project = 'test')

# subset by alevin whitelist norm/stabilize using good cell only
ctrl_scseq <- ctrl_scseq[, ctrl_scseq$whitelist]
ctrl_scseq <- preprocess_scseq(ctrl_scseq)

test_scseq <- test_scseq[, test_scseq$whitelist]
test_scseq <- preprocess_scseq(test_scseq)

# perform integration
anchors <- FindIntegrationAnchors(object.list = list(ctrl_scseq, test_scseq))
                                  # anchor.features = intersect(row.names(ctrl_scseq), row.names(test_scseq)))
combined <- IntegrateData(anchors)
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)

# get clusters/tsne/markers
combined <- add_scseq_clusters(combined)
combined <- run_tsne(combined)
markers <- get_scseq_markers(combined)

# explore
explore_scseq_clusters(combined, markers)

combined$cluster.group <- paste(Idents(combined), combined$orig.ident, sep = "_")
Idents(combined) <- "cluster.group"
mono.sjia <- FindMarkers(combined, ident.1 = "0_test", ident.2 = "0_ctrl")
head(mono.sjia, n = 15)

FeaturePlot(combined, features = c("APOBEC3A", "MT1E", "CCL3L1"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

mono.sjia$avg_logFC

