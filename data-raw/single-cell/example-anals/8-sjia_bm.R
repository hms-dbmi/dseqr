library(Seurat)
library(drugseqr)

# bone marrow analysis -----
data_dir <- 'data-raw/single-cell/example-data/Run2622-10X-BoneMarrow-STA1/10X_STA-1_3hg'

# check whitelisting
scseq <- load_scseq(data_dir, project = 'sjia_bm')
scseq <- preprocess_scseq(scseq)
scseq <- add_scseq_clusters(scseq)
scseq <- run_umap(scseq)
hist_scseq_whitelist(scseq)
tsne_scseq_whitelist(scseq)

# subset by whitelist and redo analysis
scseq <- scseq[, scseq$whitelist]
scseq <- preprocess_scseq(scseq)
scseq <- add_scseq_clusters(scseq)
scseq <- run_umap(scseq)
scseq <- jitter_umap(scseq)

# get markers and explore
markers <- get_scseq_markers(scseq)
anal <- list(scseq = scseq, markers = markers)
saveRDS(anal, 'data-raw/single-cell/example-data/bm_anal.rds')


anal <- readRDS('data-raw/single-cell/example-data/bm_anal.rds')
explore_scseq_clusters(anal$scseq, anal$markers)


# Normal Annotations
# put in order that want reports to appear
cluster_markers <- list('B-cells#1'       = c('IGHD', 'TNFRSF13C', 'MS4A1'),
                        'B-cells#2'       = c('CD79B', 'BCL7A'),
                        'Pre B-cells'     = c('IGLL1', 'DNTT'),
                        'Dendritic Cells' = c('SCT', 'PTPRS', 'CSF2RB'),
                        'T-cells'         = c('CD3E', 'MAL', 'CD3D'),
                        'NK-cells'        = c('CTSW', 'GZMH', 'CST7'),
                        'Monocytes'       = c('RAB31', 'FCN1', 'HMOX1', 'LST1'),
                        'Erythroid'       = c('SELENBP1', 'TMCC2', 'ALAS2'))

# rename based on identification and save reports
# need to be in same order as clusters
annot <- c('B-cells#1',
           'T-cells',
           'Monocytes',
           'B-cells#2',
           'NK-cells',
           'Pre B-cells',
           'Dendritic Cells',
           'Erythroid')

anal$annot <- annot
saveRDS(anal, 'data-raw/single-cell/example-data/bm_anal.rds')

levels(anal$scseq$seurat_clusters) <- annot
save_scseq_reports(anal$scseq, markers = cluster_markers, pt.size = 3,
                   fname = 'data-raw/single-cell/example-data/bm_sjia_markers.pdf')
