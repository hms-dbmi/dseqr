library(Seurat)
library(drugseqr)

# bone marrow analysis -----
data_dir <- 'data-raw/single-cell/example-data/10X-Schulert-PID200673-20190218-3v3hg/fastqs'

# check whitelisting
scseq <- load_scseq(data_dir, project = 'sjia')
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

# get markers and explore
markers <- get_scseq_markers(scseq)
anal <- list(scseq = scseq, markers = markers)
# saveRDS(anal, 'data-raw/single-cell/example-data/bm_anal.rds')


anal <- readRDS('data-raw/single-cell/example-data/bm_anal.rds')
explore_scseq_clusters(anal$scseq, anal$markers)


# Normal Annotations
# put in order that want reports to appear
cluster_markers <- list('NK-cells'    = c('GZMA', 'PRF1', 'CD247'),
                        'Monocytes'   = c('NCF2', 'CLEC7A', 'C5AR1', 'AIF1'),
                        'B-cells'     = c('CD19', 'MS4A1', 'IGHD', 'BANK1'))

# rename based on identification and save reports
# need to be in same order as clusters
cell_ids <- c('NK-cells', 'Monocytes', 'B-cells')


levels(anal$scseq$seurat_clusters) <- cell_ids
Idents(anal$scseq) <- 'seurat_clusters'
save_scseq_reports(anal$scseq, markers = cluster_markers, pt.size = 3,
                   fname = 'data-raw/single-cell/example-data/bm_sjia_markers.pdf')
