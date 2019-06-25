library(Seurat)
library(drugseqr)

# sample1 -----
data_dir <- 'data-raw/single-cell/example-data/10X-Schulert-PID200673-20190218-3v3hg/fastqs'

# check whitelisting
scseq <- load_scseq(data_dir, project = 'sjia_blood1')
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
saveRDS(anal, 'data-raw/single-cell/example-data/sjia_blood1.rds')


anal <- readRDS('data-raw/single-cell/example-data/sjia_blood1.rds')
explore_scseq_clusters(anal1$scseq, anal1$markers)


# Normal Annotations
# put in order that want reports to appear
cluster_markers <- list('NK-cells'    = c('GZMA', 'PRF1', 'CD247'),
                        'Monocytes'   = c('NCF2', 'CLEC7A', 'C5AR1', 'AIF1'),
                        'B-cells'     = c('CD19', 'MS4A1', 'IGHD', 'BANK1'))

# rename based on identification and save reports
# need to be in same order as clusters
cell_ids <- c('NK-cells', 'Monocytes', 'B-cells')


levels(anal$scseq$seurat_clusters) <- cell_ids
save_scseq_reports(anal$scseq, markers = cluster_markers, pt.size = 3,
                   fname = 'data-raw/single-cell/example-data/sjia_blood1_markers.pdf')

# sample 2 ------
data_dir <- 'data-raw/single-cell/example-data/10X-Schulert-PID200396-20190218-3v3hg/fastqs'

# check whitelisting
scseq <- load_scseq(data_dir, project = 'sjia_blood2')
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

# cells that cluster with B-cells but mainly because of Ribosome genes
select.cells <- Seurat::CellSelector(DimPlot(scseq))
Idents(scseq, cells = select.cells) <- "3"
scseq@active.ident <- scseq$seurat_clusters <- factor(Idents(scseq), levels = sort(levels(Idents(scseq))))

# get markers and explore
markers <- get_scseq_markers(scseq)
anal <- list(scseq = scseq, markers = markers)
saveRDS(anal, 'data-raw/single-cell/example-data/sjia_blood2.rds')

anal <- readRDS('data-raw/single-cell/example-data/sjia_blood2.rds')
explore_scseq_clusters(anal$scseq, anal$markers)


# Normal Annotations
# put in order that want reports to appear
cluster_markers <- list('NK-cells'    = c('GNLY', 'NKG7', 'GZMB'),
                        'Monocytes'   = c('S100A12', 'PLBD1', 'MS4A6A', 'MNDA', 'AIF1'),
                        'B-cells'     = c('LTB', 'MS4A1', 'IGHD', 'BANK1'),
                        'T-cells'     = c('IL7R', 'CD3D', 'CD3G', 'TRAC'))

# rename based on identification and save reports
# need to be in same order as clusters
cell_ids <- c('NK-cells', 'Monocytes', 'B-cells', 'T-cells')


levels(anal$scseq$seurat_clusters) <- cell_ids
save_scseq_reports(anal$scseq, markers = cluster_markers, pt.size = 3,
                   fname = 'data-raw/single-cell/example-data/sjia_blood2_markers.pdf')


# integration ----

anal1 <- readRDS('data-raw/single-cell/example-data/sjia_blood1.rds')
anal2 <- readRDS('data-raw/single-cell/example-data/sjia_blood2.rds')

combined <- integrate_scseq(list(anal1$scseq, anal2$scseq))

# get clusters/tsne/markers
combined <- add_scseq_clusters(combined)
combined <- run_umap(combined)
markers <- get_scseq_markers(combined)

anal <- list(scseq = combined, markers = markers)
saveRDS(anal, 'data-raw/single-cell/example-data/sjia_blood1&2.rds')

# look at how integration changes cell groupings
explore_scseq_clusters(anal$scseq, anal$markers)

cluster_markers <- list('NK-cells'    = c('GZMA', 'PRF1', 'CD247'),
                        'Monocytes'   = c('NCF2', 'CLEC7A', 'C5AR1', 'AIF1'),
                        'B-cells'     = c('BANK1', 'MS4A1', 'IGHD', 'VPREB3'))

# rename based on identification and save reports
# need to be in same order as clusters
annot <- c('NK-cells', 'Monocytes', 'B-cells')
anal$annot <- annot
saveRDS(anal, 'data-raw/single-cell/example-data/sjia_blood1&2.rds')

levels(combined$seurat_clusters) <- cell_ids

orig_clusters <- list('sjia_blood1' = anal1$scseq$seurat_clusters,
                      'sjia_blood2' = anal2$scseq$seurat_clusters)

save_combined_scseq_reports(combined, markers = cluster_markers, orig_clusters = orig_clusters, pt.size = 3,
                            fname = 'data-raw/single-cell/example-data/blood_combined.pdf')
