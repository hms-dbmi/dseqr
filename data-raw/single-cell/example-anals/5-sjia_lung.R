library(Seurat)
library(drugseqr)


# control lung analysis ----
ctrl_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'
ctrl_scseq <- load_scseq(ctrl_dir, project = 'ctrl')
ctrl_scseq <- ctrl_scseq[, ctrl_scseq$whitelist]
ctrl_scseq <- preprocess_scseq(ctrl_scseq)

# add clusters and run tsne
ctrl_scseq <- add_scseq_clusters(ctrl_scseq, resolution = 1.4)
ctrl_scseq <- run_umap(ctrl_scseq)

# get markers
ctrl_markers <- get_scseq_markers(ctrl_scseq)
ctrl_anal <- list(scseq = ctrl_scseq, markers = ctrl_markers)
saveRDS(ctrl_anal, 'data-raw/single-cell/example-data/ctrl_lung_anal.rds')


ctrl_anal <- readRDS('data-raw/single-cell/example-data/ctrl_lung_anal.rds')
explore_scseq_clusters(ctrl_anal$scseq, ctrl_anal$markers)


# Normal Annotations
# put in order that want reports to appear
ctrl_markers <- list('Macrophages#1'         = c('PTAFR', 'CD68', 'APOC1'),
                     'Macrophages#2'         = c('FGL2', 'CFP', 'MPEG1'),
                     'T-cells'               = c('CD3E', 'CD3D', 'IL2RG', 'TRBC2'),
                     'Alveolar Epithelium'   = c('NKX2-1', 'SFTA2', 'EMP2', 'KRT19'),
                     'RBCs'                  = c('HBA1', 'HBA2', 'HBB'))

# rename based on identification and save reports
# need to be in same order as clusters
ctrl_anal$annot <- c('T-cells',
                     'RBCs#1',
                     'Alveolar Epithelium',
                     'RBCs#2',
                     'Macrophages#1',
                     'Macrophages#2',
                     '6')

saveRDS(ctrl_anal, 'data-raw/single-cell/example-data/sjia_lung_healthy.rds')

ctrl_anal <- readRDS('data-raw/single-cell/example-data/sjia_lung_healthy.rds')
levels(ctrl_anal$scseq$seurat_clusters) <- ctrl_anal$annot
Idents(ctrl_anal$scseq) <- ctrl_anal$scseq$seurat_clusters
names(ctrl_anal$markers) <- ctrl_anal$annot

save_scseq_reports(ctrl_anal$scseq, markers = ctrl_markers, pt.size = 3,
                   fname = 'data-raw/single-cell/example-data/lung_normal_markers.pdf')


# test lung analysis -----
test_dir <- 'data-raw/single-cell/example-data/Run2643-10X-Lung/10X_FID12518_Diseased_3hg'
test_scseq <- load_scseq(test_dir, project = 'test')
test_scseq <- test_scseq[, test_scseq$whitelist]
test_scseq <- preprocess_scseq(test_scseq)

# add clusters and run tsne
test_scseq <- add_scseq_clusters(test_scseq)
test_scseq <- run_umap(test_scseq)
test_scseq <- jitter_umap(test_scseq)
test_scseq <- jitter_umap(test_scseq)

# get markers
test_markers <- get_scseq_markers(test_scseq)
test_anal <- list(scseq = test_scseq, markers = test_markers)
# saveRDS(test_anal, 'data-raw/single-cell/example-data/test_lung_anal.rds')

test_anal <- readRDS('data-raw/single-cell/example-data/test_lung_anal.rds')
explore_scseq_clusters(test_anal$scseq, test_anal$markers, pt.size = 2.5)

# Diseased Annotations
# put in order that want reports to appear
test_markers <- list('Macrophages#1'         = c('CD68', 'MRC1', 'CTSL', 'MARCO'),
                     'Macrophages#2'         = c('FGL2','FCN1', 'S100A9', 'S100A8'),
                     'T-cells'               = c('IL7R', 'CD3E', 'CD3D'),
                     'NK-cells'              = c('GZMH', 'GZMA', 'CCL4'),
                     'Endothelial'           = c('CLEC14A', 'VWF', 'EMCN', 'CLDN5'),
                     'Alveolar Epithelium'   = c('NKX2-1', 'SFTPB', 'SFTPA2', 'SLPI'),
                     'Mast cells'            = c('MS4A2', 'CPA3', 'TPSAB1', 'HPGDS'),
                     'B-cells#1'             = c('IGHG1', 'IGHG3', 'CD79A'),
                     'B-cells#2'             = c('JCHAIN', 'IGHD', 'SLITRK5'),
                     'SMC/Adipocyte'         = c('COL1A2', 'COL3A1', 'COL6A2', 'PCOLCE'),
                     'Respiratory Cilia'     = c('PIFO', 'RSPH1', 'SNTN', 'RSPH4A'))



# rename based on identification and save reports
# need to be in same order as clusters
test_annot <- c('Macrophages#1',
                'T-cells',
                'NK-cells',
                'Alveolar Epithelium',
                'Macrophages#2',
                'Mast cells',
                'B-cells#1',
                'B-cells#2',
                'Endothelial',
                'SMC/Adipocyte',
                'Respiratory Cilia')

test_anal$annot <- test_annot
# saveRDS(test_anal, 'data-raw/single-cell/example-data/test_lung_anal.rds')

test_anal <- readRDS('data-raw/single-cell/example-data/test_lung_anal.rds')
levels(test_anal$scseq$seurat_clusters) <- test_anal$annot
Idents(test_anal$scseq) <- test_anal$scseq$seurat_clusters
names(test_anal$markers) <- test_anal$annot

data_dir <- '/home/alex/Documents/Batcave/zaklab/drugseqr/data-raw/single-cell/example-anals/sjia/'
explore_scseq_clusters(data_dir, pt.size = 2.5)


save_scseq_reports(test_anal$scseq, markers = test_markers, pt.size = 3,
                   fname = 'data-raw/single-cell/example-data/lung_diseased_markers.pdf')


# Integrated analysis  ----

# perform integration
combined <- integrate_scseq(list(ctrl_anal$scseq, test_anal$scseq))

# get clusters/tsne/markers
combined <- add_scseq_clusters(combined)
combined <- run_umap(combined)
combined <- jitter_umap(combined)
markers <- get_scseq_markers(combined)

anal <- list(scseq = combined, markers = markers)
# look at how integration changes cell groupings
explore_scseq_clusters(anal$scseq, anal$markers)

# Combined Annotations
# put in order that want reports to appear
combined_markers <- list('Macrophages#1'         = c('CD68', 'MRC1', 'IFI30', 'GPNMB', 'CTSB', 'OLR1'),
                         'Macrophages#2'         = c('FGL2','P2RY13', 'CFP', 'FCN1'),
                         'T-cells'               = c('IL7R', 'TRAC', 'LTB', 'CXCR4'),
                         'Endothelial'           = c('CLEC14A', 'ECSCR', 'EMCN', 'CLDN5'),
                         'Alveolar Epithelium#1' = c('S100A14', 'C4BPA', 'C16ORF89', 'ALPL'),
                         'Alveolar Epithelium#2' = c('COL1A2', 'COL1A1', 'UPK3B', 'CRYAB', 'FSTL1'),
                         'RBCs'                  = c('ALAS2', 'HBM', 'HBA1'),
                         'Mast cells'            = c('RGS13', 'MS4A2', 'HPGDS'),
                         'NK/T-cells'            = c('CXCR6', 'GZMA', 'CD3D'),
                         'NK-cells'              = c('FGFBP2', 'KLRD1', 'S1PR5', 'CTSW'),
                         'B-cells'               = c('CD79A', 'DERL3', 'MZB1'),
                         'Respiratory Cilia'     = c('DNAH12', 'HYDIN', 'RSPH4A'))


# rename based on identification and save reports
# need to be in same order as clusters
combined_annot <- c('Macrophages#1',
                    'T-cells',
                    'NK/T-cells',
                    'Alveolar Epithelium#1',
                    'Alveolar Epithelium#2',
                    'B-cells',
                    'Macrophages#2',
                    'RBCs',
                    'Mast cells',
                    'Endothelial',
                    'NK-cells',
                    'Respiratory Cilia')

levels(combined$seurat_clusters) <- combined_annot_not_smooth
orig_clusters <- list('ctrl' = ctrl_anal$scseq$seurat_clusters,
                      'test' = test_anal$scseq$seurat_clusters)

save_combined_scseq_reports(combined, markers = combined_markers, orig_clusters = orig_clusters, point_size = 3,
                            fname = 'data-raw/single-cell/example-data/lung_combined_markers.pdf')

# Integrated analysis without smooth muscle cells  ----

# remove smooth muscle cells because they screw up alveolar epithelial cell clustering
not_smooth <- test_anal$scseq$seurat_clusters != 'SMC/Adipocyte'
not_smooth <- test_anal$scseq$seurat_clusters != 7
test_anal$scseq <- test_anal$scseq[, not_smooth]


# perform integration
scseqs <- list(ctrl_anal$scseq, test_anal$scseq)
combined <- integrate_scseq(scseqs)

# get clusters/tsne/markers
combined <- add_scseq_clusters(combined)
combined <- run_tsne(combined)
markers <- get_scseq_markers(combined)

combined_anal <- list(scseq = combined, markers = markers)
saveRDS(combined_anal, 'data-raw/single-cell/example-data/combined_not_smooth.rds')

# look at how integration changes cell groupings
combined_anal <- readRDS('data-raw/single-cell/example-data/combined_not_smooth.rds')
explore_scseq_clusters(combined_anal$scseq, combined_anal$markers)

Idents(combined) <- paste0(Idents(combined), c(as.character(ctrl_anal$scseq$orig.ident), as.character(test_anal$scseq$orig.ident)))


# Combined Annotations
# put in order that want reports to appear
combined_markers <- list('Macrophages#1'         = c('CD68', 'MRC1', 'IFI30', 'GPNMB', 'CTSB', 'OLR1'),
                         'Macrophages#2'         = c('FGL2','P2RY13', 'CFP', 'FCN1'),
                         'T-cells'               = c('IL7R', 'TRAC', 'LTB', 'CXCR4'),
                         'Endothelial'           = c('CLEC14A', 'ECSCR', 'EMCN', 'CLDN5'),
                         'Alveolar Epithelium#1' = c('NKX2-1', 'SFTA2', 'SFTA3', 'SLPI'),
                         'RBCs'                  = c('ALAS2', 'HBM', 'HBA1'),
                         'Mast cells'            = c('RGS13', 'MS4A2', 'HPGDS'),
                         'NK/T-cells'            = c('CXCR6', 'GZMA', 'CD3D'),
                         'NK-cells'              = c('FGFBP2', 'KLRD1', 'S1PR5', 'CTSW'),
                         'B-cells'               = c('CD79A', 'DERL3', 'MZB1'),
                         'Respiratory Cilia'     = c('DNAH12', 'HYDIN', 'RSPH4A'))

# rename based on identification and save reports
# need to be in same order as clusters
# cluster labels if exclude SMCs
combined_annot <- c('Macrophages#1',
                    'T-cells',
                    'Alveolar Epithelium#1',
                    'NK/T-cells',
                    'RBCs',
                    'Macrophages#2',
                    'B-cells',
                    'Mast cells',
                    'NK-cells',
                    'Endothelial',
                    'Respiratory Cilia')

levels(combined$seurat_clusters) <- combined_annot
orig_clusters <- list('ctrl' = ctrl_anal$scseq$seurat_clusters,
                      'test' = test_anal$scseq$seurat_clusters)

save_combined_scseq_reports(combined, markers = combined_markers, orig_clusters = orig_clusters, point_size = 3,
                            fname = 'data-raw/single-cell/example-data/lung_combined_markers_not_smooth.pdf')


saveRDS(combined, 'data-raw/single-cell/example-data/combined_not_smooth.rds')


# compare test vs control groups where there are enough cells -----

# merge.data = TRUE doesn't seem to merge the scale.data as expected
combined <- merge(ctrl_scseq, test_scseq, add.cell.ids = c("ctrl", "test"), merge.data = TRUE)

# do so manually
rows <- intersect(row.names(ctrl_scseq@assays$SCT@scale.data), row.names(test_scseq@assays$SCT@scale.data))
combined@assays$SCT@scale.data <- cbind(ctrl_scseq@assays$SCT@scale.data[rows, ], test_scseq@assays$SCT@scale.data[rows, ])
colnames(combined@assays$SCT@scale.data) <- colnames(combined)

combined <- add_scseq_clusters(combined)

combined_markers <- FindAllMarkers(combined, assay = 'SCT', slot = 'scale.data')



combined <- readRDS('data-raw/single-cell/example-data/combined_not_smooth.rds')

# remove subgroup distinctions for macrophages
levels(combined$seurat_clusters) <- gsub('#\\d$', '', levels(combined$seurat_clusters))

# make idents group aware
Idents(combined) <- factor(paste(combined$seurat_clusters, combined$orig.ident, sep='_'))

# get markers for all test vs control clusters
cell_groups <- levels(combined$seurat_clusters)
test_markers <- list()

for (cell_group in cell_groups) {

  test_markers[[cell_group]] <- tryCatch(get_scseq_markers(combined,
                                                           ident.1 = paste0(cell_group, '_test'),
                                                           ident.2 = paste0(cell_group, '_ctrl'))[[1]],
                                         error = function(e) return(NULL))
}

explore_scseq_clusters(combined, test_markers)














