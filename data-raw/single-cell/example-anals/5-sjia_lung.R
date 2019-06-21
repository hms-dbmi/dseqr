library(Seurat)
library(drugseqr)

run_scseq_anal <- function(data_dir, project, resolution = 0.8) {
    scseq.full <- load_scseq(data_dir, project = project)

    # qc plots
    hist_scseq_whitelist(scseq.full)
    # run sctransform
    scseq.full <- preprocess_scseq(scseq.full)
    tsne_scseq_whitelist(scseq.full)

    # subset by alevin whitelist
    scseq <- scseq.full[, scseq.full$whitelist]
    scseq <- preprocess_scseq(scseq)

    # add clusters
    scseq <- add_scseq_clusters(scseq, resolution = resolution)

    # run tnse
    scseq <- run_tsne(scseq)
    return(scseq)
}

# control lung analysis ----
ctrl_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'
ctrl_scseq <- run_scseq_anal(ctrl_dir, project = 'ctrl', resolution = 1.6)
ctrl_markers <- get_scseq_markers(ctrl_scseq)
ctrl_anal <- list(scseq = ctrl_scseq, markers = ctrl_markers)
saveRDS(ctrl_anal, 'data-raw/single-cell/example-data/ctrl_anal.rds')


ctrl_anal <- readRDS('data-raw/single-cell/example-data/ctrl_anal.rds')
explore_scseq_clusters(ctrl_anal$scseq, ctrl_anal$markers)


# Normal Annotations
# put in order that want reports to appear
ctrl_markers <- list('Macrophages#1'         = c('TNFSF13', 'ADAMTSL4', 'TREM1'),
                     'Macrophages#2'         = c('FGL2', 'CFP', 'MPEG1', 'TLR2'),
                     'T-cells'               = c('CD3E', 'CD3D', 'IL2RG', 'TRBC2'),
                     'Endothelial'           = c('ROBO4', 'PODXL', 'RAMP3', 'ACVRL1'),
                     'Alveolar Epithelium#1' = c('NKX2-1', 'SFTA2', 'EMP2', 'KRT19'),
                     'Alveolar Epithelium#2' = c('IRX3', 'SCGB3A1', 'ALPL'),
                     'RBCs'                  = c('HBA1', 'HBA2', 'HBB'))

# rename based on identification and save reports
# need to be in same order as clusters
ctrl_annot <- c('T-cells',
                'RBCs',
                'Alveolar Epithelium#1',
                'RBCs',
                'Macrophages#1',
                'Macrophages#2',
                'Endothelial',
                'Alveolar Epithelium#2')

levels(ctrl_anal$scseq$seurat_clusters) <- ctrl_annot
save_scseq_reports(ctrl_anal$scseq, markers = ctrl_markers, point_size = 3,
                   fname = 'data-raw/single-cell/example-data/lung_normal_markers.pdf')


# test lung analysis -----
test_dir <- 'data-raw/single-cell/example-data/Run2643-10X-Lung/10X_FID12518_Diseased_3hg'
test_scseq <- run_scseq_anal(test_dir, 'test')
test_markers <- get_scseq_markers(test_scseq)
test_anal <- list(scseq = test_scseq, markers = test_markers)
saveRDS(test_anal, 'data-raw/single-cell/example-data/test_anal.rds')

test_anal <- readRDS('data-raw/single-cell/example-data/test_anal.rds')
explore_scseq_clusters(test_anal$scseq, test_anal$markers)

# Diseased Annotations
# put in order that want reports to appear
test_markers <- list('Macrophages#1'         = c('CD68', 'MRC1', 'IFI30', 'LST1'),
                     'Macrophages#2'         = c('FGL2','FCN1', 'S100A9', 'S100A8'),
                     'T-cells'               = c('IL7R', 'CD3E', 'CD3D'),
                     'Endothelial'           = c('CLEC14A', 'VWF', 'EMCN', 'CLDN5'),
                     'Alveolar Epithelium#1' = c('NKX2-1', 'SFTPB', 'SFTPA2', 'SLPI'),
                     'NK-cells'              = c('GZMH', 'GZMA', 'CCL4'),
                     'Mast cells'            = c('MS4A2', 'CPA3', 'TPSAB1', 'HPGDS'),
                     'B-cells'               = c('CD79A', 'IGKC', 'MZB1'),
                     'SMC/Adipocyte'         = c('COL1A2', 'COL3A1', 'COL6A2', 'PCOLCE'),
                     'Respiratory Cilia'     = c('PIFO', 'RSPH1', 'SNTN', 'RSPH4A'))



# rename based on identification and save reports
# need to be in same order as clusters
test_annot <- c('Macrophages#1',
                'T-cells',
                'NK-cells',
                'Alveolar Epithelium#1',
                'B-cells',
                'Mast cells',
                'Endothelial',
                'SMC/Adipocyte',
                'Macrophages#2',
                'Respiratory Cilia')

levels(test_anal$scseq$seurat_clusters) <- test_annot
save_scseq_reports(test_anal$scseq, markers = test_markers, point_size = 3,
                   fname = 'data-raw/single-cell/example-data/lung_diseased_markers.pdf')


# Integrated analysis  ----

# perform integration
anchors <- FindIntegrationAnchors(object.list = list(ctrl_anal$scseq, test_anal$scseq))
combined <- IntegrateData(anchors)
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)

# get clusters/tsne/markers
combined <- add_scseq_clusters(combined)
combined <- run_tsne(combined)
markers <- get_scseq_markers(combined)

# look at how integration changes cell groupings
explore_scseq_clusters(combined, markers)

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














