library(drugseqr)
library(Seurat)

run_scseq_anal <- function(data_dir, project, resolution = 0.8) {
    scseq.full <- load_scseq(data_dir, project = project)

    # qc plots
    # hist_scseq_whitelist(scseq.full)
    # run sctransform
    # scseq.full <- preprocess_scseq(scseq.full)
    # tsne_scseq_whitelist(scseq.full)

    # subset by alevin whitelist
    scseq <- scseq.full[, scseq.full$whitelist]
    scseq <- preprocess_scseq(scseq)

    # add clusters
    scseq <- add_scseq_clusters(scseq, resolution = resolution)

    # run tnse and explore clusters
    scseq <- run_tsne(scseq)
    markers <- get_scseq_markers(scseq)

    return(list(scseq=scseq, markers=markers))
}

# control lung analysis ----
ctrl_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'
ctrl_anal <- run_scseq_anal(ctrl_dir, project = 'ctrl', resolution = 2.1)
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
ctrl_annot <- c('T-cells', 'RBCs', 'Alveolar Epithelium#1',
                'RBCs', 'Macrophages#1', 'Macrophages#2',
                'Endothelial', 'Alveolar Epithelium#2')
levels(ctrl_anal$scseq$seurat_clusters) <- ctrl_annot
save_scseq_reports(ctrl_anal$scseq, markers = ctrl_markers, point_size = 3,
                   fname = 'data-raw/single-cell/example-data/lung_normal_markers.pdf')


# test lung analysis -----
test_dir <- 'data-raw/single-cell/example-data/Run2643-10X-Lung/10X_FID12518_Diseased_3hg'
test_anal <- run_scseq_anal(test_dir, 'test')
explore_scseq_clusters(test_anal$scseq, test_anal$markers)

# Diseased Annotations
# put in order that want reports to appear
test_markers <- list('Macrophages#1'         = c('CD68', 'MRC1', 'IFI30', 'LST1'),
                     'Macrophages#2'         = c('FGL2','FCN1', 'S100A12', 'CFP'),
                     'T-cells'               = c('IL7R', 'CD3E', 'CD3D'),
                     'Endothelial'           = c('CLEC14A', 'ECSCR', 'EMCN', 'CLDN5'),
                     'Alveolar Epithelium#1' = c('NKX2-1', 'SFTA2', 'SFTA3', 'SFTPD'),
                     'NK-cells'              = c('GZMH', 'GZMA', 'CCL4'),
                     'Mast cells'            = c('MS4A2', 'SLC18A2', 'TPSD1', 'HPGDS'),
                     'B-cells'               = c('CD79A', 'FCRL5', 'MZB1'),
                     'Smooth Muscle'         = c('COL1A2', 'COL3A1', 'COL6A2', 'PCOLCE'),
                     'Respiratory Cilia'     = c('DNAH12', 'HYDIN', 'RSPH4A'))



# rename based on identification and save reports
# need to be in same order as clusters
test_annot <- c('Macrophages#1', 'T-cells', 'NK-cells', 'Alveolar Epithelium#1',
                'B-cells', 'Mast cells', 'Endothelial', 'Smooth Muscle',
                'Macrophages#2', 'Respiratory Cilia')

levels(test_anal$scseq$seurat_clusters) <- test_annot
save_scseq_reports(test_anal$scseq, markers = test_markers, point_size = 3,
                   fname = 'data-raw/single-cell/example-data/lung_diseased_markers.pdf')


# Integrated analysis ----

not_smooth <- test_anal$scseq$seurat_clusters != 'Smooth Muscle'

# perform integration
anchors <- FindIntegrationAnchors(object.list = list(ctrl_anal$scseq, test_anal$scseq))
# anchors <- FindIntegrationAnchors(object.list = list(ctrl_anal$scseq, test_anal$scseq[, not_smooth]))
# anchor.features = intersect(row.names(ctrl_scseq), row.names(test_scseq)))
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
                         'Alveolar Epithelium#2' = c('UPK3B', 'CRYAB', 'COL1A2', 'FSTL1', 'COL1A1'),
                         'RBCs'                  = c('ALAS2', 'HBM', 'HBA1'),
                         'Mast cells'            = c('RGS13', 'MS4A2', 'HPGDS'),
                         'NK/T-cells'            = c('CXCR6', 'GZMA', 'CD3D'),
                         'NK-cells'              = c('FGFBP2', 'KLRD1', 'S1PR5', 'CTSW'),
                         'B-cells'               = c('CD79A', 'DERL3', 'MZB1'),
                         'Respiratory Cilia'     = c('DNAH12', 'HYDIN', 'RSPH4A'))

# rename based on identification and save reports
# need to be in same order as clusters
combined_annot <- c('Macrophages#1', 'T-cells', 'NK/T-cells', 'Alveolar Epithelium#1', 'Alveolar Epithelium#2',
                    'B-cells', 'Macrophages#2', 'RBCs', 'Mast cells', 'Endothelial', 'NK-cells', 'Respiratory Cilia')

levels(combined$seurat_clusters) <- combined_annot
orig_clusters <- list('ctrl' = ctrl_anal$scseq$seurat_clusters,
                      'test' = test_anal$scseq$seurat_clusters)

save_combined_scseq_reports(combined, markers = combined_markers, orig_clusters = orig_clusters, point_size = 3,
                            fname = 'data-raw/single-cell/example-data/lung_combined_markers.pdf')


