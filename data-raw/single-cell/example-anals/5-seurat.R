library(drugseqr)
library(Seurat)
library(cowplot)

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
ctrl_markers <- list('Alveolar Epithelium#1' = c('NKX2-1', 'SFTA2', 'EMP2', 'KRT19'),
                     'Alveolar Epithelium#2' = c('IRX3', 'SCGB3A1', 'ALPL'),
                     'T-cells'             = c('CD3E', 'CD3D', 'IL2RG', 'TRBC2'),
                     'RBCs'                = c('HBA1', 'HBA2', 'HBB'),
                     'Macrophages#1'       = c('TNFSF13', 'ADAMTSL4', 'TREM1'),
                     'Macrophages#2'       = c('FGL2', 'CFP', 'MPEG1', 'TLR2'),
                     'Endothelial'         = c('ROBO4', 'PODXL', 'RAMP3', 'ACVRL1'))

# rename based on identification and save reports
ctrl_annot <- c('T-cells', 'RBCs', 'Alveolar Epithelium#1', 'RBCs', 'Macrophages#1', 'Macrophages#2', 'Endothelial', 'Alveolar Epithelium#2')
levels(ctrl_anal$scseq$seurat_clusters) <- ctrl_annot
save_scseq_reports(ctrl_anal$scseq, fname = 'control.pdf', markers = ctrl_markers, point_size = 2)


# test lung analysis -----
test_dir <- 'data-raw/single-cell/example-data/Run2643-10X-Lung/10X_FID12518_Diseased_3hg'
test_anal <- run_scseq_anal(test_dir, 'test')
explore_scseq_clusters(test_anal$scseq, test_anal$markers)

# Diseased Annotations
test_markers <- list('Macrophages#1'       = c('CD68', 'MRC1', 'IFI30', 'LST1'),
                     'T-cells'             = c('IL7R', 'CD3E', 'CD3D'),
                     'NK-cells'            = c('GZMH', 'GZMA', 'CCL4'),
                     'Alveolar Epithelium' = c('NKX2-1', 'SFTA2', 'SFTA3', 'SFTPD'),
                     'B-cells'             = c('CD79A', 'FCRL5', 'MZB1'),
                     'Mast cells'          = c('MS4A2', 'SLC18A2', 'TPSD1', 'HPGDS'),
                     'Endothelial'         = c('CLEC14A', 'ECSCR', 'EMCN', 'CLDN5'),
                     'Smooth Muscle'       = c('COL1A2', 'COL3A1', 'COL6A2', 'PCOLCE'),
                     'Macrophages#2'       = c('FGL2','FCN1', 'S100A12', 'CFP'),
                     'Respiratory Cilia'   = c('DNAH12', 'HYDIN', 'RSPH4A'))



# rename based on identification and save reports
test_annot <- c('Macrophages#1', 'T-cells', 'NK-cells',
                'Alveolar Epithelium', 'B-cells', 'Mast cells',
                'Endothelial', 'Smooth Muscle Cells',
                'Macrophages#2', 'Respiratory Cilia')

levels(test_anal$scseq$seurat_clusters) <- test_annot
save_scseq_reports(test_anal$scseq, fname = 'test.pdf', markers = test_markers, point_size = 3)


# Integrated analysis ----

# label by cell type


# perform integration
anchors <- FindIntegrationAnchors(object.list = list(ctrl_anal$scseq, test_anal$scseq))
# anchor.features = intersect(row.names(ctrl_scseq), row.names(test_scseq)))
combined <- IntegrateData(anchors)
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)

# get clusters/tsne/markers
combined <- add_scseq_clusters(combined)
combined <- run_tsne(combined)
markers <- get_scseq_markers(combined)

# look at how integration changes cell groupings
combined$labels <- c(paste0(as.character(Idents(ctrl_anal$scseq)), '_ctrl'), paste0(as.character(Idents(test_anal$scseq)), '_test'))
explore_scseq_clusters(combined, markers)

# compare test and control
