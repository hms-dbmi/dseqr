#create example app
drugseqr::init_drugseqr('example', local_dir = 'data-raw/patient_data')

# data for example app was downloaded into example app folder:

# single-cell: cellranger files from GSM2560248, GSM2560249
# Bulk: used GEOfastq to get WT treated and untreated from GSE128113
# ran quantification and differential expression through app

# make single cell datasets smaller for example purposes

# for IFNb stimulated
scseq <- readRDS('data-raw/patient_data/example/single-cell/GSM2560249_pbmc_ifnb/scseq.rds')
scseq <- drugseqr::downsample_scseq(scseq, 3000)
scseq[['SCT']]@misc$vst.out$cell_attr <- scseq[['SCT']]@misc$vst.out$cell_attr[colnames(scseq), ]

saveRDS(scseq, 'data-raw/patient_data/example/single-cell/GSM2560249_pbmc_ifnb/scseq.rds')



scseq <- readRDS('data-raw/patient_data/example/single-cell/GSM2560248_pbmc_ctrl/scseq.rds')
scseq <- drugseqr::downsample_scseq(scseq, 3000)
scseq[['SCT']]@misc$vst.out$cell_attr <- scseq[['SCT']]@misc$vst.out$cell_attr[colnames(scseq), ]

saveRDS(scseq, 'data-raw/patient_data/example/single-cell/GSM2560248_pbmc_ctrl/scseq.rds')
