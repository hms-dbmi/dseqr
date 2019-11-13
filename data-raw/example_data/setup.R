#create example app
drugseqr::init_drugseqr('example', local_dir = 'data-raw/patient_data')

# data for example app was downloaded into example app folder:

# single-cell: cellranger files from GSM2560248, GSM2560249
# Bulk: used GEOfastq to get WT treated and untreated from GSE128113
data_dir <- 'data-raw/patient_data/example/bulk'
gse_name <- 'GSE128113'
srp_meta <- GEOfastq::get_srp_meta(gse_name, data_dir)
srp_meta <- srp_meta[grepl('^WT', srp_meta$title), ]
GEOfastq::get_fastqs(gse_name, srp_meta, data_dir)

# ran quantification and differential expression through app
library(drugseqr)
app_dir <- 'inst/app'
data_dir <- 'data-raw/patient_data/example'
pert_query_dir <- 'data-raw/drug_gene_queries/data'
pert_signature_dir <- 'data-raw/drug_es/signatures'
run_drugseqr(data_dir, app_dir, pert_query_dir, pert_signature_dir, test_data = FALSE)

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


# sync to s3 after integration, labeling, etc (whatever want pre-done for app)
# tar cvzf example_data.tar.gz --exclude='*.fastq.gz' --exclude='kallisto*/*' --directory=/home/alex/Documents/Batcave/zaklab/drugseqr/data-raw/patient_data example
# aws s3 cp example_data.tar.gz s3://drugseqr/example_data.tar.gz
