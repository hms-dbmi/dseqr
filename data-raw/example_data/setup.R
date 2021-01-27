#create example app
drugseqr::init_drugseqr(app_name = 'example', data_dir = 'data-raw/patient_data')

# data for example app was downloaded into example app folder:

# single-cell: cellranger files
# Bulk: used GEOfastq to get progress and healthy from GSE93624
data_dir <- 'data-raw/patient_data/example/bulk'
gse_name <- 'GSE93624'
srp_meta <- GEOfastq::get_srp_meta(gse_name, data_dir)
srp_meta <- srp_meta[grepl('complication: yes', srp_meta$characteristics)
                     | srp_meta$source_name == 'Non-IBD control_Ileal biopsy', ]
GEOfastq::get_fastqs(gse_name, srp_meta, data_dir)


# just eight single-cell samples for example purposes
gsms <- c('GSM3972009', 'GSM3972010', 'GSM3972011', 'GSM3972012', 'GSM3972013', 'GSM3972014', 'GSM3972016', 'GSM3972015')
snums <- c('69', '68', '122', '123', '128', '129', '138', '135')

for (i in seq_along(gsms)) {
    gsm <- gsms[i]
    snum <- snums[i]

    # samples ordered above so that involved then uninvolved
    stype <- ifelse(i %% 2 == 0, 'uninvolved', 'involved')
    gsm_snum <- paste0(gsm, '_', snum, '_')

    data_dir <- paste0('data-raw/patient_data/example/single-cell/', gsm, '_illeal_', stype, snum)
    dir.create(data_dir)
    download.file(paste0('ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3972nnn/', gsm, '/suppl/', gsm_snum, 'barcodes.tsv.gz'), file.path(data_dir, paste0(gsm_snum, 'barcodes.tsv.gz')))
    download.file(paste0('ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3972nnn/', gsm, '/suppl/', gsm_snum, 'genes.tsv.gz'), file.path(data_dir, paste0(gsm_snum, 'genes.tsv.gz')))
    download.file(paste0('ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3972nnn/', gsm, '/suppl/', gsm_snum, 'matrix.mtx.gz'), file.path(data_dir, paste0(gsm_snum, 'matrix.mtx.gz')))
}

# load single-cell datasets
indices_dir <- '/srv/drugseqr/indices'
sc_dir <- 'data-raw/patient_data/example/single-cell'
dataset_names <- list.dirs(sc_dir, recursive = FALSE, full.names = FALSE)

for (i in seq_along(dataset_names)) {
  dataset_name <- dataset_names[i]
  cat('Working on', dataset_name, '...\n')
  fastq_dir <- file.path(sc_dir, dataset_name)

  files <- list.files(fastq_dir)
  keep <- grep('matrix.mtx|barcodes.tsv|genes.tsv', files, value = TRUE)
  files <- files[!files %in% keep]
  unlink(file.path(fastq_dir, files), recursive = TRUE, force = TRUE)


  # exclude QC0 for size
  # load_raw_scseq(paste0(dataset_name, '_QC0'), fastq_dir, sc_dir, indices_dir, recount = FALSE, metrics = NULL, founder = dataset_name)
  drugseqr::load_raw_scseq(paste0(dataset_name, '_QC1'), fastq_dir, sc_dir, indices_dir, recount = FALSE, founder = dataset_name)
}


# ran quantification and differential expression through app
# dataset name: GSE93624_ped_crohns_illeum
app_name <- 'example'
app_dir <- 'inst/app'
data_dir <- 'data-raw/patient_data'
pert_query_dir <- 'data-raw/drug_gene_queries/data'
pert_signature_dir <- 'data-raw/drug_es/signatures'
indices_dir <- '/srv/drugseqr/indices'

drugseqr::run_drugseqr(app_name, data_dir, app_dir, pert_query_dir, pert_signature_dir, indices_dir, port = 3838)


# sync to s3 after integration, labeling, etc (whatever want pre-done for app)
# tar cvzf example_data.tar.gz --exclude='*.fastq.gz' --exclude='kallisto*/*' --exclude '*/barcodes.tsv' --exclude '*/matrix.mtx' --exclude '*/genes.tsv' --directory=/home/alex/patient_data example
# aws s3 cp example_data.tar.gz s3://drugseqr/example_data.tar.gz
