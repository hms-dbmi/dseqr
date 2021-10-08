# data for example app was downloaded into example app folder:

# single-cell: cellranger files
# Bulk: used GEOfastq to get progress and healthy from GSE93624
data_dir <- '/mnt/12tb/Data/dseqr/example/bulk/GSE93624'
gse_name <- 'GSE93624'
meta_path <- file.path(data_dir, 'srp_meta.rds')
# gse_text <- GEOfastq::crawl_gse(gse_name)
# gsm_names <- GEOfastq::extract_gsms(gse_text)
# srp_meta <- GEOfastq::crawl_gsms(gsm_names)
# saveRDS(srp_meta, meta_path)
srp_meta <- readRDS(meta_path)

is.ctrl <- srp_meta$source_name == 'Non-IBD control_Ileal biopsy'
is.test <- grepl('complication: yes', srp_meta$characteristics_ch1.5)

# just two each
srp_meta <- srp_meta[c(which(is.ctrl)[1:2], which(is.test)[1:2]), ]
# srp_meta <- srp_meta[is.ctrl | is.test, ]


GEOfastq::get_fastqs(srp_meta[3:4,], data_dir)


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
indices_dir <- '/srv/dseqr/indices'
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
  dseqr::load_raw_scseq(paste0(dataset_name, '_QC1'), fastq_dir, sc_dir, indices_dir, recount = FALSE, founder = dataset_name)
}


# ran quantification and differential expression through app
# dataset name: GSE93624_ped_crohns_illeum
app_name <- 'example'
app_dir <- 'inst/app'
data_dir <- '~/patient_data'
pert_query_dir <- '~/drug_gene_queries/data'
pert_signature_dir <- '~/drug_es/signatures'
indices_dir <- '/srv/dseqr/indices'

dseqr::run_dseqr(app_name, data_dir, app_dir, pert_query_dir, pert_signature_dir, indices_dir, port = 3838)


# sync to s3 after integration, labeling, etc (whatever want pre-done for app)
# tar cvzf example_data.tar.gz --exclude='*.fastq.gz' --exclude='kallisto*/*' --exclude '*/barcodes.tsv' --exclude '*/matrix.mtx' --exclude '*/genes.tsv' --directory=/home/alex/patient_data example
# aws s3 cp example_data.tar.gz s3://dseqr/example_data.tar.gz
