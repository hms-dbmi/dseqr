# studies with both RNA-Seq and Microarray assays of CMAP02 treated drugs

# RNA-Seq    | Microarray    | cell line/tissue  | drugs
# -----------|---------------|-------------------|-----
# GSE115609  |  GSE115564    | MCF7              | fulvestrant
# GSE55347   |  GSE47875     | liver             | estradiol, bezafibrate, clofibrate, clotrimazole, econazole, gemfibrozil, ifosfamide, leflunomide, lovastatin, miconazole, pirinixic acid, rosiglitazone, simvastatin
# GSE41586   |  GSE41364     | HT29              | azacitidine
# GSE43526   |  GSE43899     | cortical neurons  | topotecan (structurally similar to irinotecan and camptothecin)
# GSE122315  |  GSE122184    |liver              | paracetamol (acetaminophen), diclofenac [!LIVER TOXICANTS!]


library(crossmeta)
library(drugseqr)
library(Biobase)


# microarray GSE47875 ----
data_dir <- 'data-raw/benchmarks/rnaseq_microarray_cmap'
gse_name <- 'GSE47875'

crossmeta::get_raw(gse_name, data_dir)
eset <- crossmeta::load_raw(gse_name, data_dir)


# RNA-Seq GSE55347 ----
gse_name <- 'GSE55347'
data_dir <- file.path('/mnt/shared', gse_name)

srp_meta <- get_srp_meta(gse_name, '/mnt/shared')
species <- unique(srp_meta$organism)

# download the data
# get_fastqs(gse_name, srp_meta, '/mnt/shared')

# setup the pairs
indices_dir <- '~/Documents/Batcave/zaklab/drugseqr.data/inst/indices'
fastqs <- list.files(data_dir, '.fastq.gz')
pdata <- tibble::tibble('File Name' = fastqs)
pdata <- tibble::add_column(pdata, Pair = NA, Replicate = NA, .before = 1)
pdata$Pair <- rep(seq_len(nrow(pdata)/2), each = 2)

# run quantification
run_kallisto_bulk(indices_dir, data_dir, pdata = pdata, species = species)
run_salmon_bulk(indices_dir, data_dir, pdata = pdata, species = species)

# load quants (run for both kallisto and salmon)
type <- 'kallisto'
eset <- drugseqr::load_seq(data_dir, type = type, species = species)

# merge with srp_meta
srr_names <- stringr::str_extract(colnames(eset), 'SRR\\d+')
colnames(eset) <- srr_names
pData(eset) <- cbind(pData(eset), srp_meta[colnames(eset), ])
colnames(eset) <- pData(eset)$gsm_name

# merge with GEO pdata
esets <- crossmeta:::getGEO(gse_name)
pdatas <- lapply(esets, function(eset) pData(eset))
pdata <- rbindlist(pdatas)
pdata <- as.data.frame(pdata)
row.names(pdata) <- pdata$geo_accession
pdata <- pdata[pData(eset)$gsm_name, ]
pdata <- cbind(pData(eset), pdata)
pdata <- pdata[, !duplicated(colnames(pdata))]
pData(eset) <- data.frame(pdata)

# setup for differential expression analyses
chemicals <- as.character(pdata$characteristics_ch1.3)
chemicals <- gsub('^chemical: ', '', chemicals)
unique_chemicals <- setdiff(chemicals, 'Vehicle')

vehicles <- as.character(pdata$characteristics_ch1.5)

# just run for drugs also assayed in CMAP
# estradiol, bezafibrate, clofibrate, clotrimazole, econazole, gemfibrozil, ifosfamide, leflunomide, lovastatin, miconazole, pirinixic acid, rosiglitazone, simvastatin
cmap_drugs <- c('BETA-ESTRADIOL', 'BEZAFIBRATE', 'CLOFIBRIC ACID', 'ECONAZOLE', 'GEMFIBROZIL', 'IFOSFAMIDE',
                'LEFLUNOMIDE', 'LOVASTATIN', 'MICONAZOLE', 'PIRINIXIC ACID', 'ROSIGLITAZONE', 'SIMVASTATIN')

pkg_version <- get_pkg_version(type)
prefix <- paste(type, pkg_version, sep = '_')

for (drug in cmap_drugs) {
  anal_pdata <- pdata
  anal_pdata$group <- NA

  # use as control samples treated with the vehicle for the drug
  is.drug <- chemicals == drug
  is.ctrl <- vehicles == unique(vehicles[is.drug]) & !is.drug

  anal_pdata$group[is.drug] <- 'test'
  anal_pdata$group[is.ctrl] <- 'ctrl'
  anal_pdata <- anal_pdata[!is.na(anal_pdata$group), ]

  # run differential expression analysis
  anal_name <- paste(prefix, drug, sep = '_')
  anal <- diff_expr(eset, data_dir, anal_name = anal_name, prev_anal = list(pdata = anal_pdata))
}
