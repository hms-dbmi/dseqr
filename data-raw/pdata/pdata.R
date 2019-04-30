# this script sets up CMAP02_pdata.rds and L1000_pdata.rds
# which are used to match Pubchem CIDS to signatures

library(Biobase)
data_dir <- file.path('inst', 'extdata')

# setup pdata for CMAP02 -----

# load CMAP from GEOdb
cmap_eset <- readRDS('/home/alex/Documents/Batcave/GEO/GEOdb/data-raw/avpick/data_dir/CMAP02/CMAP02_eset.rds')
cmap_pdata <- pData(cmap_eset[[1]])

# setup so that matches drug_es
cmap_es <- readRDS(file.path(data_dir, 'cmap_es_ind.rds'))

# remove columns that don't need right now
# smiles and trt_nums (from rnama.com to match GSE1_p1-p2 contrasts from predicted groups)
cmap_pdata$trt_nums <- cmap_pdata$characteristics_ch1.3 <- NULL
row.names(cmap_pdata) <- NULL

# remove columns that can reconstruct from title (reduces package memory)
cmap_nsamp <- gsub('^nsamples: (\\d+)$', '\\1', cmap_pdata$characteristics_ch1.4)
cmap_title_nsamp <- gsub('^.+?_(\\d+)$', '\\1', cmap_pdata$title)

if (all(cmap_nsamp == cmap_title_nsamp)) {
  cat('Removing nsamples\n')
  cmap_pdata$characteristics_ch1.4 <- NULL
}

cmap_cmp <- gsub('^compound: (.+?)$', '\\1', cmap_pdata$characteristics_ch1.1)
cmap_title_cmp <- gsub('^([^_]+)_.+?$', '\\1', cmap_pdata$title)

if (all(cmap_cmp == cmap_title_cmp)) {
  cat('Removing compound\n')
  cmap_pdata$characteristics_ch1.1 <- NULL
}

# check that matches cmap_es
cmap_title_no_nsamp <- gsub('^(.+?)_\\d+$', '\\1', cmap_pdata$title)
table(colnames(cmap_es) %in% cmap_title_no_nsamp)

# save
saveRDS(cmap_pdata, file.path(data_dir, 'CMAP02_pdata.rds'))



