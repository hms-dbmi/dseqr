# this script sets up CMAP02_pdata.rds and L1000_pdata.rds
# which are used to match Pubchem CIDS to signatures

library(Biobase)
data_dir <- file.path('inst', 'extdata')

# setup pdata for CMAP02 -----

# load CMAP from GEOdb
cmap_eset <- readRDS('/home/alex/Documents/Batcave/GEO/GEOdb/data-raw/avpick/data_dir/CMAP02/CMAP02_eset.rds')
cmap_pdata <- pData(cmap_eset[[1]])

# setup so that matches drug_es
cmap_es_path <- '/home/alex/Documents/Batcave/GEO/ccdata/data-raw/cmap_es/cmap_es_ind.rds'
cmap_es <- readRDS(cmap_es_path)

# remove columns that don't need right now
# smiles and trt_nums (from rnama.com to match GSE1_p1-p2 contrasts from predicted groups)
cmap_pdata$trt_nums <- cmap_pdata$characteristics_ch1.3 <- NULL

# pretty up pubchem cids
cmap_pdata$`Pubchem CID` <- gsub('^cid: ', '', cmap_pdata$characteristics_ch1.2)
cmap_pdata$characteristics_ch1.2 <- NULL

# remove columns that can reconstruct from title (reduces package memory)
cmap_nsamp <- gsub('^nsamples: (\\d+)$', '\\1', cmap_pdata$characteristics_ch1.4)
cmap_title_nsamp <- gsub('^.+?_(\\d+)$', '\\1', cmap_pdata$title)

if (all(cmap_nsamp == cmap_title_nsamp))  {
  cat('nsamples are correct\n')
  cmap_pdata$`Samples(n)` <- cmap_nsamp
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
all(colnames(cmap_es) %in% cmap_title_no_nsamp)

# replace title with match and put in order of cmap_es
row.names(cmap_pdata) <- cmap_pdata$title <- cmap_title_no_nsamp
cmap_pdata <- cmap_pdata[colnames(cmap_es), ]
all(cmap_pdata$title == colnames(cmap_es))
row.names(cmap_pdata) <- NULL

# save
file.copy(cmap_es_path, file.path(data_dir, 'cmap_es_ind.rds'))
saveRDS(cmap_pdata, file.path(data_dir, 'CMAP02_pdata.rds'))

# setup pdata for L1000 ------

l1000_eset <- readRDS('/home/alex/Documents/Batcave/GEO/GEOdb/data-raw/avpick/data_dir/L1000/L1000_eset.rds')
l1000_pdata <- pData(l1000_eset[[1]])

# setup so that matches drug_es
l1000_es <- readRDS('/home/alex/Documents/Batcave/GEO/l1000/level1/6-limma/l1000_es.rds')

# remove columns that don't need right now
# smiles and trt_nums (from rnama.com to match GSE1_p1-p2 contrasts from predicted groups)
l1000_pdata$trt_nums <- l1000_pdata$characteristics_ch1.3 <- NULL

# pretty up pubchem cids
l1000_pdata$`Pubchem CID` <- gsub('^cid: ', '', l1000_pdata$characteristics_ch1.2)
l1000_pdata$characteristics_ch1.2 <- NULL

# remove columns that can reconstruct from title (reduces package memory)
if (all(row.names(l1000_pdata) == l1000_pdata$title)) {
  cat('Removing row.names\n')
  row.names(l1000_pdata) <- NULL
}

l1000_title_nsamp <- gsub('^.+?_(\\d+)$', '\\1', l1000_pdata$title)
l1000_nsamp <- gsub('^nsamples: ', '', l1000_pdata$characteristics_ch1.4)
if (all(l1000_title_nsamp == l1000_nsamp)) {
  cat('nsamples are correct\n')
  l1000_pdata$`Samples(n)` <- l1000_nsamp
  l1000_pdata$characteristics_ch1.4 <- NULL
}

# make l1000_es column names match l1000_pdata titles
l1000_titles_no_nsamp <-  gsub('^(.+?)_\\d+$', '\\1', l1000_pdata$title)

# -sh, -lig, and -oe to _sh, _lig, and _oe (from rnama.com legacy)
l1000_titles_no_nsamp <- gsub('-sh_', '_sh_', l1000_titles_no_nsamp)
l1000_titles_no_nsamp <- gsub('-oe_', '_oe_', l1000_titles_no_nsamp)
l1000_titles_no_nsamp <- gsub('-lig_', '_lig_', l1000_titles_no_nsamp)
table(colnames(l1000_es) %in% l1000_titles_no_nsamp)

# pdata in same order as es
row.names(l1000_pdata) <- l1000_pdata$title <- l1000_titles_no_nsamp
l1000_pdata <- l1000_pdata[colnames(l1000_es), ]

# fix compound name mangling in titles from data.frame check.names issues
l1000_cmp <- gsub('compound: (.+?)$', '\\1', l1000_pdata$characteristics_ch1.1)
l1000_cmp <- gsub('_', '-', l1000_cmp) # makes easy to extract (no '_' within compound name)

l1000_titles_no_cmp <- gsub('.+?_([^_]+_[^_]+_[^_]+)$', '\\1', l1000_pdata$title)
l1000_titles_fixed <- paste(l1000_cmp, l1000_titles_no_cmp, sep='_')

colnames(l1000_es) <- l1000_pdata$title <- l1000_titles_fixed

# check
l1000_titles_cmp <- gsub('^([^_]+)_.+?$', '\\1', l1000_pdata$title)
all(l1000_titles_cmp == l1000_cmp)

row.names(l1000_pdata) <- l1000_pdata$characteristics_ch1.1 <- NULL

# remove string 'NA' pubchem cids
l1000_pdata$`Pubchem CID`[l1000_pdata$`Pubchem CID` == 'NA'] <- NA

# save pdata and overwrite existing l1000_es data with fixed names
# ALSO UPDATE S3 DATA in drugseqr BUCKET IF ANYTHING CHANGES
saveRDS(l1000_es, file.path(data_dir, 'l1000_es.rds'))
saveRDS(l1000_pdata, file.path(data_dir, 'L1000_pdata.rds'))
