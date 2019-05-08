library(drugseqr)
setwd("~/Documents/Batcave/zaklab/drugseqr/data-raw/drug_annot/pug_view")

cids <- read.table('cids.csv', stringsAsFactors = FALSE)$V1

pug_annot <- tibble::tibble(pubchem_cid = cids, drugbank = NA_character_, gras = FALSE)

for (i in 1:length(cids)) {
  cat('Working on', i, 'of', length(cids), '\n')
  cid <- cids[i]

  # load pug view
  pug_file <- file.path('views', paste0(cid, '.json'))
  if (!file.exists(pug_file)) next()
  pug_view <- rjson::fromJSON(file=pug_file)

  pug_annot[i, 'drugbank'] <- get_drugbank(pug_view)
  pug_annot[i, 'gras'] <- check_gras(pug_view)
}

sum(!is.na(pug_annot$drugbank))
# 1230

sum(pug_annot$gras)
# 7

saveRDS(pug_annot, 'pug_annot.rds')
