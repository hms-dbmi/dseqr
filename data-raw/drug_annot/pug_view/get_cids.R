# get unique compound ids for 2-get_views.sh
setwd("~/Documents/Batcave/zaklab/drugseqr/data-raw/drug_annot/pug_view")

cmap_pdata <- readRDS(system.file('extdata', 'CMAP02_pdata.rds', package = 'drugseqr'))
l1000_pdata <- readRDS(system.file('extdata', 'L1000_pdata.rds', package = 'drugseqr'))

cids <- c(cmap_pdata$`Pubchem CID`, l1000_pdata$`Pubchem CID`)
names(cids) <- c(cmap_pdata$title, l1000_pdata$title)
names(cids) <- gsub('^([^_]+)_.+?$', '\\1', names(cids))
cids <- cids[!duplicated(cids) & !is.na(cids) & cids != '-666']

write.table(cids, 'cids.csv', row.names = FALSE, col.names = FALSE, quote = FALSE)
saveRDS(cids, 'cids.rds')
