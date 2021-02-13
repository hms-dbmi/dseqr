setwd("~/Documents/Batcave/zaklab/dseqr/data-raw/drug_annot/pug_view")

check_sider <- function(cid) {
  url <- paste0('http://sideeffects.embl.de/drugs/', cid, '/')
  is_url <- RCurl::url.exists(url)
  return(is_url)
}

cids <- read.table('cids.csv', stringsAsFactors = FALSE)$V1

sider <- rep(FALSE, length(cids))
names(sider) <- cids

for (i in 1:length(cids)) {
  cid <- cids[i]
  sider[cid] <- check_sider(cid)

  if (i %% 100 == 0 | i == length(cids)) {
    cat('Working on', i, 'of', length(cids), '\n')
    saveRDS(sider, 'sider.rds')
  }
}

sum(sider)
# 465
