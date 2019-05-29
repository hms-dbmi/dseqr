# this script generate a transcript --> HGNC map for --tgMap in alevin
# avoids subsequest remapping to HGNC symbol and resolution of many-to-one relationships
# using tx2gene.rds provides closest match with l1000_es and cmap_es annotations

library(data.table)
library(dplyr)

txp2gene <- fread('inst/extdata/txp2gene.tsv', col.names = c('tx_id', 'gene_id'))
tx2gene <- readRDS('data-raw/tx2gene/tx2gene.rds')

txp2hgnc <- txp2gene %>%
  mutate(tx_id_sub = tx_id,
         tx_id = gsub('\\.\\d+$', '', tx_id)) %>%
  left_join(tx2gene, by='tx_id') %>%
  select(tx_id_sub, gene_name) %>%
  rename(tx_id = tx_id_sub)

write.table(txp2hgnc, 'inst/extdata/txp2hgnc.tsv', quote = FALSE, row.names = FALSE, col.names = FALSE)


