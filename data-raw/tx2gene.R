# get tx2gene
library(drugseqr)
library(tidyr)
library(dplyr)

tx2gene <- get_tx2gene()

# get entrezid --> HGNC map used by cmap_es_ind and l1000_es
load("~/Documents/Batcave/GEO/crossmeta/R/sysdata.rda")
hs$ENTREZID <- as.integer(hs$ENTREZID)

# unnest and join on entrezid
tx2gene_unnest <- unnest(tx2gene, entrezid)
tx2gene_unnest <- left_join(tx2gene_unnest, hs, by = c('entrezid' = 'ENTREZID'))

summarise_symbol <- function(gnames, syms) {
  # use SYMBOL_9606 if available
  syms <- unique(na.omit(syms))
  if (length(syms)) return(list(syms))

  # otherwise use gene_name
  gnames <- unique(na.omit(gnames))
  return(list(gnames))
}


# re-nest
tx2gene_nest <- group_by(tx2gene_unnest, tx_id) %>%
  summarise(entrezid = list(unique(entrezid)),
            gene_name = unique(gene_name),
            SYMBOL_9606 = summarise_symbol(gene_name, SYMBOL_9606))
