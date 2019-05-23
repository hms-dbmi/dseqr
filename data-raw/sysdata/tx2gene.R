# get tx2gene
library(drugseqr)
library(tidyr)
library(dplyr)
library(tibble)

tx2gene <- get_tx2gene()

# get entrezid --> HGNC map used by cmap_es_ind and l1000_es
load("~/Documents/Batcave/GEO/crossmeta/R/sysdata.rda")
hs$ENTREZID <- as.integer(hs$ENTREZID)

# unnest and join on entrezid
tx2gene_unnest <- unnest(tx2gene, entrezid)
tx2gene_unnest <- left_join(tx2gene_unnest, hs, by = c('entrezid' = 'ENTREZID'))



summarise_symbol <- function(gnames, syms, entrezid) {

  # use SYMBOL_9606 if available
  syms <- unique(toupper(na.omit(syms)))
  if (length(syms) == 1) return(syms)

  # otherwise gene_name
  gnames <- unique(toupper(na.omit(gnames)))
  return(gnames)
}

# re-nest
tx2gene_nest <- group_by(tx2gene_unnest, tx_id) %>%
  summarise(gene_name = summarise_symbol(gene_name, SYMBOL_9606, entrezid),
            entrezid = list(unique(entrezid)))


# check concordance with l1000_es/cmap_es
cmap_es_ind <- readRDS("inst/extdata/cmap_es_ind.rds")
l1000_es <- readRDS("inst/extdata/l1000_es.rds")

sum(!row.names(cmap_es_ind) %in% tx2gene_nest$gene_name)
# EnsDb89: 207
# EnsDb90: 203
# EnsDb91: 201
# EnsDb92: 220
# EnsDb94: 169
sum(!row.names(l1000_es) %in% tx2gene_nest$gene_name)
# EnsDb89: 4
# EnsDb90: 4
# EnsDb91: 4
# EnsDb92: 4
# EnsDb94: 3


saveRDS(tx2gene_nest, 'data-raw/sysdata/tx2gene.rds')
