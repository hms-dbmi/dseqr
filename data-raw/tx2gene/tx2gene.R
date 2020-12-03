# get tx2gene
library(GEOkallisto)
library(tidyr)
library(dplyr)
library(tibble)

tx2gene_human <- get_tx2gene()
tx2gene_mouse <- get_tx2gene(species = 'Mus musculus', release = '98')

# get entrezid --> HGNC map used by cmap_es_ind and l1000_es
load("~/Documents/Batcave/GEO/crossmeta/R/sysdata.rda")
hs$ENTREZID <- as.integer(hs$ENTREZID)

# TODO: figure out homologene stuff for mouse so that can use drugs tab
tx2gene <- tx2gene_mouse

# unnest and join on entrezid
tx2gene_unnest <- unnest(tx2gene, entrezid)
tx2gene_unnest <- left_join(tx2gene_unnest, hs, by = c('entrezid' = 'ENTREZID'))


summarise_symbol <- function(gnames, syms, species = 'human') {

  if (species == 'human') {
    # use SYMBOL_9606 if available
    syms <- unique(toupper(na.omit(syms)))
    if (length(syms) == 1) return(syms)

    # otherwise gene_name
    gnames <- unique(toupper(na.omit(gnames)))
    return(gnames)

  } else {
    unique(na.omit(gnames))
  }
}

# re-nest
tx2gene <- group_by(tx2gene_unnest, tx_id) %>%
  summarise(gene_name = summarise_symbol(gene_name, SYMBOL_9606, species = 'mouse'),
            entrezid = entrezid[1],
            gene_id = unique(gene_id),
            seq_name = unique(seq_name),
            description = unique(description))

# need tx_id, gene_name and entrezid for load_seq
# need gene_id for annotation
# need seq_name for mitochondrial genes
# need description for app
saveRDS(tx2gene, 'data-raw/tx2gene/tx2gene_mouse.rds')

# check concordance with l1000_es/cmap_es ----
data_dir <- system.file('extdata', package = 'drugseqr.data')
cmap_es_ind <- readRDS(file.path(data_dir, "inst/extdata/cmap_es_ind.rds"))
l1000_es <- readRDS(file.path(data_dir, "inst/extdata/l1000_es.rds"))

sum(!row.names(cmap_es_ind) %in% tx2gene$gene_name)
# EnsDb89: 207
# EnsDb90: 203
# EnsDb91: 201
# EnsDb92: 220
# EnsDb94: 169
# EnsDb95: 197
sum(!row.names(l1000_es) %in% tx2gene$gene_name)
# EnsDb89: 4
# EnsDb90: 4
# EnsDb91: 4
# EnsDb92: 4
# EnsDb94: 3
# EnsDb95: 4

