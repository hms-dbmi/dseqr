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
tx2gene <- group_by(tx2gene_unnest, tx_id) %>%
  summarise(gene_name = summarise_symbol(gene_name, SYMBOL_9606, entrezid),
            entrezid = list(unique(entrezid)),
            gene_id = unique(gene_id),
            seq_name = unique(seq_name))

# only need tx_id and gene_name for salmon quant
# need gene_id for annotation
# need seq_name to for mitochondrial genes
tx2gene %>%
  dplyr::select(tx_id, gene_name, gene_id, seq_name) %>%
  saveRDS('data-raw/tx2gene/tx2gene.rds')

# check concordance with l1000_es/cmap_es ----
cmap_es_ind <- readRDS("inst/extdata/cmap_es_ind.rds")
l1000_es <- readRDS("inst/extdata/l1000_es.rds")

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

# get mitochondrial and ribosomal genes for alevin ----
#  for --mrna and --rrna flag of alevin
# described here: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1670-y#Sec16

library('org.Hs.eg.db')
library(data.table)

# use gene-ontology to retrieve ribo genes
# sample approach as scPipe calculate_QC_metrics: https://github.com/LuyiTian/scPipe/blob/master/R/qc.R
ribo_go <- 'GO:0005840'

AnnotationDbi::columns(org.Hs.eg.db)

ribo_genes <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=ribo_go, column='ENSEMBLTRANS', keytype="GO",multiVals = "list")[[1]]
ribo_genes <- tibble(tx_id = unique(na.omit(ribo_genes)))

mito_genes <- tx2gene %>%
  dplyr::select(tx_id, seq_name) %>%
  group_by(tx_id) %>%
  summarise_all(unique) %>%
  dplyr::filter(seq_name == 'MT')

# make sure genes include subversion as used by alevin
txp2gene <- fread('inst/extdata/txp2gene.tsv', col.names = c('tx_id', 'gene_id'))
txp2gene <- as_tibble(txp2gene) %>%
  mutate(tx_id = gsub('\\.[0-9]+$', '', tx_id))

mito_genes <- mito_genes %>%
  left_join(txp2gene, by = 'tx_id') %>%
  dplyr::filter(!is.na(gene_id)) %>%
  pull(gene_id)

ribo_genes <- ribo_genes %>%
  left_join(txp2gene, by = 'tx_id') %>%
  dplyr::filter(!is.na(gene_id)) %>%
  dplyr::select(gene_id) %>%
  dplyr::distinct() %>%
  pull(gene_id)


write.table(mito_genes, 'inst/extdata/mrna.csv', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(ribo_genes, 'inst/extdata/rrna.csv', quote = FALSE, row.names = FALSE, col.names = FALSE)



