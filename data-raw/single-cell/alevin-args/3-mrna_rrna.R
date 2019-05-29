# get mitochondrial and ribosomal genes for alevin ----
#  for --mrna and --rrna flag of alevin
# described here: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1670-y#Sec16

library('org.Hs.eg.db')
library(data.table)

tx2gene <- readRDS('data-raw/tx2gene/tx2gene.rds')

# use gene-ontology to retrieve ribo genes
# sample approach as scPipe calculate_QC_metrics: https://github.com/LuyiTian/scPipe/blob/master/R/qc.R
ribo_go <- 'GO:0005840'

AnnotationDbi::columns(org.Hs.eg.db)

ribo_genes <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=ribo_go, column='ENSEMBLTRANS', keytype="GO",multiVals = "list")[[1]]

ribo_genes <- tibble(tx_id = unique(na.omit(ribo_genes))) %>%
  left_join(tx2gene, by='tx_id') %>%
  dplyr::select(gene_name) %>%
  distinct() %>%
  pull(gene_name)

mito_genes <- tx2gene %>%
  dplyr::select(seq_name, gene_name) %>%
  dplyr::filter(seq_name == 'MT') %>%
  dplyr::select(gene_name) %>%
  distinct() %>%
  pull(gene_name)

write.table(mito_genes, 'inst/extdata/mrna.csv', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(ribo_genes, 'inst/extdata/rrna.csv', quote = FALSE, row.names = FALSE, col.names = FALSE)

