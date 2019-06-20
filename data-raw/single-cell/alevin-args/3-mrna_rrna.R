# get mitochondrial and ribosomal genes for alevin ----
#  for --mrna and --rrna flag of alevin
# described here: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1670-y#Sec16

library(data.table)

tx2gene <- readRDS('data-raw/tx2gene/tx2gene.rds')

# ribo genes from: https://www.genenames.org/data/genegroup/#!/group/1054
ribo_genes <- fread('data-raw/single-cell/alevin-args/ribo.txt')
ribo_genes <- c(ribo_genes$`Approved symbol`, ribo_genes$`Previous symbols`)
ribo_genes <- tx2gene$gene_name[tx2gene$gene_name %in% toupper(ribo_genes)]
ribo_genes <- unique(ribo_genes)

mito_genes <- tx2gene %>%
  dplyr::select(seq_name, gene_name) %>%
  dplyr::filter(seq_name == 'MT') %>%
  dplyr::select(gene_name) %>%
  distinct() %>%
  pull(gene_name)

write.table(mito_genes, 'inst/extdata/mrna.csv', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(ribo_genes, 'inst/extdata/rrna.csv', quote = FALSE, row.names = FALSE, col.names = FALSE)

