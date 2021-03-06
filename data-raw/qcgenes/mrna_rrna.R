# get mitochondrial and ribosomal genes for QC ----
library(data.table)

tx2gene <- dseqr.data::load_tx2gene()

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

# mouse ribo genes from: http://ribosome.med.miyazaki-u.ac.jp/rpg.cgi?mode=search&org=Mus+musculus&gene=
tx2gene <- dseqr.data::load_tx2gene('Mus musculus')

ribo_genes <- read.csv('data-raw/single-cell/alevin-args/ribo_mouse.csv')
ribo_genes <- ribo_genes$Gene.Name
ribo_genes <- tx2gene$gene_name[tx2gene$gene_name %in% ribo_genes]
ribo_genes <- unique(ribo_genes)

mito_genes <- tx2gene %>%
    dplyr::select(seq_name, gene_name) %>%
    dplyr::filter(seq_name == 'MT') %>%
    dplyr::select(gene_name) %>%
    distinct() %>%
    pull(gene_name)


write.table(mito_genes, 'inst/extdata/mrna_mouse.csv', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(ribo_genes, 'inst/extdata/rrna_mouse.csv', quote = FALSE, row.names = FALSE, col.names = FALSE)


