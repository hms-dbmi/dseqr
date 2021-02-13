# from supplemental PMC6542613
library(dplyr)
library(org.Mm.eg.db)

psgs <- readRDS('data-raw/macspectrum/psgs.rds')
amdsgs <- readRDS('data-raw/macspectrum/amdsgs.rds')

# use maps that using for quantification
tx2gene_mouse <- dseqr.data::load_tx2gene('Mus musculus')

psg_map <- tx2gene_mouse[tx2gene_mouse$gene_id %in% psgs, c('gene_id', 'entrezid')]
psg_map <- unique(psg_map)
psg_map$ENTREZID <- as.character(psg_map$entrezid)

amdsg_map <- tx2gene_mouse[tx2gene_mouse$gene_id %in% amdsgs, c('gene_id', 'entrezid')]
amdsg_map <- unique(amdsg_map)
amdsg_map$ENTREZID <- as.character(amdsg_map$entrezid)

# map from mouse entrezids to human
homologene <- readRDS(system.file('extdata', 'homologene.rds', package = 'dseqr.data'))

tx2gene <- dseqr.data::load_tx2gene()
tx2gene$ENTREZID <- as.character(tx2gene$entrezid)
tx2gene <- tx2gene[, c('ENTREZID', 'gene_name')]
tx2gene <- unique(tx2gene)
tx2gene <- tx2gene[!is.na(tx2gene$ENTREZID), ]

psg_hgnc <- psg_map %>%
    left_join(homologene, by = 'ENTREZID') %>%
    left_join(tx2gene, by = c('ENTREZID_HS' = 'ENTREZID')) %>%
    filter(!is.na(gene_name)) %>%
    distinct(gene_name) %>%
    pull(gene_name)

amdsg_hgnc <- amdsg_map %>%
    left_join(homologene, by = 'ENTREZID') %>%
    left_join(tx2gene, by = c('ENTREZID_HS' = 'ENTREZID')) %>%
    filter(!is.na(gene_name)) %>%
    distinct(gene_name) %>%
    pull(gene_name)


writeLines(psg_hgnc, 'data-raw/macspectrum/psg_hgnc.txt')
writeLines(amdsg_hgnc, 'data-raw/macspectrum/amdsg_hgnc.txt')

writeLines(c(psg_hgnc, amdsg_hgnc), 'data-raw/macspectrum/psg_amdsg_hgnc.txt')

