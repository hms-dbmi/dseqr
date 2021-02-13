data_dir <- system.file('extdata', package = 'dseqr.data')

l1000_es <- readRDS(file.path(data_dir, 'l1000_genes_es.rds'))
cmap_es <- readRDS(file.path(data_dir, 'cmap_es_ind.rds'))

cmap_genes <- row.names(cmap_es)
l1000_genes <- row.names(l1000_es)
common <- intersect(cmap_genes, l1000_genes)

# all l1000 genes are also in CMAP
genes <- list(cmap_only = setdiff(cmap_genes, common), common = common)

saveRDS(genes, 'data-raw/genes/genes.rds')
