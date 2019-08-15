data_dir <- file.path('inst', 'extdata')

l1000_es <- readRDS(file.path(data_dir, 'l1000_genes_es.rds'))
cmap_es <- readRDS(file.path(data_dir, 'cmap_es_ind.rds'))

genes <- unique(c(row.names(cmap_es), row.names(l1000_es)))
saveRDS(genes, 'data-raw/genes/genes.rds')
