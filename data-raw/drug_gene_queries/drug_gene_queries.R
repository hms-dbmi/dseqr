library(drugseqr)
library(data.table)

# drug by genetic pert searches

# load data
cmap_path <- system.file('extdata', 'cmap_es_ind.rds', package = 'drugseqr', mustWork = TRUE)
cmap_es <- readRDS(cmap_path)

l1000_genes_path <- system.file('extdata', 'l1000_genes_es.rds', package = 'drugseqr', mustWork = TRUE)
l1000_genes <- readRDS(l1000_genes_path)

l1000_drugs_path <- system.file('extdata', 'l1000_drugs_es.rds', package = 'drugseqr', mustWork = TRUE)
l1000_drugs <- readRDS(l1000_drugs_path)

cmap_compounds <- gsub('^([^_]+)_.+?$', '\\1', colnames(cmap_es))
l1000_compounds <- gsub('^([^_]+)_.+?$', '\\1', colnames(l1000_drugs))


# for each genetic signature query cmap_es and l1000_drugs_es
for (i in seq_len(ncol(l1000_genes))) {
  query_genes <- l1000_genes[, i]
 cat('Working on', i, 'of', ncol(l1000_genes), '\n')

  # get query results
  res_cmap <- query_drugs(query_genes, cmap_es)
  res_l1000 <- query_drugs(query_genes, l1000_drugs)

  # not feasible to save everything
  res_cmap <- data.table(Compound = cmap_compounds,
                         Correlation = res_cmap, key = 'Compound')

  res_l1000 <- data.table(Compound = l1000_compounds,
                         Correlation = res_l1000, key = 'Compound')

  res_cmap <- get_top_cors(res_cmap, is_genetic = TRUE)
  res_l1000 <- get_top_cors(res_l1000, is_genetic = TRUE)
}


