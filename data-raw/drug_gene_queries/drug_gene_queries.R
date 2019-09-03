library(drugseqr)
library(data.table)
library(foreach)
library(doSNOW)

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

#setup parallel backend to use many processors
cl <- makeCluster(4) #not to overload your computer
registerDoSNOW(cl)

# progress
iterations <- ncol(l1000_genes)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# for each genetic signature query drug signatures (cmap_es)
res_cmap <- foreach(i = seq_len(iterations), .options.snow = opts) %dopar% {
  require(drugseqr)
  require(rlang)
  query_genes <- l1000_genes[, i]
  query_name <- colnames(l1000_genes)[i]

  # get query results
  res <- drugseqr::query_drugs(query_genes, cmap_es)

  # not feasible to save everything
  # store vector of top results
  res_table <- data.table::data.table(Compound = cmap_compounds,
                    Correlation = res, key = 'Compound')

  top_cors <- drugseqr::get_top_cors(res_table, is_genetic = TRUE)
  res[res_table$Compound %in% top_cors$Compound]
}
close(pb)
stopCluster(cl)
names(res_cmap) <- colnames(l1000_genes)

# for each genetic signature query drug signatures (l1000_drugs_es)
res_l1000 <- foreach(i = seq_len(ncol(l1000_genes))) %dopar% {
  require(drugseqr)
  require(rlang)
  query_genes <- l1000_genes[, i]
  query_name <- colnames(l1000_genes)[i]

  # get query results
  res <- drugseqr::query_drugs(query_genes, l1000_drugs)

  # not feasible to save everything
  # store vector of top results
  res_table <- data.table::data.table(Compound = l1000_compounds,
                                      Correlation = res, key = 'Compound')

  top_cors <- drugseqr::get_top_cors(res_table, is_genetic = TRUE)
  res[res_table$Compound %in% top_cors$Compound]
}
close(pb)
stopCluster(cl)
names(res_cmap) <- colnames(l1000_genes)


