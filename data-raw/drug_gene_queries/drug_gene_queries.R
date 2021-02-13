library(dseqr)
library(data.table)
library(rlang)

# drug by genetic pert searches

# load data
# cmap_path <- system.file('extdata', 'cmap_es_ind.rds', package = 'dseqr.data', mustWork = TRUE)
# cmap_es <- readRDS(cmap_path)

l1000_genes_path <- system.file('extdata', 'l1000_genes_es.rds', package = 'dseqr.data', mustWork = TRUE)
l1000_genes <- readRDS(l1000_genes_path)

l1000_drugs_path <- system.file('extdata', 'l1000_drugs_es.rds', package = 'dseqr.data', mustWork = TRUE)
l1000_drugs <- readRDS(l1000_drugs_path)

# cmap_compounds <- gsub('^([^_]+)_.+?$', '\\1', colnames(cmap_es))
l1000_compounds <- gsub('^([^_]+)_.+?$', '\\1', colnames(l1000_drugs))
l1000_genetic <- gsub('^([^_]+)_.+?$', '\\1', colnames(l1000_genes))

# check that l1000 and cmap compounds are unique
# length(intersect(colnames(cmap_es), colnames(l1000_drugs)))
# [1] 0

save_dir <- 'data-raw/drug_gene_queries/data'
dir.create(save_dir)

#' Run queries for drug by gene signatures
#'
#' @param drug_es matrix that is queried against.
#' @param query_es matrix of signatures to query with.
#' @param prefix Appended to start of saved query results.
#'   Either 'cmap_res_', 'l1000_genes_res_', or 'l1000_drugs_res_'. Use based on \code{drug_es}.
#' @param compounds The compound names for \code{drug_es}.
run_pert_queries <- function(drug_es, query_es, prefix, compounds) {

  is_genetic <- prefix == 'l1000_genes_res_'

  # progress
  iterations <- ncol(query_es)
  pb <- txtProgressBar(max = iterations, style = 3)

  for (i in seq_len(iterations)) {
    # remove illegal path characters
    query_name <- colnames(query_es)[i]
    query_name <- fs::path_sanitize(query_name)

    query_path <- file.path(save_dir, paste0(prefix, query_name, '.rds'))

    # get query results
    query_genes <- query_es[, i]
    res <- dseqr::query_drugs(query_genes, drug_es)

    # not feasible to save everything
    # save vector of top results
    res_table <- data.table::data.table(Compound = compounds,
                                        Correlation = res, key = 'Compound')

    # always run with is_genetic TRUE so that can show most similar absolute results for pert queries
    top_cors <- dseqr::get_top_cors(res_table, is_genetic = TRUE)
    res <- res[res_table$Compound %in% top_cors$Compound]

    saveRDS(res, query_path)
    setTxtProgressBar(pb, i)
  }
}

# run each query_es against each drug_es (cmap_es, l1000_drugs, l1000_genes)
# run_pert_queries(cmap_es, cmap_es, 'cmap_res_', cmap_compounds)
# run_pert_queries(cmap_es, l1000_genes, 'cmap_res_', cmap_compounds)
# run_pert_queries(cmap_es, l1000_drugs, 'cmap_res_', cmap_compounds)

# run_pert_queries(l1000_genes, cmap_es, 'l1000_genes_res_', l1000_genetic)
# run_pert_queries(l1000_genes, l1000_genes, 'l1000_genes_res_', l1000_genetic)
# run_pert_queries(l1000_genes, l1000_drugs, 'l1000_genes_res_', l1000_genetic)

# run_pert_queries(l1000_drugs, cmap_es, 'l1000_drugs_res_', l1000_compounds)
# run_pert_queries(l1000_drugs, l1000_genes, 'l1000_drugs_res_', l1000_compounds)
# run_pert_queries(l1000_drugs, l1000_drugs, 'l1000_drugs_res_', l1000_compounds)


# sync to s3
# aws s3 sync ~/Documents/Batcave/zaklab/dseqr/data-raw/drug_gene_queries/data s3://dseqr/pert_query_dir
