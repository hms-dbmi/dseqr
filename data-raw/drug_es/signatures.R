# saves all CMAP02 and L1000 signatures as individual files for quick loading in Pathways tab
library(dseqr)
# load data
cmap_path <- system.file('extdata', 'cmap_es_ind.rds', package = 'dseqr.data', mustWork = TRUE)
cmap_es <- readRDS(cmap_path)

l1000_genes_path <- system.file('extdata', 'l1000_genes_es.rds', package = 'dseqr.data', mustWork = TRUE)
l1000_genes <- readRDS(l1000_genes_path)

l1000_drugs_path <- system.file('extdata', 'l1000_drugs_es.rds', package = 'dseqr.data', mustWork = TRUE)
l1000_drugs <- readRDS(l1000_drugs_path)

# save each signature seperately
save_signatures <- function(es, data_dir) {
  sig_names <- colnames(es)
  dir.create(data_dir, recursive = TRUE)


  for (sig_name in sig_names) {
    sig <- es[,sig_name]
    fname <- paste0(sig_name, '.rds')
    fpath <- file.path(data_dir, fname)
    if (file.exists(fpath)) next()

    tryCatch(saveRDS(sig, fpath),
             error = function(e) {
               new_fname <- fs::path_sanitize(fname)
               new_fpath <- file.path(data_dir, new_fname)

               cat(fname, 'invalid, changing to:', new_fname, '\n')
               saveRDS(sig, new_fpath)
             })
  }
}

save_signatures(cmap_es, 'data-raw/drug_es/signatures/cmap')
save_signatures(l1000_genes, 'data-raw/drug_es/signatures/l1000_genes')
save_signatures(l1000_drugs, 'data-raw/drug_es/signatures/l1000_drugs')

# sync to s3
# aws s3 sync ~/Documents/Batcave/zaklab/dseqr/data-raw/drug_es/signatures s3://dseqr/drug_es_dir
