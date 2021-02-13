# saves all CMAP02 and L1000 adjusted pvalues as individual files for quick loading for dseqr::heatplot
library(dseqr)
# load data
l1000_pvals_path <- system.file('extdata', 'l1000_pval.adj.rds', package = 'dseqr.data', mustWork = TRUE)
cmap_pvals_path <- system.file('extdata', 'cmap_pval.adj_ind.rds', package = 'dseqr.data', mustWork = TRUE)
l1000_pvals <- readRDS(l1000_pvals_path)
cmap_pvals <- readRDS(cmap_pvals_path)

# split l1000 by compounds/genetic perts & ligands
is.genetic <- grepl('-oe_|-sh_|-lig_', colnames(l1000_pvals))

l1000_genes_pvals <- l1000_pvals[, is.genetic]
l1000_drugs_pvals <- l1000_pvals[, !is.genetic]

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

save_signatures(l1000_genes_pvals, 'data-raw/drug_es/pvals/l1000_genes')
save_signatures(l1000_drugs_pvals, 'data-raw/drug_es/pvals/l1000_drugs')
save_signatures(cmap_pvals, 'data-raw/drug_es/pvals/cmap')

# sync to s3
# aws s3 sync ~/Documents/Batcave/zaklab/dseqr/data-raw/drug_es/pvals s3://dseqr/drug_pvals_dir
