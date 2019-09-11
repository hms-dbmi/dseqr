pert_names <- list.files('data-raw/drug_gene_queries/data')

pert_names <- gsub('.rds$', '', pert_names)
pert_names <- gsub('^.+?_res_', '', pert_names)
pert_names <- unique(pert_names)

saveRDS(pert_names, 'data-raw/drug_gene_queries/pert_names.rds')
