biogps <- readRDS('data-raw/single-cell/biogps/biogps.rds')
cell_info <- readRDS('data-raw/cell_info/cell_info.rds')
genes <- readRDS('data-raw/genes/genes.rds')
pert_names <- readRDS('data-raw/drug_gene_queries/pert_names.rds')


usethis::use_data(biogps, cell_info, genes, pert_names, internal = TRUE, overwrite = TRUE)
