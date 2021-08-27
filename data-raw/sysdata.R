biogps <- readRDS('data-raw/biogps/biogps.rds')
cell_info <- readRDS('data-raw/cell_info/cell_info.rds')
genes <- readRDS('data-raw/genes/genes.rds')
pert_names <- readRDS('data-raw/drug_gene_queries/pert_names.rds')
azimuth_refs <- c('human_pbmc', 'human_lung', 'human_motorcortex', 'mouse_motorcortex')
ensmap <- readRDS('data-raw/ensmap/ensmap.rds')

usethis::use_data(biogps, cell_info, genes, pert_names, azimuth_refs, ensmap, internal = TRUE, overwrite = TRUE)
