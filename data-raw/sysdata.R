setwd("~/Documents/Batcave/zaklab/drugseqr/")

tx2gene <- readRDS('data-raw/tx2gene/tx2gene.rds')
biogps <- readRDS('data-raw/single-cell/biogps/biogps.rds')
cell_info <- readRDS('data-raw/cell_info/cell_info.rds')
gslist <- readRDS('data-raw/padog/gslist.rds')
gs.names <- readRDS('data-raw/padog/gs.names.rds')
genes <- readRDS('data-raw/genes/genes.rds')

usethis::use_data(tx2gene, biogps, cell_info, gslist, gs.names, genes, internal = TRUE, overwrite = TRUE)
