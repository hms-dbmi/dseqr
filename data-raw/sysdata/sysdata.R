setwd("~/Documents/Batcave/zaklab/drugseqr/")

tx2gene <- readRDS('data-raw/tx2gene/tx2gene.rds')
biogps <- readRDS('data-raw/single-cell/biogps/biogps.rds')
cell_info <- readRDS('data-raw/cell_info/cell_info.rds')

usethis::use_data(tx2gene, biogps, cell_info, internal = TRUE, overwrite = TRUE)
