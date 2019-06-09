setwd("~/Documents/Batcave/zaklab/drugseqr/")

tx2gene <- readRDS('data-raw/tx2gene/tx2gene.rds')
biogps <- readRDS('data-raw/single-cell/biogps/biogps.rds')

usethis::use_data(tx2gene, biogps, internal = TRUE, overwrite = TRUE)
