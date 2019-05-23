setwd("~/Documents/Batcave/zaklab/drugseqr/")

tx2gene <- readRDS('data-raw/tx2gene/tx2gene.rds')

usethis::use_data(tx2gene, internal = TRUE, overwrite = TRUE)
