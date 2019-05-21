setwd("~/Documents/Batcave/zaklab/drugseqr/data-raw/sysdata")

tx2gene <- readRDS('tx2gene.rds')

usethis::use_data(tx2gene, internal = TRUE, overwrite = TRUE)
