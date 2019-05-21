tx2gene <- readRDS("~/Documents/Batcave/zaklab/drugseqr/data-raw/sysdata/tx2gene.rds")

tgmap <- tx2gene[, 1:2]
write.table(tgmap, 'inst/extdata/tgmap.tsv', row.names = FALSE, col.names = FALSE, quote = FALSE)


