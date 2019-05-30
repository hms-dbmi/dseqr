# this script is used to get reference expression signatures for single cell identification
# GSE1133 is the same dataset used by BioGPS for human gene expression by tissue

library(crossmeta)
library(Biobase)
library(celaref)
library(tibble)
library(data.table)

setwd('~/Documents/Batcave/zaklab/drugseqr/data-raw/single-cell/cell-id/')

gse_name <- 'GSE1133'
# get_raw(gse_name)
# ensql <- '/home/alex/Documents/Batcave/GEO/crossmeta/data-raw/entrezdt/ensql.sqlite'
eset <- load_raw(gse_name, ensql = ensql)[[1]]


# clean up pdata
pData(eset)[] <- lapply(pData(eset), as.character)
pData(eset)$description[1] <- 'Colorectal Adenocarcinoma'

# annotate to symbol level ----

# keep highest iqr when duplicated symbol
iqr_rows <- which_max_iqr(eset, annot)
eset <- eset[iqr_rows, ]

# use annot for feature names
featureNames(eset) <- fData(eset)[, annot]


# run 1 vs rest DE with limma ----

de_table <- contrast_each_group_to_the_rest_for_norm_ma_with_limma(
  exprs(eset),
  as_tibble(pData(eset)),
  'GSE1133',
  'geo_accession',
  'description'
)
