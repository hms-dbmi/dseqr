# This script explores using celaref to identify cell clusters in single-cell RNA-seq datasets

library(crossmeta)
library(Biobase)
library(celaref)
library(tibble)
library(data.table)
library(drugseqr)

setwd('~/Documents/Batcave/zaklab/drugseqr/')

# get reference expression signatures for single cell identification ----

# GSE1133 is the same dataset used by BioGPS for human gene expression by tissue
gse_name <- 'GSE1133'
data_dir <- 'data-raw/single-cell/example-anals'
# get_raw(gse_name)
# ensql <- '/home/alex/Documents/Batcave/GEO/crossmeta/data-raw/entrezdt/ensql.sqlite'
eset <- load_raw(gse_name, data_dir, ensql = ensql)[[1]]


# clean up pdata
pData(eset)[] <- lapply(pData(eset), as.character)
pData(eset)$description[1] <- 'Colorectal Adenocarcinoma'

# annotate to symbol level

# keep highest iqr when duplicated symbol
iqr_rows <- which_max_iqr(eset, 'SYMBOL')
eset <- eset[iqr_rows, ]

# use annot for feature names
featureNames(eset) <- fData(eset)[, 'SYMBOL']

# run 1 vs rest DE with limma ----

de_table.ref <- contrast_each_group_to_the_rest_for_norm_ma_with_limma(
  exprs(eset),
  as_tibble(pData(eset)),
  'GSE1133',
  'geo_accession',
  'description'
)

# load 10X lung data ----

# folder with 10X fastq files
data_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'

# import alevin quants
sce <- load_scseq(data_dir)

# run qc
sce <- qc_scseq(sce)

# subset by alevin whitelist norm/stabilize using good cell only
sce <- sce[, sce$whitelist]
sce <- norm_scseq(sce)
sce <- stabilize_scseq(sce)

# clusters -----

sce <- add_scseq_clusters(sce)

se <- SummarizedExperiment(SimpleList(counts = exprs(sce)),
                           colData = data.frame(group = sce$cluster),
                           rowData = data.frame(ID = row.names(sce), stringsAsFactors = FALSE))

se <- trim_small_groups_and_low_expression_genes(se)
de_table.lung <- contrast_each_group_to_the_rest(se, dataset_name="Lung", num_cores=7)

label_table <- make_ref_similarity_names(de_table.test=de_table.lung, de_table.ref=de_table.ref)
label_table$shortlab

# [1] "2:(bone marrow-CD33Myeloid|peripheral blood-CD14Monocytes|WHOLEBLOOD)"
# [2] "3:(Prostate)"
# [3] "4:No similarity
