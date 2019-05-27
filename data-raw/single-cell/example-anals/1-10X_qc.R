# adapted from: https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/tenx.html

library(scater)
library(dplyr)
library(DropletUtils)

# folder with 10X fastq files
data_dir <- file.path('data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg')

# this is LONG RUNNING
# drugseqr::run_alevin(data_dir)

# import alevin quants
alevin_dir <- file.path(data_dir, 'alevin_output', 'alevin')
alevin <- tximport::tximport(file.path(alevin_dir, 'quants_mat.gz'), type = 'alevin')$counts

# Subset to only cells in whitelist
whitelist <- read.delim(file.path(alevin_dir, 'whitelist.txt'), header = FALSE, as.is = TRUE)$V1
alevin <- alevin[, colnames(alevin) %in% whitelist]

# Convert to SingleCellExperiment
sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=alevin))

# annotate sce using tx2gene (best match with cmap_es/l1000_es)
tx2gene <- readRDS("data-raw/tx2gene/tx2gene.rds")
tx2gene <- tx2gene %>%
  as_tibble() %>%
  dplyr::select(-c(tx_id, entrezid)) %>%
  group_by(gene_id) %>%
  summarise_all(unique) %>%
  ungroup()

rowData(sce) <- tibble(gene_id = gsub('\\.[0-9]+$', '', row.names(sce))) %>%
  left_join(tx2gene, by = 'gene_id')

# appends _ID to any non-unique values of names
rownames(sce) <- scater::uniquifyFeatureNames(ID = rowData(sce)$gene_id, names = rowData(sce)$gene_name)

# mitochondrial gene
summary(rowData(sce)$seq_name=="MT")
# Mode   FALSE    TRUE
# logical   20271      13


