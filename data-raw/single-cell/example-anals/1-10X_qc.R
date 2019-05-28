# adapted from: https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/tenx.html

library(scater)
library(dplyr)
library(DropletUtils)

# load/annotate data from alevin ----

# folder with 10X fastq files
data_dir <- file.path('data-raw/single-cell/example-data/Run2644-10X-Lung')

# this is LONG RUNNING
# drugseqr::run_alevin(data_dir)

# import alevin quants
alevin_dir <- file.path(data_dir, 'alevin_output', 'alevin')
alevin <- tximport::tximport(file.path(alevin_dir, 'quants_mat.gz'), type = 'alevin')$counts

# final alevin whitelist
whitelist <- read.delim(file.path(alevin_dir, 'whitelist.txt'), header = FALSE, as.is = TRUE)$V1

# Convert to SingleCellExperiment
sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=alevin))

# annotate sce using tx2gene (best match with cmap_es/l1000_es)
tx2gene <- readRDS("data-raw/tx2gene/tx2gene.rds")
tx2gene <- tx2gene %>%
  as_tibble() %>%
  dplyr::select(-c(tx_id, entrezid)) %>%
  group_by(gene_id) %>%
  summarise_all(unique)

rowData(sce) <- tibble(gene_id = gsub('\\.[0-9]+$', '', row.names(sce))) %>%
  left_join(tx2gene, by = 'gene_id')

# appends _ID to any non-unique values of names
rownames(sce) <- scater::uniquifyFeatureNames(ID = rowData(sce)$gene_id, names = rowData(sce)$gene_name)

# alevin does initial knee-based whitelisting then refines based on mtRNA, rRNA, etc
# see what things look like before/after final whitelisting
sce <- calculateQCMetrics(sce, feature_controls=list(Mito=which(rowData(sce)$seq_name=="MT")))
par(mfrow=c(1,4))
hist(sce$log10_total_counts, breaks=20, col="grey80",
     xlab="Log-total UMI count")
hist(sce$log10_total_features_by_counts, breaks=20, col="grey80",
     xlab="Log-total number of expressed features")
hist(sce$pct_counts_Mito, breaks=20, col="grey80",
     xlab="Percent mito reads")
hist(sce$pct_counts_Mito[colnames(sce) %in% whitelist], breaks=20, col="grey80",
     xlab="Percent mito reads (whitelist)")

# subset by alevin whitelist
sce <- sce[, colnames(sce) %in% whitelist]

# look at most highly expressed genes
# expect dominated by mitochondrial/ribosomal
ave <- calcAverage(sce)
rowData(sce)$AveCount <- ave
hist(log10(ave), col="grey80", breaks = 100)
plotHighestExprs(sce)

# normalizing for cells-specific biases ----
library(scran)
library(BiocSingular)
set.seed(1000)
clusters <- quickCluster(sce, use.ranks=FALSE, BSPARAM=IrlbaParam())
table(clusters)

# sizefactors should correlated well with total counts
# indicates capture efficiency/sequencing depth are the major biases
plot(sce$total_counts, sizeFactors(sce), log="xy")

# normalize based on size-factors
sce <- normalize(sce)


