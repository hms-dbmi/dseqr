# adapted from: https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/tenx.html

library(scater)

# folder with 10X fastq files
data_dir <- file.path('data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg')

# this is LONG RUNNING
drugseqr::run_alevin(data_dir)

# import alevin quants
alevin_dir <- file.path(data_dir, 'alevin_output', 'alevin')
alevin <- tximport::tximport(file.path(alevin_dir, 'quants_mat.gz'), type = 'alevin')

# Subset to only cells in whitelist
whitelist <- read.delim(file.path(alevin_dir, 'whitelist.txt'), header = FALSE, as.is = TRUE)$V1
alevin$whitelist <- colnames(alevin$counts) %in% whitelist
alevin$counts <- alevin$counts[, alevin$whitelist]

# Convert to SingleCellExperiment
sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=alevin$counts))

# annotate sce using EnsDB package




