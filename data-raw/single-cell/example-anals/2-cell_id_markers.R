# This script explores using marker genes and expert domain knowledge to id cell clusters for single-cell RNA-seq
library(drugseqr)

# folder with 10X fastq files
data_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'

# this is LONG RUNNING
# drugseqr::run_alevin(data_dir)

# import alevin quants
sce <- load_scseq(data_dir)

# run qc
sce <- qc_scseq(sce)

# alevin does initial knee-based whitelisting then refines based on mtRNA, rRNA, etc
# compare metrics for whitelisted and non-whitelisted cells
hist_scseq_whitelist(sce)

# normalize and stabilize using all cells
sce.full <- sce
sce.full <- norm_scseq(sce.full)
sce.full <- stabilize_scseq(sce.full)

# show tSNE of all cells, colored by whitelist metrics
tsne_scseq_whitelist(sce.full)

# subset by alevin whitelist norm/stabilize using good cell only
sce <- sce[, sce$whitelist]
sce <- norm_scseq(sce)
sce <- stabilize_scseq(sce)

# clusters and marker gene detection -----

# Note: tSNE plot differs from above because only using whitelisted cells
# it is possible to maintain the above cell coordinates by subsetting sce.full
# however clusters should still be found on sce as norm/stabilization on whitelisted only

# get clusters and run tSNE
sce <- add_scseq_clusters(sce)

set.seed(1000)
sce <- scater::runTSNE(sce, use_dimred="PCA")
scater::plotTSNE(sce, colour_by = "cluster")

# only upregulated as more useful for positive id of cell type
markers <- scran::findMarkers(sce, clusters=sce$cluster, direction="up")

# inspect markers for a group
View(as.data.frame(markers[[1]]))

# tSNE plots colored by different markers
sce$`HBA2` <- exprs(sce)['HBA2', ]
scater::plotTSNE(sce, colour_by = "HBA2", point_size = 2, point_alpha = 1, theme_size = 12)

