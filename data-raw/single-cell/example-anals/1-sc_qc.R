# first several parts adapted from: https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/tenx.html

library(scater)
library(dplyr)
library(DropletUtils)
library(scran)
library(BiocSingular)
library(ggplot2)
library(cowplot)

# load/annotate data from alevin ----

# folder with 10X fastq files
data_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'

# this is LONG RUNNING
# drugseqr::run_alevin(data_dir)

# import alevin quants
sce <- load_scseq(data_dir)

# run qc
sce <- qc_scseq(sce)

# alevin does initial knee-based whitelisting then refines based on mtRNA, rRNA, etc
# see what things look like before/after final whitelisting
hist_scseq_whitelist(sce)

# first normalize and stabilize
sce.full <- sce
sce.full <- norm_scseq(sce.full)
sce.full <- stabilize_scseq(sce.full)

# show tSNE based on whitelist and whitelist metrics
tsne_scseq_whitelist(sce.full)

# subset by alevin whitelist and re-norm/stabilize
sce <- sce[, sce$whitelist]
sce <- norm_scseq(sce)
sce <- stabilize_scseq(sce)

# look at most highly expressed genes ----

# expect dominated by mitochondrial/ribosomal
ave <- scater::calcAverage(sce)
rowData(sce)$AveCount <- ave
hist(log10(ave), col="grey80", breaks = 100)
scater::plotHighestExprs(sce)

# sizefactors should correlated well with total counts
# indicates capture efficiency/sequencing depth are the major biases
plot(sce$total_counts, sizeFactors(sce), log="xy")

# variance stabilization ----
# variance goes up with mean in count-based data
# note that this is same goal as limma::voom for bulk rna-seq data

# assume technical variation is poisson
new.trend <- makeTechTrend(x=sce)

# estimate observed variances
fit <- trendVar(sce, use.spikes=FALSE, loess.args=list(span=0.05))
plot(fit$mean, fit$var, pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE)
curve(new.trend(x), col="red", add=TRUE)

# decompose variance using poisson as technical variance
fit$trend <- new.trend # overwrite trend.
dec <- decomposeVar(fit=fit) # use per-gene variance estimates in 'fit'.
top.dec <- dec[order(dec$bio, decreasing=TRUE),]
head(top.dec)

# genes with the largest biological component in variance
plotExpression(sce, features=rownames(top.dec)[1:15])

# dimensionality reduction -----

# tSNE version
set.seed(1000)
sce <- runTSNE(sce, use_dimred="PCA")

# clustering with graphs ----
snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
clusters <- igraph::cluster_walktrap(snn.gr)
sce$cluster <- factor(clusters$membership)
table(sce$cluster)

# using original tSNE plot show clusters of whitelisted cells
sce.full$cluster <- NA
sce.full$cluster[sce.full$whitelist] <- sce$cluster
sce.full$cluster <- as.factor(sce.full$cluster)
plotTSNE(sce.full[, sce.full$whitelist], colour_by="cluster")

# this re-does the tSNE so can't compare to original anymore
plotTSNE(sce, colour_by="cluster")

# Marker gene detection of clusters -----
# only upregulated as more useful for positive id of cell type
markers <- findMarkers(sce, clusters=sce$Cluster, direction="up")

# look at e.g. cluster 1 with genes with top 10 minimum ranks
marker.set <- markers[["1"]]
chosen <- rownames(marker.set)[marker.set$Top <= 10]
plotHeatmap(sce, features=chosen, exprs_values="logcounts",
            zlim=5, center=TRUE, symmetric=TRUE, cluster_cols=FALSE,
            colour_columns_by="Cluster", columns=order(sce$Cluster),
            show_colnames=FALSE)


# check for discarded cell types -----
# based on https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/qc.html
# doesn't automate well, just a visual check
# personal note: may be useful to not conclude that one sample has cluster and other doesn't

pass_qc <- colnames(sce.full) %in% sce.full$whitelist

lost <- calcAverage(counts(sce.full)[, !pass_qc])
kept <- calcAverage(counts(sce.full)[, pass_qc])

# Avoid loss of points where either average is zero.
capped.lost <- pmax(lost, min(lost[lost>0]))
capped.kept <- pmax(kept, min(kept[kept>0]))

# look for shift to bottom right (genes that are high in discarded cells, low in retained cells)
plot(capped.lost, capped.kept, xlab="Average count (discarded)",
     ylab="Average count (retained)", log="xy", pch=16)

# check for doublets -----
# based on https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/doublets.html
# doesn't automate well, good to check that not relying on possibly double population for downstream analysis

# relies on good clusters
dbl.out <- doubletCluster(sce, sce$cluster)
dbl.out

dbl.markers <- markers[["1"]]
chosen <- rownames(dbl.markers)[dbl.markers$Top <= 10]
plotHeatmap(sce, columns=order(sce$Cluster), colour_columns_by="Cluster",
            features=chosen, cluster_cols=FALSE, center=TRUE, symmetric=TRUE,
            zlim=c(-5, 5), show_colnames=FALSE)


# simulates doublets from data
set.seed(100)
dbl.dens <- doubletCells(sce, BSPARAM=IrlbaParam())
summary(dbl.dens)
sce$DoubletScore <- dbl.dens
plotTSNE(sce, colour_by="DoubletScore")

# can check doublet score for each cluster
# is the worst offender also the worst identified by doubletCluster?
plotColData(sce, x="Cluster", y="DoubletScore", colour_by="Cluster")
