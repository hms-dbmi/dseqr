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
sce <- drugseqr::load_scseq(data_dir)

# alevin does initial knee-based whitelisting then refines based on mtRNA, rRNA, etc
# see what things look like before/after final whitelisting
sce <- scater::calculateQCMetrics(sce, feature_controls=list(mito = which(row.names(sce) %in% sce@metadata$mrna),
                                                             ribo = which(row.names(sce) %in% sce@metadata$rrna)))

sce_df <- tibble(log10_total_counts = sce$log10_total_counts,
                 log10_total_features_by_counts = sce$log10_total_features_by_counts,
                 pct_counts_ribo = sce$pct_counts_ribo,
                 pct_counts_mito = sce$pct_counts_mito,
                 whitelist = sce$whitelist)

ribo <- ggplot(sce_df, aes(x=pct_counts_ribo, fill=whitelist)) +
  theme_minimal() +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  geom_histogram(position="identity", colour="black", alpha = 0.5, size = 0.05) +
  xlab('Proportion of reads in ribosomal genes') +
  theme(legend.position = "none")

mito <- ggplot(sce_df, aes(x=pct_counts_mito, fill=whitelist)) +
  theme_minimal() +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  geom_histogram(position="identity", colour="black", alpha = 0.5, size = 0.05) +
  xlab('Proportion of reads in mitochondrial genes') +
  theme(legend.position = "none") +
  ylab('')

ncounts <- ggplot(sce_df, aes(x=log10_total_counts, fill=whitelist)) +
  theme_minimal() +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  geom_histogram(position="identity", colour="black", alpha = 0.5, size = 0.05) +
  xlab('Log-total UMI count') +
  theme(legend.position = "none")

nexp <- ggplot(sce_df, aes(x=log10_total_features_by_counts, fill=whitelist)) +
  theme_minimal() +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  geom_histogram(position="identity", colour="black", alpha = 0.5, size = 0.05) +
  xlab('Log-total number of expressed features') +
  theme(legend.position = "none") +
  ylab('')

legend <- get_legend(nexp + theme(legend.position="top"))
plot_grid(legend, NULL, ncounts, nexp, ribo, mito, nrow = 3, rel_heights = c(0.2, 1, 1))

# subset by alevin whitelist
sce.full <- sce
sce <- sce[, sce$whitelist]

# look at most highly expressed genes
# expect dominated by mitochondrial/ribosomal
ave <- scater::calcAverage(sce)
rowData(sce)$AveCount <- ave
hist(log10(ave), col="grey80", breaks = 100)
scater::plotHighestExprs(sce)

# normalizing for cells-specific biases ----

sce <- norm_scseq(sce)

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
set.seed(1000)
sce <- denoisePCA(sce, technical=new.trend, BSPARAM=IrlbaParam())
ncol(reducedDim(sce, "PCA"))

plot(attr(reducedDim(sce), "percentVar"), xlab="PC",
     ylab="Proportion of variance explained")
abline(v=ncol(reducedDim(sce, "PCA")), lty=2, col="red")

plotPCA(sce, ncomponents=3, colour_by="log10_total_features_by_counts")

# tSNE version
set.seed(1000)
sce <- runTSNE(sce, use_dimred="PCA", perplexity=30)
plotTSNE(sce, colour_by="log10_total_features_by_counts")

# clustering with graphs ----
snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
clusters <- igraph::cluster_walktrap(snn.gr)
sce$Cluster <- factor(clusters$membership)
table(sce$Cluster)

cluster.mod <- clusterModularity(snn.gr, sce$Cluster, get.values=TRUE)
log.ratio <- log2(cluster.mod$observed/cluster.mod$expected + 1)

library(pheatmap)
pheatmap(log.ratio, cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("white", "blue"))(100))

plotTSNE(sce, colour_by="Cluster")

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

pass_qc <- colnames(sce.full) %in% whitelist

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
dbl.out <- doubletCluster(sce, sce$Cluster)
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
