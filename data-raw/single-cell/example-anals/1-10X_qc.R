# adapted from: https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/tenx.html

library(scater)
library(dplyr)
library(DropletUtils)
library(scran)
library(BiocSingular)

# load/annotate data from alevin ----

# folder with 10X fastq files
data_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'

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
  dplyr::select(-tx_id) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise_all(unique)

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

par(mfrow=c(1,1))
# subset by alevin whitelist
sce.full <- sce
sce <- sce[, colnames(sce) %in% whitelist]

# look at most highly expressed genes
# expect dominated by mitochondrial/ribosomal
ave <- calcAverage(sce)
rowData(sce)$AveCount <- ave
hist(log10(ave), col="grey80", breaks = 100)
plotHighestExprs(sce)

# normalizing for cells-specific biases ----
set.seed(1000)
clusters <- quickCluster(sce, use.ranks=FALSE, BSPARAM=IrlbaParam())
table(clusters)

sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters)

# sizefactors should correlated well with total counts
# indicates capture efficiency/sequencing depth are the major biases
plot(sce$total_counts, sizeFactors(sce), log="xy")

# normalize based on size-factors
sce <- normalize(sce)

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
# based on https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/qc.html#32_in_the_pbmc_data_set

pass_qc <- colnames(sce.full) %in% whitelist

lost <- calcAverage(counts(sce.full)[, !pass_qc])
kept <- calcAverage(counts(sce.full)[, pass_qc])

# Avoid loss of points where either average is zero.
capped.lost <- pmax(lost, min(lost[lost>0]))
capped.kept <- pmax(kept, min(kept[kept>0]))

# look for shift to bottom right (genes that are high in discarded cells, low in retained cells)
plot(capped.lost, capped.kept, xlab="Average count (discarded)",
     ylab="Average count (retained)", log="xy", pch=16)
is.spike <- isSpike(sce.full.416b)
points(capped.lost[is.spike], capped.kept[is.spike], col="red", pch=16)
is.mito <- rowData(sce.full.416b)$is_feature_control_Mt
points(capped.lost[is.mito], capped.kept[is.mito], col="dodgerblue", pch=16)
