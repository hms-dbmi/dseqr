library(Seurat)
library(drugseqr)
library(dplyr)

data_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'
srt <- load_scseq(data_dir)

# subset by alevin whitelist
srt <- srt[, srt$whitelist]


# store mitochondrial percentage
qcgenes <- load_scseq_qcgenes()
srt <- PercentageFeatureSet(srt, pattern = paste0('^', scqc_genes$mrna, '$', collapse = '|'), col.name = "percent.mt")

# run sctransform
srt <- SCTransform(srt, vars.to.regress = "percent.mt", verbose = FALSE)
srt <- SCTransform(srt, verbose = FALSE)


# These are now standard steps in the Seurat workflow for visualization and clustering
srt <- RunPCA(srt, verbose = FALSE)
srt <- RunUMAP(srt, dims = 1:30, verbose = FALSE)

srt <- FindNeighbors(srt, dims = 1:30, verbose = FALSE)
srt <- FindClusters(srt, verbose = FALSE)
DimPlot(srt, label = TRUE) + NoLegend()



markers <- FindAllMarkers(srt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

VlnPlot(srt, features = c("C1QC", "APOC1"))

FeaturePlot(srt, features = top_markers$gene)
