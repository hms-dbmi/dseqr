library(Seurat)
library(SingleCellExperiment)

# read in exported data from dseqr
scseq <- qs::qread('scseq.qs')
colnames(scseq) <- make.unique(colnames(scseq))

sdata <- as.Seurat(scseq)

# add cluster to appropriate column
sdata$seurat_clusters <- sdata$cluster
Idents(sdata) <- 'seurat_clusters'
sdata$samples <- sdata$batch

# add groupings (will be auto-detected by Cellenics)
meta <- scseq@metadata$meta
sdata$Group <- meta[sdata$batch, 'group']

# add variable features and pca (needed for Cellenics)
sdata <- SeuratObject::RenameAssays(sdata, 'originalexp' = 'RNA')

sdata <- FindVariableFeatures(sdata)
sdata <- ScaleData(sdata)
sdata <- RunPCA(sdata)


# save for Cellenics
saveRDS(sdata, 'scseq_cellenics.rds')
