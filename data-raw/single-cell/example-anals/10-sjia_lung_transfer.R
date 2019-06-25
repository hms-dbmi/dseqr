library(drugseqr)
library(Seurat)

# use lung as reference
lung_anal <- readRDS('data-raw/single-cell/example-data/test_lung_anal.rds')
levels(lung_anal$scseq$seurat_clusters) <- lung_anal$annot
names(lung_anal$markers) <- lung_anal$annot
Idents(lung_anal$scseq) <- lung_anal$scseq$seurat_clusters

explore_scseq_clusters(lung_anal$scseq, lung_anal$markers)

# use bone marrow as query -----
bm_anal <- readRDS('data-raw/single-cell/example-data/bm_anal.rds')
levels(bm_anal$scseq$seurat_clusters) <- bm_anal$annot
names(bm_anal$markers) <- bm_anal$annot
Idents(bm_anal$scseq) <- bm_anal$scseq$seurat_clusters

explore_scseq_clusters(bm_anal$scseq, bm_anal$markers)

bm_preds <- transfer_labels(lung_anal$scseq, bm_anal$scseq)
table(bm_preds$predicted.id)

# explore bone marrow with new labels
bm_anal$scseq$seurat_clusters <- Idents(bm_anal$scseq) <- factor(preds$predicted.id)
new_markers <- get_scseq_markers(bm_anal$scseq)

explore_scseq_clusters(bm_anal$scseq, new_markers)


# use combined blood samples as query -----
# blood_anal <- readRDS('data-raw/single-cell/example-data/sjia_blood1.rds')
blood_anal <- readRDS('data-raw/single-cell/example-data/sjia_blood2.rds')

blood_preds <- transfer_labels(lung_anal$scseq, blood_anal$scseq)
table(blood_preds$predicted.id)

# Macrophages#2      NK-cells       T-cells
#            31            44            14

# explore bone marrow with new labels
blood_anal$scseq$seurat_clusters <- Idents(blood_anal$scseq) <- factor(blood_preds$predicted.id)
new_markers <- get_scseq_markers(blood_anal$scseq)

explore_scseq_clusters(blood_anal$scseq, new_markers)

