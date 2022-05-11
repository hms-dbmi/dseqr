library(Seurat)
library(SingleCellExperiment)

get_seurat_object_list <- function() {
  # Tests for integration related fxns
  set.seed(42)
  pbmc_small <- suppressWarnings(UpdateSeuratObject(pbmc_small))

  # Setup test objects
  ref <- pbmc_small
  ref <- FindVariableFeatures(object = ref, verbose = FALSE, nfeatures = 100)
  ref$sample <- 'ref'

  query <- CreateSeuratObject(
    counts = GetAssayData(object = pbmc_small[['RNA']], slot = "counts") + rpois(n = ncol(pbmc_small), lambda = 1)
  )
  query$sample <- 'query'

  query.list <- list(ref, query)
  return(query.list)
}

mock_default_integration <- function() {

  query.list <- get_seurat_object_list()

  query.list <- lapply(X = query.list, FUN = NormalizeData, verbose = FALSE)
  query.list <- lapply(X = query.list, FUN = FindVariableFeatures, verbose = FALSE, nfeatures = 100)
  query.list <- lapply(X = query.list, FUN = ScaleData, verbose = FALSE)
  query.list <- suppressWarnings(lapply(X = query.list, FUN = RunPCA, verbose = FALSE, npcs = 20))

  anchors <- suppressMessages(suppressWarnings(FindIntegrationAnchors(object.list = query.list, k.filter = NA, verbose = FALSE)))
  int <- IntegrateData(anchorset = anchors, k.weight = 50, verbose = FALSE)

  # run the standard workflow for visualization and clustering
  int <- ScaleData(int, verbose = FALSE)
  int <- RunPCA(int, npcs = 30, verbose = FALSE)
  int <- suppressWarnings(RunUMAP(int, reduction = "pca", dims = 1:30))
  int <- FindNeighbors(int, reduction = "pca", dims = 1:30)
  int <- FindClusters(int, resolution = 0.5)

  # add annotation
  levels(Idents(int)) <- letters[seq_along(levels(Idents(int)))]

  return(int)
}

mock_sct_integration <- function() {
  query.list <- get_seurat_object_list()

  query.list <- lapply(X = query.list, FUN = SCTransform, verbose = FALSE)
  features <- SelectIntegrationFeatures(object.list = query.list, nfeatures = 100, verbose = FALSE)
  query.list <- PrepSCTIntegration(object.list = query.list, anchor.features = features, verbose = FALSE)

  anchors <- suppressWarnings(FindIntegrationAnchors(
    object.list = query.list,
    normalization.method = "SCT",
    anchor.features = features,
    k.filter = NA,
    verbose = FALSE))

  int <- IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT",
    k.weight = 50,
    verbose = FALSE)

  int <- RunPCA(int, npcs = 30, verbose = FALSE)
  int <- suppressWarnings(RunUMAP(int, reduction = "pca", dims = 1:30, verbose = FALSE))
  int <- FindNeighbors(int, reduction = "pca", dims = 1:30)
  int <- FindClusters(int, resolution = 0.5)

  # add cluster annotation
  levels(Idents(int)) <- letters[seq_along(levels(Idents(int)))]
  return(int)
}

test_that("a unisample Seurat object can be imported", {

  # setup
  dataset_name <- 'test'
  sc_dir <- file.path(tempdir(), 'single-cell')
  uploaded_data_dir <- file.path(sc_dir, dataset_name)
  dir.create(uploaded_data_dir, recursive = TRUE)

  tx2gene_dir <- file.path(tempdir(), 'tx2gene')
  dir.create(tx2gene_dir)

  pbmc_raw <- read.table(
    file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
    as.is = TRUE
  )

  pbmc_small <- Seurat::CreateSeuratObject(counts = pbmc_raw)
  qs::qsave(pbmc_small, file.path(uploaded_data_dir, 'test.qs'))

  expect_error(
    suppressWarnings(import_scseq(dataset_name, uploaded_data_dir, sc_dir, tx2gene_dir, species = 'Homo sapiens')),
    NA
  )

  # cleanup
  unlink(c(sc_dir, uploaded_data_dir, tx2gene_dir), recursive = TRUE)
})


test_that("a multisample Seurat object can be imported", {

  int <- mock_default_integration()

  # setup
  dataset_name <- 'test'
  sc_dir <- file.path(tempdir(), 'single-cell')
  uploaded_data_dir <- file.path(sc_dir, dataset_name)
  dir.create(uploaded_data_dir, recursive = TRUE)

  tx2gene_dir <- file.path(tempdir(), 'tx2gene')
  dir.create(tx2gene_dir)

  qs::qsave(int, file.path(uploaded_data_dir, 'test.qs'))

  expect_error(
    suppressWarnings(import_scseq(dataset_name, uploaded_data_dir, sc_dir, tx2gene_dir, species = 'Homo sapiens')),
    NA
  )

  # load imported scseq and attach annotation
  dataset_dir <- file.path(sc_dir, dataset_name)
  scseq <- load_scseq_qs(dataset_dir, with_logs = TRUE, with_counts = TRUE)
  annot <- qs::qread(file.path(dataset_dir, 'snn1', 'annot.qs'))
  levels(scseq$cluster) <- annot

  # have UMAP, corrected, counts, and logcounts
  expect_setequal(c('corrected', 'UMAP'), reducedDimNames(scseq))
  expect_setequal(c('logcounts', 'counts'), assayNames(scseq))

  # UMAP preserved
  umap.imported <- reducedDim(scseq, 'UMAP')
  umap.original <- int[["umap"]]@cell.embeddings
  expect_equal(umap.imported, umap.original)

  # PCA preserved
  pca.imported <- reducedDim(scseq, 'corrected')
  pca.original <- int[['pca']]@cell.embeddings
  expect_equal(pca.imported, pca.original)

  # counts preserved
  expect_equal(int[['RNA']]@counts, counts(scseq))

  # have clusters and annotation
  expect_length(levels(scseq$cluster), length(levels(int$seurat_clusters)))

  # clusters preserved
  expect_equal(table(Idents(int)), table(scseq$cluster))

  # added samples to batch
  expect_equal(unname(int$sample), scseq$batch)

  # have same rows as RNA assay (integrated just has variable features)
  expect_equal(row.names(scseq), row.names(int[['RNA']]))

  # cleanup
  unlink(c(sc_dir, uploaded_data_dir, tx2gene_dir), recursive = TRUE)
})


test_that("a multisample Seurat object without samples specified can be imported", {

  # remove sample column
  int <- mock_default_integration()
  int$sample <- NULL

  # setup
  dataset_name <- 'test'
  sc_dir <- file.path(tempdir(), 'single-cell')
  uploaded_data_dir <- file.path(sc_dir, dataset_name)
  dir.create(uploaded_data_dir, recursive = TRUE)

  tx2gene_dir <- file.path(tempdir(), 'tx2gene')
  dir.create(tx2gene_dir)

  qs::qsave(int, file.path(uploaded_data_dir, 'test.qs'))

  expect_error(
    suppressWarnings(import_scseq(dataset_name, uploaded_data_dir, sc_dir, tx2gene_dir, species = 'Homo sapiens')),
    NA
  )

  # load imported scseq and attach annotation
  dataset_dir <- file.path(sc_dir, dataset_name)
  scseq <- load_scseq_qs(dataset_dir, with_logs = TRUE, with_counts = TRUE)
  annot <- qs::qread(file.path(dataset_dir, 'snn1', 'annot.qs'))
  levels(scseq$cluster) <- annot

  # have UMAP, corrected, counts, and logcounts
  expect_setequal(c('corrected', 'UMAP'), reducedDimNames(scseq))
  expect_setequal(c('logcounts', 'counts'), assayNames(scseq))

  # UMAP preserved
  umap.imported <- reducedDim(scseq, 'UMAP')
  umap.original <- int[["umap"]]@cell.embeddings
  expect_equal(umap.imported, umap.original)

  # PCA preserved
  pca.imported <- reducedDim(scseq, 'corrected')
  pca.original <- int[['pca']]@cell.embeddings
  expect_equal(pca.imported, pca.original)

  # counts preserved
  expect_equal(int[['RNA']]@counts, counts(scseq))

  # have clusters and annotation
  expect_length(levels(scseq$cluster), length(levels(int$seurat_clusters)))

  # clusters preserved
  expect_equal(table(Idents(int)), table(scseq$cluster))

  # have same rows as RNA assay (integrated just has variable features)
  expect_equal(row.names(scseq), row.names(int[['RNA']]))

  # cleanup
  unlink(c(sc_dir, uploaded_data_dir, tx2gene_dir), recursive = TRUE)
})


test_that("a multisample Seurat object integrated with SCTransform can be imported", {

  # run SCTransform integration
  int <- mock_sct_integration()

  # setup
  dataset_name <- 'test'
  sc_dir <- file.path(tempdir(), 'single-cell')
  uploaded_data_dir <- file.path(sc_dir, dataset_name)
  dir.create(uploaded_data_dir, recursive = TRUE)

  tx2gene_dir <- file.path(tempdir(), 'tx2gene')
  dir.create(tx2gene_dir)

  qs::qsave(int, file.path(uploaded_data_dir, 'test.qs'))

  expect_error(
    suppressWarnings(import_scseq(dataset_name, uploaded_data_dir, sc_dir, tx2gene_dir, species = 'Homo sapiens')),
    NA
  )

  # load imported scseq and attach annotation
  dataset_dir <- file.path(sc_dir, dataset_name)
  scseq <- load_scseq_qs(dataset_dir, with_logs = TRUE, with_counts = TRUE)
  annot <- qs::qread(file.path(dataset_dir, 'snn1', 'annot.qs'))
  levels(scseq$cluster) <- annot

  # have UMAP, corrected, counts, and logcounts
  expect_setequal(c('corrected', 'UMAP'), reducedDimNames(scseq))
  expect_setequal(c('logcounts', 'counts'), assayNames(scseq))

  # UMAP preserved
  umap.imported <- reducedDim(scseq, 'UMAP')
  umap.original <- int[["umap"]]@cell.embeddings
  expect_equal(umap.imported, umap.original)

  # PCA preserved
  pca.imported <- reducedDim(scseq, 'corrected')
  pca.original <- int[['pca']]@cell.embeddings
  expect_equal(pca.imported, pca.original)

  # counts preserved
  expect_equal(int[['RNA']]@counts, counts(scseq))

  # have clusters and annotation
  expect_length(levels(scseq$cluster), length(levels(int$seurat_clusters)))

  # clusters preserved
  expect_equal(table(Idents(int)), table(scseq$cluster))

  # added samples to batch
  expect_equal(unname(int$sample), scseq$batch)

  # have same rows as RNA assay (integrated just has variable features)
  expect_equal(row.names(scseq), row.names(int[['RNA']]))

  # cleanup
  unlink(c(sc_dir, uploaded_data_dir, tx2gene_dir), recursive = TRUE)
})
