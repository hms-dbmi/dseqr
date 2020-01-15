# conversion of legacy Seurat to new SingleCellExperiment stuff
library(drugseqr)
setwd("~/Documents/Batcave/zaklab/drugseqr/data-raw/patient_data/sjia/single-cell")
dirs <- setdiff(list.files(), c('integrated.rds', '10X-data'))
sc_dir <- '.'

for (dataset_name in dirs[-1]) {
  cat('working on', dataset_name, '...\n')
  scseq_path <- file.path(dataset_name, 'scseq.rds')
  sce <- readRDS(scseq_path)

  if (class(sce) == 'Seurat') {
    # make copy of saved scseq in case mess up
    scseq_path <- scseq_part_path(sc_dir, dataset_name, 'scseq')
    scseq_copy <- gsub('scseq.rds$', 'scseq_copy.rds', scseq_path)
    file.copy(scseq_path, scseq_copy)

    # restore counts slot
    sce[['RNA']]@counts <- sce[['RNA']]@data

    Seurat::DefaultAssay(sce) <- 'RNA'
    sce <- Seurat::as.SingleCellExperiment(sce)
  }

  cat('normalizing...\n')
  # sce <- normalize_scseq(sce)
  cat('clustering...\n')
  # sce <- add_hvgs(sce)
  sce <- add_scseq_clusters(sce)
  cat('reducing...\n')
  sce <- run_tsne(sce)

  cat('getting markers...\n')
  wilcox_tests <- pairwise_wilcox(sce)
  markers <- get_scseq_markers(wilcox_tests)

  # top markers are for SingleR
  top_markers <- scran::getTopMarkers(wilcox_tests$statistics, wilcox_tests$pairs)
  annot <- names(markers)

  cat('saving...\n')
  scseq_data <- list(scseq = sce, markers = markers, top_markers = top_markers, annot = annot)

  data_dir <- file.path(sc_dir, dataset_name)
  unlink(data_dir, recursive = TRUE)
  dir.create(data_dir)

  save_scseq_data(scseq_data, dataset_name, sc_dir)
}
