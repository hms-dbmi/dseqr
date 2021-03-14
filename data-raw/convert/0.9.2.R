# this script moves everything to leiden clustering
# it is necessary to transition data created prior to version 0.9.2
# it uses dseqr version 0.9.2
library(SingleCellExperiment)

# fill in path to app data
app_dir <- "/home/alex/patient_data/example"


sc_dir <- file.path(app_dir, 'single-cell')
dataset_names <- list.dirs(sc_dir, full.names = FALSE, recursive = FALSE)

resoln = 1
cluster_alg = 'leiden'
npcs = 30

transition_dataset <- function(dataset_name) {
  cat('Working on:', dataset_name, '...\n')
  dataset_dir <- file.path(sc_dir, dataset_name)
  unlink(list.files(dataset_dir, 'snn\\d', full.names = TRUE), recursive = TRUE)

  scseq <- dseqr:::load_scseq_qs(dataset_dir)
  scseq@metadata$npcs <- 30
  dimred <- ifelse('corrected' %in% reducedDimNames(scseq), 'corrected', 'PCA')
  scseq <- dseqr::run_tsne(scseq, dimred)

  # get snn graph
  snn_graph <- dseqr:::get_snn_graph(scseq, npcs)

  saveRDS(resoln, file.path(dataset_dir, 'resoln.rds'))
  qs::qsave(snn_graph, file.path(dataset_dir, 'snn_graph.qs'))

  # update saved
  ref_cluster <- scseq$cluster
  if (is.null(ref_cluster)) {
    ref_cluster <- as.character(scseq$seurat_clusters)
    ref_cluster <- factor(as.numeric(ref_cluster)+1)
  }
  scseq$cluster <- NULL
  scseq$seurat_clusters <- NULL
  qs::qsave(scseq, file.path(dataset_dir, 'scseq.qs'))
  dseqr:::save_scle(scseq, file.path(sc_dir, dataset_name))

  # run what depends on resolution
  query_cluster <- dseqr:::get_clusters(snn_graph, cluster_alg, resoln)
  scseq$cluster <- query_cluster

  dseqr:::run_post_cluster(scseq, dataset_name, sc_dir, resoln)

  # transfer labels
  tab <- table(assigned = ref_cluster, cluster = query_cluster)
  pred <- row.names(tab)[apply(tab, 2, which.max)]
  annot <- dseqr:::get_pred_annot(pred, dataset_name, file.path(dataset_name, 'snn1'), sc_dir)
  annot_nums <- as.character(seq_along(annot))

  # keep ordered nums where prediction is numeric
  is.num <- !is.na(as.numeric(gsub('_\\d+$', '', annot)))
  annot[is.num] <- annot_nums[is.num]

  saveRDS(annot, file.path(dataset_dir, 'snn1', 'annot.rds'))
}

remove_old <- function(dataset_name) {
  cat('Working on:', dataset_name, '...\n')
  dataset_dir <- file.path(sc_dir, dataset_name)

  to.remove <- list.files(
    dataset_dir,
    paste0(
      'aggr_ref.rds|tests|cluster_stats.rds|scseq_sample.rds|l1000_.+.rds|cmap_.+?.rds|',
      'kegg_.+?.rds|goana_.+?.rds|preds.rds|annot.rds|kegga_.+?.rds|go_.+?.rds|top_tables.rds|',
      'pbulk_esets.rds|top_markers.rds|lm_fit_.+?.rds|markers_.+?.rds|summed.rds|plots'
    ))

  unlink(file.path(dataset_dir, to.remove), recursive = TRUE)

}

for (dataset_name in dataset_names) {
  try(transition_dataset(dataset_name))
}

for (dataset_name in dataset_names) {
  remove_old(dataset_name)
}
