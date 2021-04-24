app_dir <- "/home/alex/patient_data/sjia"
sc_dir <- file.path(app_dir, 'single-cell')

integrated <- qs::qread(file.path(sc_dir, 'integrated.qs'))

dataset_dirs <- file.path(sc_dir, integrated)
deleted <- c()

for (i in seq_along(dataset_dirs)) {
  dataset_dir <- dataset_dirs[i]
  cat('working on', integrated[i], '...\n')
  scseq_path <- file.path(dataset_dir, 'scseq.qs')

  if (!file.exists(scseq_path)) {
    deleted <- c(deleted, i)
    next()
  }

  combined <- qs::qread(scseq_path)
  dataset_names <- unique(combined$batch)

  # add ambient info for each sample
  scseqs <- dseqr:::load_scseq_subsets(dataset_names, sc_dir)
  combined <- dseqr:::add_combined_ambience(combined, scseqs)

  # remove 'test' and 'ctrl' ambient
  SummarizedExperiment::rowData(combined)$ctrl_ambient <- NULL
  SummarizedExperiment::rowData(combined)$test_ambient <- NULL

  # create meta & prev_groups from test/ctrl orig.ident
  if (all(levels(combined$orig.ident) %in% c('test', 'ctrl'))) {
    samples <- unique(combined$batch)
    test <- unique(combined$batch[combined$orig.ident == 'test'])
    ctrl <- unique(combined$batch[combined$orig.ident == 'ctrl'])

    group <- rep(NA, length(samples))
    group[samples %in% test] <- 'test'
    group[samples %in% ctrl] <- 'ctrl'

    meta <- data.frame(group,
                       pair = NA,
                       row.names = samples,
                       check.names = FALSE, stringsAsFactors = FALSE)

    prev_groups <- c('test', 'ctrl')
    qs::qsave(meta, file.path(dataset_dir, 'meta.qs'))
    qs::qsave(prev_groups, file.path(dataset_dir, 'prev_groups.qs'))
  }

  # use sample for orig.ident
  combined$orig.ident <- factor(combined$batch)


  # overwrite loom/qs objects
  dseqr::save_scle(combined, dataset_dir, overwrite = TRUE)
  qs::qsave(combined, scseq_path)

  # remove unnecessary
  unlink(file.path(dataset_dir, c('ambient.qs', 'has_replicates.qs')))
  snn_dirs <- list.files(dataset_dir, 'snn\\d', full.names = TRUE)
  unlink(file.path(snn_dirs, 'tests'), recursive = TRUE)
  unlink(file.path(snn_dirs, 'plots'), recursive = TRUE)
  unlink(file.path(snn_dirs, 'lm_fit_0svs.qs'), recursive = TRUE)
  unlink(file.path(snn_dirs, 'has_replicates.qs'), recursive = TRUE)
  unlink(file.path(snn_dirs, 'top_tables.qs'), recursive = TRUE)
  unlink(file.path(snn_dirs, 'top_markers.qs'), recursive = TRUE)
  unlink(file.path(snn_dirs, 'cluster_stats.qs'), recursive = TRUE)

  # calculate ambient
  resoln <- dseqr:::load_resoln(dataset_dir)
  snn_dir <- file.path(dataset_dir, resoln)
  cat(' --- working on snn_dir', basename(snn_dir), '..\n')

  ambient_path <- file.path(snn_dir, 'ambient.qs')
  summed <- qs::qread(file.path(snn_dir, 'summed.qs'))

  clusters <- levels(summed$cluster)
  ambience <- dseqr:::get_ambience(combined)

  amb <- list()
  for (i in seq_along(clusters)) {
    clus <- clusters[i]
    cat('working on cluster', clus, '..\n')
    amb[[clus]] <- dseqr:::calc_cluster_ambience(summed, ambience, clus)
  }
  qs::qsave(amb, ambient_path)
}

integrated <- integrated[-deleted]

