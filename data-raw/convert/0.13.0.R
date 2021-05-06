app_dir <- "/home/alex/patient_data/sjia"
sc_dir <- file.path(app_dir, 'single-cell')

integrated <- qs::qread(file.path(sc_dir, 'integrated.qs'))

dataset_dirs <- file.path(sc_dir, integrated)
deleted <- c()


for (i in seq_along(dataset_dirs)) {
  dataset_name <- integrated[i]
  dataset_dir <- dataset_dirs[i]
  cat('working on', integrated[i], '...\n')
  scseq_path <- file.path(dataset_dir, 'scseq.qs')

  if (!file.exists(scseq_path)) {
    deleted <- c(deleted, i)
    next()
  }

  # update test/ctrl from args json
  args <- dseqr:::load_args(sc_dir, dataset_name)
  if (is.null(args$dataset_names)) {
    args$dataset_names <- c(args$test, args$ctrl)
    args$test <- args$ctrl <- NULL
    dseqr:::save_scseq_args(args, dataset_name, sc_dir)
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
  qs::qsave(combined, scseq_path, preset = 'fast')

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
}

# remove deleted
integrated <- integrated[-deleted]
qs::qsave(integrated, file.path(sc_dir, 'integrated.qs'))

# change to presto markers
# and replace looms with fast qs files
dataset_names <- list.files(sc_dir)
dataset_names <- dataset_names[!grepl('.qs$', dataset_names)]


for (dataset_name in dataset_names) {
  cat('working on', dataset_name, '...\n')
  dataset_dir <- file.path(sc_dir, dataset_name)
  scseq_path <- file.path(dataset_dir, 'scseq.qs')
  if (!file.exists(scseq_path)) next()

  scseq <- qs::qread(scseq_path)
  qs::qsave(scseq, scseq_path, preset = 'fast')
  unlink(file.path(dataset_dir, 'scle.loom'))

  snn_names <- list.files(dataset_dir, '^snn\\d')

  marker_files <- list.files(dataset_dir, '^markers_\\d+.qs$', recursive = TRUE, full.names = TRUE)
  unlink(marker_files)

  for (snn_name in snn_names) {
    snn_dir <- file.path(dataset_dir, snn_name)

    # remove uneeded
    unlink(file.path(snn_dir, 'tests'), recursive = TRUE)
    unlink(list.files(snn_dir, 'l1000_.+?.qs', full.names = TRUE))
    unlink(list.files(snn_dir, 'cmap_.+?.qs', full.names = TRUE))
    unlink(list.files(snn_dir, 'kegg_.+?.qs', full.names = TRUE))
    unlink(list.files(snn_dir, 'kegga_.+?.qs', full.names = TRUE))
    unlink(list.files(snn_dir, 'go_.+?.qs', full.names = TRUE))
    unlink(list.files(snn_dir, 'goana_.+?.qs', full.names = TRUE))
    unlink(file.path(snn_dir, 'applied.qs'))
    unlink(file.path(snn_dir, 'top_markers.qs'))
    unlink(file.path(snn_dir, 'pbulk_esets.qs'))
  }
}
