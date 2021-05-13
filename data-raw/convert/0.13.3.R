# remove hardcoded is_test and is_ctrl

app_dir <- "/home/alex/patient_data/sjia"
app_dirs <- list.files('/srv/dseqr')
app_dirs <- setdiff(app_dirs,
                    c('indices', 'node_modules', 'pert_signature_dir',
                      'pert_query_dir', 'example_data.tar.gz'))

convert <- function(app_dir) {
  cat('==== Working on', app_dir, '====\n')

  sc_dir <- file.path(app_dir, 'single-cell')

  integrated <- dseqr:::qread.safe(file.path(sc_dir, 'integrated.qs'))
  if (!length(integrated)) next()

  dataset_dirs <- file.path(sc_dir, integrated)

  for (i in seq_along(dataset_dirs)) {
    dataset_name <- integrated[i]
    dataset_dir <- dataset_dirs[i]
    cat('working on', integrated[i], '...\n')
    scseq_path <- file.path(dataset_dir, 'scseq.qs')

    if (!file.exists(scseq_path)) next()

    combined <- qs::qread(scseq_path)

    # remove 'is_test' and 'is_ctrl'
    combined$is_test <- combined$is_ctrl <- NULL

    # overwrite qs
    qs::qsave(combined, scseq_path, preset = 'fast')
  }
}
