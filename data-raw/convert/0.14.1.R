library(SingleCellExperiment)

app_dir <- "/home/alex/patient_data/sjia"
app_dirs <- list.files('/srv/dseqr')
app_dirs <- setdiff(app_dirs,
                    c('indices', 'node_modules', 'pert_signature_dir',
                      'pert_query_dir', 'example_data.tar.gz'))

convert <- function(app_dir) {
  cat('==== Working on', app_dir, '====\n')

  sc_dir <- file.path(app_dir, 'single-cell')

  dataset_names <- list.files(sc_dir)
  dataset_dirs <- file.path(sc_dir, dataset_names)

  for (i in seq_along(dataset_dirs)) {
    dataset_name <- dataset_names[i]
    dataset_dir <- dataset_dirs[i]
    cat('working on', dataset_names[i], '...\n')

    scseq_path <- file.path(dataset_dir, 'shell.qs')
    logs_path <- file.path(dataset_dir, 'dgclogs.qs')

    if (!file.exists(scseq_path)) next()
    scseq <- qs::qread(scseq_path)
    logcounts(scseq) <- qs::qread(logs_path)
    decs <- scran::modelGeneVar(scseq, block=scseq$batch)
    SummarizedExperiment::rowData(scseq)$bio <- decs$bio

    logcounts(scseq) <- NULL

    # orig.ident (test/ctrl) can be updated
    scseq$orig.ident <- scseq$orig.cluster <- scseq$orig.resoln <- NULL

    qs::qsave(scseq, scseq_path, preset = 'fast')
  }
}

for  (app_dir in app_dirs) convert(app_dir)
