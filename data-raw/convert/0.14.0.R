library(SingleCellExperiment)

app_dir <- "/home/alex/patient_data/example"
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

    scseq_path <- file.path(dataset_dir, 'scseq.qs')

    if (!file.exists(scseq_path)) next()
    scseq <- qs::qread(scseq_path)

    # save as seperate parts
    dseqr:::split_save_scseq(scseq, dataset_dir)

    # delete previous
    unlink(scseq_path)
    remove <- list.files(dataset_dir,
                       '^lm_fit|^top_tables|^cluster_stats|plots|ambient.qs|^l1000_|^cmap_|^kegg|^go',
                       recursive = TRUE, full.names = TRUE, include.dirs = TRUE)
    unlink(remove, recursive = TRUE)

  }
}

for  (app_dir in app_dirs) convert(app_dir)
