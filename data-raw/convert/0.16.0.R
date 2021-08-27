# create HDF5Array TENxMatrix for fast
library(SingleCellExperiment)

app_dir <- '/home/alex/patient_data/example'
app_dirs <- list.files('/srv/dseqr')
app_dirs <- setdiff(app_dirs,
                    c('indices', 'node_modules', 'pert_signature_dir',
                      'pert_query_dir', 'gs_dir', 'example_data.tar.gz', 'tx2gene'))


convert <- function(app_dir) {
    cat('==== Working on', app_dir, '====\n')

    sc_dir <- file.path(app_dir, 'single-cell')

    dataset_names <- list.files(sc_dir)
    dataset_dirs <- file.path(sc_dir, dataset_names)

    for (i in seq_along(dataset_dirs)) {
        dataset_name <- dataset_names[i]
        dataset_dir <- dataset_dirs[i]
        cat('working on', dataset_names[i], '...\n')

        dgclogs_path <- file.path(dataset_dir, 'dgclogs.qs')
        dgrlogs_path <- file.path(dataset_dir, 'dgrlogs.qs')
        if (!file.exists(dgclogs_path)) next()
        if (!file.exists(dgrlogs_path)) next()

        logs <- qs::qread(dgclogs_path)
        tlogs <- Matrix::t(logs)

        tlogs_path <- file.path(dataset_dir, 'tlogs.tenx')
        HDF5Array::writeTENxMatrix(tlogs, tlogs_path, group = 'mm10')
        unlink(dgrlogs_path)
    }
}

for  (app_dir in app_dirs) convert(app_dir)
