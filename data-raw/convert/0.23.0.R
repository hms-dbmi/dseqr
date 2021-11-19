# azimuth_ref.qs to ref_name.qs
# remove sample metrics from cdata

data_dir <- '/srv/dseqr'
app_dirs <- list.files(data_dir)
app_dirs <- setdiff(app_dirs,
                    c('indices', 'node_modules', 'pert_signature_dir',
                      'pert_query_dir', 'gs_dir', 'example_data.tar.gz', 'tx2gene'))


convert <- function(data_dir, app_dir) {
    cat('==== Working on', app_dir, '====\n')

    sc_dir <- file.path(data_dir, app_dir, 'single-cell')

    dataset_names <- list.files(sc_dir)
    dataset_dirs <- file.path(sc_dir, dataset_names)

    for (i in seq_along(dataset_dirs)) {
        dataset_name <- dataset_names[i]
        dataset_dir <- dataset_dirs[i]
        cat('working on', dataset_names[i], '...\n')

        # remove sample columns from metrics
        shell_path <- file.path(dataset_dir, 'shell.qs')
        if (!file.exists(shell_path)) next()
        shell <- qs::qread(shell_path)

        samples <- unique(shell$batch)
        if (length(samples) > 1) {
            shell@colData <- shell@colData[, !colnames(shell@colData) %in% samples]
            qs::qsave(shell, shell_path)
        }

        # azimuth_ref to ref_name
        azimuth_path <- file.path(dataset_dir, 'azimuth_ref.qs')
        if (!file.exists(azimuth_path)) next()

        ref_name <- qs::qread(azimuth_path)
        ref_path <- file.path(dataset_dir, 'ref_name.qs')
        qs::qsave(ref_name, ref_path)
        unlink(azimuth_path)
    }
}

for  (app_dir in app_dirs) convert(data_dir, app_dir)
