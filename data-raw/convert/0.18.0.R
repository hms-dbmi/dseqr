# add species.qs for individual datasets
library(SingleCellExperiment)

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

        species_path <- file.path(dataset_dir, 'species.qs')
        if (file.exists(species_path)) next()

        shell_path <- file.path(dataset_dir, 'shell.qs')
        if (!file.exists(shell_path)) next()

        scseq <- qs::qread(shell_path)

        species <- scseq@metadata$species
        qs::qsave(species, species_path)
    }
}

for  (app_dir in app_dirs) convert(app_dir)
