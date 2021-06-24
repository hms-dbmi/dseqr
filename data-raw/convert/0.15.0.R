# convert azimuth ref to more portable format
library(SingleCellExperiment)

app_dir <- "/home/alex/patient_data/example"
app_dirs <- list.files('/srv/dseqr')
app_dirs <- setdiff(app_dirs,
                    c('indices', 'node_modules', 'pert_signature_dir',
                      'pert_query_dir', 'example_data.tar.gz'))

is.numstring <- function(x) !is.na(suppressWarnings(as.numeric(x)))

convert <- function(app_dir) {
    cat('==== Working on', app_dir, '====\n')

    sc_dir <- file.path(app_dir, 'single-cell')

    dataset_names <- list.files(sc_dir)
    dataset_dirs <- file.path(sc_dir, dataset_names)
    azi_paths <- file.path(dataset_dirs, 'azimuth_ref.qs')
    is.azi <- file.exists(azi_paths)
    azi_paths <- azi_paths[is.azi]
    dataset_dirs <- dataset_dirs[is.azi]
    dataset_names <- dataset_names[is.azi]

    for (i in seq_along(dataset_dirs)) {
        dataset_name <- dataset_names[i]
        dataset_dir <- dataset_dirs[i]
        cat('working on', dataset_names[i], '...\n')

        scseq_path <- file.path(dataset_dir, 'shell.qs')
        if (!file.exists(scseq_path)) next()

        # skip if already in correct format
        resoln_path <- file.path(dataset_dir, 'resoln.qs')
        resoln <- qs::qread(resoln_path)
        if (!is.numstring(resoln)) next()

        scseq <- qs::qread(scseq_path)

        # revert names of azimuth slots
        scseq$predicted.celltype.l1 <- scseq$cluster_l1
        scseq$predicted.celltype.l2 <- scseq$cluster_l2
        scseq$predicted.celltype.l3 <- scseq$cluster_l3
        scseq$predicted.celltype.l1.score <- scseq$cluster_l1_score
        scseq$predicted.celltype.l2.score <- scseq$cluster_l2_score
        scseq$predicted.celltype.l3.score <- scseq$cluster_l3_score

        scseq$cluster_l1 <- scseq$cluster_l2 <- scseq$cluster_l3 <- NULL
        scseq$cluster_l1_score <- scseq$cluster_l2_score <- scseq$cluster_l3_score <- NULL
        scseq$mapping.score <- scseq$mapping_score
        scseq$mapping_score <- NULL

        qs::qsave(scseq, scseq_path, preset = 'fast')

        # change resoln
        resoln <- switch(resoln,
                         '1' = 'predicted.celltype.l1',
                         '2' = 'predicted.celltype.l2',
                         '3' = 'predicted.celltype.l3')

        qs::qsave(resoln, resoln_path)

        # change directory names
        file.rename(file.path(dataset_dir, 'snn1'), file.path(dataset_dir, 'predicted.celltype.l1'))
        file.rename(file.path(dataset_dir, 'snn2'), file.path(dataset_dir, 'predicted.celltype.l2'))
        file.rename(file.path(dataset_dir, 'snn3'), file.path(dataset_dir, 'predicted.celltype.l3'))
    }
}

for  (app_dir in app_dirs) convert(app_dir)
