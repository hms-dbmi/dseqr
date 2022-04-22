context("subsetting single cell datasets works")

mock_scseq_files <- function(sc_dir, dataset_name, sample_names = 'a') {
    dataset_dir <- file.path(sc_dir, dataset_name)
    dir.create(dataset_dir, recursive = TRUE)

    sce <- scuttle::mockSCE()
    sce <- scuttle::logNormCounts(sce)

    # add clusters and batch (sample)
    sce$cluster <- factor(scran::quickCluster(sce, min.size=50))
    sce$batch <- sample(sample_names, size = ncol(sce), replace = TRUE)

    # add some doublets
    set.seed(0)
    doublet_idx <- sample(ncol(sce), size = ncol(sce)*0.1)
    sce$high_doublet_score <- FALSE
    sce$high_doublet_score[doublet_idx] <- TRUE

    # save clusters
    anal <- list(clusters = sce$cluster)
    dataset_subname <- file.path(dataset_name, 'snn1')
    save_scseq_data(anal, dataset_subname, sc_dir, overwrite = FALSE)

    # save scseq
    suppressWarnings(split_save_scseq(sce, dataset_dir))

    return(sce)
}

mock_integration_files <- function(dataset_name, sc_dir, integration_type = 'harmony') {

    # need args with previous integration type
    args <- list(integration_type = integration_type)
    save_scseq_args(args, dataset_name, sc_dir)

    # need integrated.qs
    qs::qsave(dataset_name, file.path(sc_dir, 'integrated.qs'))
}


test_that("a uni-sample dataset can be subsetted", {
    # setup
    from_dataset <- 'test'
    sc_dir <- file.path(tempdir(), 'single-cell')
    scseq <- mock_scseq_files(sc_dir, from_dataset)

    # remove doublets
    dataset_name <- paste0(from_dataset, '_QC1')

    expect_true(suppressWarnings(
        subset_saved_scseq(sc_dir,
                           founder = from_dataset,
                           from_dataset = from_dataset,
                           dataset_name = dataset_name,
                           subset_metrics = 'high_doublet_score')
    ))

    # check that removed doublets and nothing else
    expected_cells <- colnames(scseq)[!scseq$high_doublet_score]

    subset_scseq <- load_scseq_qs(file.path(sc_dir, dataset_name))
    expect_setequal(colnames(subset_scseq), expected_cells)

    # cleanup
    unlink(sc_dir, recursive = TRUE)
})



test_that("an integrated dataset can be subsetted", {
    # setup
    from_dataset <- 'test'
    sc_dir <- file.path(tempdir(), 'single-cell')
    scseq <- mock_scseq_files(sc_dir, from_dataset, sample_names = c('a', 'b'))

    integration_type <- 'harmony'
    mock_integration_files(from_dataset, sc_dir, integration_type)

    # remove doublets
    dataset_name <- paste0(from_dataset, '_QC1')

    expect_true(suppressWarnings(
        subset_saved_scseq(sc_dir,
                           founder = from_dataset,
                           from_dataset = from_dataset,
                           dataset_name = dataset_name,
                           subset_metrics = 'high_doublet_score',
                           is_integrated = TRUE)
    ))

    # check that removed doublets and nothing else
    expected_cells <- colnames(scseq)[!scseq$high_doublet_score]

    subset_scseq <- load_scseq_qs(file.path(sc_dir, paste0(dataset_name, '_harmony')))
    expect_setequal(colnames(subset_scseq), expected_cells)

    # cleanup
    unlink(sc_dir, recursive = TRUE)
})
