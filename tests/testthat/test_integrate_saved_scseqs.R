mock_scseq_files <- function(sc_dir, dataset_name, sample_names = 'a', sce = NULL) {
    dataset_dir <- file.path(sc_dir, dataset_name)
    dir.create(dataset_dir, recursive = TRUE)

    if (is.null(sce)) sce <- scuttle::mockSCE()
    sce <- scuttle::logNormCounts(sce)

    # add clusters and batch (sample)
    sce$cluster <- factor(scran::quickCluster(sce, min.size=50))
    sce$batch <- sample(sample_names, size = ncol(sce), replace = TRUE)
    sce@metadata$species <- 'Homo sapiens'

    # save clusters
    anal <- list(clusters = sce$cluster)
    dataset_subname <- file.path(dataset_name, 'snn1')
    save_scseq_data(anal, dataset_subname, sc_dir, overwrite = FALSE)

    # save scseq
    suppressWarnings(split_save_scseq(sce, dataset_dir))

    return(sce)
}

test_that("multiple scseq datasets can be integrated with harmony", {
    # setup
    from_datasets <- c('test1', 'test2')
    sc_dir <- file.path(tempdir(), 'single-cell')
    dataset_dirs <- file.path(sc_dir, from_datasets)

    # two saved datasets to integrate
    scseq <- mock_scseq_files(sc_dir, from_datasets[1])
    fs::dir_copy(dataset_dirs[1], dataset_dirs[2])

    integration_name <- 'test1_vs_test2'
    integration_type <- 'harmony'

    expect_true(suppressWarnings(
        integrate_saved_scseqs(sc_dir,
                               integration_name = integration_name,
                               dataset_names = from_datasets,
                               integration_type = integration_type)
    ))

    # check that twice as many cells in new dataset
    integrated_dir <- paste0(integration_name, '_', integration_type)
    integrated_scseq <- load_scseq_qs(file.path(sc_dir, integrated_dir))

    expect_equal(ncol(integrated_scseq), ncol(scseq)*2)

    # cleanup
    unlink(sc_dir, recursive = TRUE)
})


test_that("multiple scseq datasets can be integrated with Azimuth", {
    # setup
    counts <- utils::read.table(
        file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
        as.is = TRUE
    )

    counts <- Matrix::Matrix(as.matrix(counts), sparse = TRUE)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))

    from_datasets <- c('test1', 'test2')
    sc_dir <- file.path(tempdir(), 'single-cell')
    dataset_dirs <- file.path(sc_dir, from_datasets)

    # two saved datasets to integrate
    scseq <- suppressWarnings(
        mock_scseq_files(sc_dir, from_datasets[1], sce = sce)
    )
    fs::dir_copy(dataset_dirs[1], dataset_dirs[2])

    # run integration
    integration_name <- 'test1_vs_test2'
    integration_type <- 'Azimuth'

    expect_true(suppressWarnings(
        integrate_saved_scseqs(sc_dir,
                               integration_name = integration_name,
                               dataset_names = from_datasets,
                               integration_type = integration_type,
                               ref_name = 'human_pbmc')
    ))

    # check that twice as many cells in new dataset
    integrated_dir <- paste0(integration_name, '_', integration_type)
    integrated_scseq <- load_scseq_qs(file.path(sc_dir, integrated_dir))

    expect_equal(ncol(integrated_scseq), ncol(scseq)*2)

    # cleanup
    unlink(sc_dir, recursive = TRUE)
})


