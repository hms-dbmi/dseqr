context("subsetting single cell datasets works")

mock_counts <- function(...) {
    set.seed(0)
    sce <- scDblFinder::mockDoubletSCE(...)
    counts <- sce@assays@data$counts
    row.names(counts) <- paste0('gene', seq_len(nrow(counts)))
    colnames(counts) <- paste0('cell', seq_len(ncol(counts)))

    counts <- Matrix::Matrix(counts, sparse = TRUE)
    return(counts)
}

mock_mtx_files <- function(counts, features, uploaded_data_dir) {
    DropletUtils::write10xCounts(uploaded_data_dir, counts, gene.id = features$id, gene.symbol = features$name)
}

mock_uploaded_data <- function(dataset_name, sc_dir, uploaded_data_dir) {
    dir.create(sc_dir)

    counts <- mock_counts(ncells = c(200, 300, 400, 200, 500, 300), ngenes = 1000)

    features <- data.frame(id = paste0('ENSG', seq_len(nrow(counts))),
                           name = row.names(counts))

    mock_mtx_files(counts, features, uploaded_data_dir)

}


test_that("a uni-sample dataset can be subsetted", {
    # setup
    from_dataset <- 'test'
    sc_dir <- file.path(tempdir(), 'single-cell')
    uploaded_data_dir <- file.path(tempdir(), 'uploads')

    tx2gene_dir <- file.path(tempdir(), 'tx2gene')
    dir.create(tx2gene_dir)

    mock_uploaded_data(from_dataset, sc_dir, uploaded_data_dir)

    # import and keep all cells
    suppressWarnings(
        import_scseq(from_dataset, uploaded_data_dir, sc_dir, tx2gene_dir, metrics = NULL)
    )

    # remove doublets
    dataset_name <- paste0(from_dataset, '_QC1')

    expect_true(suppressWarnings(
        subset_saved_scseq(sc_dir,
                           founder = from_dataset,
                           from_dataset = from_dataset,
                           dataset_name = dataset_name,
                           subset_metrics = 'high_doublet_score')
    ))

    # check that removed doublet and nothing else
    from_scseq <- load_scseq_qs(file.path(sc_dir, from_dataset))
    expected_cells <- colnames(from_scseq[, from_scseq$doublet_class == 'singlet'])

    scseq <- load_scseq_qs(file.path(sc_dir, dataset_name))
    expect_setequal(colnames(scseq), expected_cells)

    # cleanup
    unlink(c(sc_dir, uploaded_data_dir, tx2gene_dir), recursive = TRUE)
})
