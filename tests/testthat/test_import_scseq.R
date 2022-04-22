context("importing single cell files works")

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

mock_h5_file <- function(counts, features, uploaded_data_dir, ...) {
    dir.create(uploaded_data_dir)

    fpath <- file.path(uploaded_data_dir, 'feature_barcode_matrix.h5')
    DropletUtils::write10xCounts(fpath, counts, gene.id = features$id, gene.symbol = features$name, ...)
}

mock_uploaded_data <- function(dataset_name, sc_dir, uploaded_data_dir, type = 'mtx', ...) {
    dir.create(sc_dir)

    counts <- mock_counts(ncells = c(200, 300, 400, 200, 500, 300), ngenes = 1000)

    features <- data.frame(id = paste0('ENSG', seq_len(nrow(counts))),
                           name = row.names(counts))

    if (type == 'mtx') mock_mtx_files(counts, features, uploaded_data_dir)
    else if (type == 'h5') mock_h5_file(counts, features, uploaded_data_dir, ...)

}

test_that("cellranger matrix.mtx, features.tsv, and barcodes.tsv files can be imported", {
    # setup
    dataset_name <- 'test'
    sc_dir <- file.path(tempdir(), 'single-cell')
    uploaded_data_dir <- file.path(tempdir(), 'uploads')

    tx2gene_dir <- file.path(tempdir(), 'tx2gene')
    dir.create(tx2gene_dir)

    mock_uploaded_data(dataset_name, sc_dir, uploaded_data_dir)

    expect_error(
        suppressWarnings(
            import_scseq(dataset_name, uploaded_data_dir, sc_dir, tx2gene_dir)),
        NA
    )

    # have bio
    scseq <- load_scseq_qs(file.path(sc_dir, dataset_name))
    bio <- SingleCellExperiment::rowData(scseq)$bio

    expect_true(!is.null(bio))

    # cleanup
    unlink(c(sc_dir, uploaded_data_dir, tx2gene_dir), recursive = TRUE)
})

test_that("cellranger .h5 files can be imported", {
    # setup
    dataset_name <- 'test'
    sc_dir <- file.path(tempdir(), 'single-cell')
    uploaded_data_dir <- file.path(tempdir(), 'uploads')

    tx2gene_dir <- file.path(tempdir(), 'tx2gene')
    dir.create(tx2gene_dir)

    mock_uploaded_data(dataset_name, sc_dir, uploaded_data_dir, type = 'h5')

    expect_error(
        suppressWarnings(import_scseq(dataset_name, uploaded_data_dir, sc_dir, tx2gene_dir)),
        NA
    )

    # cleanup
    unlink(c(sc_dir, uploaded_data_dir, tx2gene_dir), recursive = TRUE)
})

test_that("cellranger .h5 files with multiple assays can be imported", {
    # setup
    dataset_name <- 'test'
    sc_dir <- file.path(tempdir(), 'single-cell')
    uploaded_data_dir <- file.path(tempdir(), 'uploads')

    tx2gene_dir <- file.path(tempdir(), 'tx2gene')
    dir.create(tx2gene_dir)

    mock_uploaded_data(dataset_name,
                       sc_dir,
                       uploaded_data_dir,
                       type = 'h5',
                       version = '3',
                       gene.type = c(rep('Gene Expression', 199), 'Antibody Capture'))

    expect_error(
        suppressWarnings(import_scseq(dataset_name, uploaded_data_dir, sc_dir, tx2gene_dir)),
        NA
    )

    # cleanup
    unlink(c(sc_dir, uploaded_data_dir, tx2gene_dir), recursive = TRUE)
})
