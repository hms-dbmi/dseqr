mock_counts <- function() {
    set.seed(0)
    sce <- scater::mockSCE()
    counts <- sce@assays@data$counts
    counts <- Matrix::Matrix(counts, sparse = TRUE)
    return(counts)
}

mock_h5_file <- function(counts, features, uploaded_data_dir, genome) {
    dir.create(uploaded_data_dir, showWarnings = FALSE)

    fpath <- file.path(uploaded_data_dir, 'feature_barcode_matrix.h5')
    DropletUtils::write10xCounts(fpath, counts, genome=genome, gene.id = features$id, gene.symbol = features$name)
}

mock_uploaded_data <- function(sc_dir, uploaded_data_dir, genome) {
    dir.create(sc_dir, showWarnings = FALSE)
    counts <- mock_counts()

    features <- data.frame(id = paste0('ENSG', seq_len(nrow(counts))),
                           name = row.names(counts))

    mock_h5_file(counts, features, uploaded_data_dir, genome)

}

test_that("a sample with one matrix.mtx, features.tsv, and barcodes.tsv is valid", {
    up_df <- data.frame(name = c('matrix.mtx', 'features.tsv', 'barcodes.tsv'))
    samples <- rep('a', 3)

    expect_null(validate_scseq_import(up_df, samples))
})

test_that("a sample missing one matrix.mtx, features.tsv, or barcodes.tsv is invalid", {
    up_df <- data.frame(name = c('matrix.mtx', 'features.tsv', 'barcodes.tsv'), datapath = NA)
    samples <- rep('a', 3)

    expect_type(validate_scseq_import(up_df[-1, ], samples), 'character')
    expect_type(validate_scseq_import(up_df[-2, ], samples), 'character')
    expect_type(validate_scseq_import(up_df[-3, ], samples), 'character')
})

test_that("a sample with more than 3 files including an .mtx file is invalid", {

    # different type of file
    up_df <- data.frame(name = c('matrix.mtx', 'features.tsv', 'barcodes.tsv', 'blah.txt'), datapath = NA)
    samples <- rep('a', 4)

    expect_type(validate_scseq_import(up_df, samples), 'character')

    # two samples labeled as one
    names <- paste0(c('GSM1_', 'GSM2_'), rep(c('matrix.mtx', 'features.tsv', 'barcodes.tsv'), each = 2))
    up_df <- data.frame(name = names, datapath = NA)
    expect_type(validate_scseq_import(up_df, samples), 'character')
})


test_that("a sample missing a name is invalid", {
    up_df <- data.frame(name = c('matrix.mtx', 'features.tsv', 'barcodes.tsv'), datapath = NA)
    samples <- rep(NA, 3)

    expect_type(validate_scseq_import(up_df, samples), 'character')
})


test_that("an h5 file with matrix slot is valid", {
    sc_dir <- file.path(tempdir(), 'single-cell')
    uploaded_data_dir <- file.path(tempdir(), 'uploads')
    mock_uploaded_data(sc_dir, uploaded_data_dir, 'matrix')

    up_df <- data.frame(
        name = c('blah.h5'),
        datapath = file.path(uploaded_data_dir, 'feature_barcode_matrix.h5'))

    samples <- 'a'
    expect_null(validate_scseq_import(up_df, samples))

    # cleanup
    unlink(c(sc_dir, uploaded_data_dir), recursive = TRUE)
})


test_that("an h5 file without matrix slot is valid", {
    sc_dir <- file.path(tempdir(), 'single-cell')
    uploaded_data_dir <- file.path(tempdir(), 'uploads')
    mock_uploaded_data(sc_dir, uploaded_data_dir, 'other')

    up_df <- data.frame(
        name = c('blah.h5'),
        datapath = file.path(uploaded_data_dir, 'feature_barcode_matrix.h5'))

    samples <- 'a'
    expect_null(validate_scseq_import(up_df, samples))

    # cleanup
    unlink(c(sc_dir, uploaded_data_dir), recursive = TRUE)
})


test_that("multiple samples in MEX format are valid", {
    up_df <- data.frame(name = rep(c('matrix.mtx', 'features.tsv', 'barcodes.tsv'), 2))
    samples <- rep(c('a', 'b'), each = 3)

    expect_null(validate_scseq_import(up_df, samples))
})


test_that("multiple samples in h5 format are valid", {

    # create one H5 file
    sc_dir <- file.path(tempdir(), 'single-cell')
    uploaded_data_dir <- file.path(tempdir(), 'uploads')
    mock_uploaded_data(sc_dir, uploaded_data_dir, 'matrix')

    # move it then create another
    fpath1 <- file.path(uploaded_data_dir, 'feature_barcode_matrix.h5')
    fpath2 <- file.path(uploaded_data_dir, 'other_feature_barcode_matrix.h5')
    file.rename(fpath1, fpath2)
    mock_uploaded_data(sc_dir, uploaded_data_dir, 'matrix')

    up_df <- data.frame(
        name = c('a.h5', 'b.h5'),
        datapath = c(fpath1, fpath2))

    samples <- c('a', 'b')
    expect_null(validate_scseq_import(up_df, samples))

    # cleanup
    unlink(c(sc_dir, uploaded_data_dir), recursive = TRUE)
})


test_that("one sample in MEX format and one in h5 is valid", {

    # create one H5 file
    sc_dir <- file.path(tempdir(), 'single-cell')
    uploaded_data_dir <- file.path(tempdir(), 'uploads')
    mock_uploaded_data(sc_dir, uploaded_data_dir, 'matrix')

    up_df <- data.frame(
        name = c('a.h5', 'matrix.mtx', 'features.tsv', 'barcodes.tsv'),
        datapath = c(file.path(uploaded_data_dir, 'feature_barcode_matrix.h5'), rep(NA, 3)))

    samples <- c('a', rep('b', 3))
    expect_null(validate_scseq_import(up_df, samples))

    # cleanup
    unlink(c(sc_dir, uploaded_data_dir), recursive = TRUE)
})

test_that("sc samples cannot have invalid names", {

  invalid_start_period <- '.sample_name'
  msg_name <- validate_scseq_add_sample(rows = 1:3, invalid_start_period)
  expect_true(grepl('not start with a period', msg_name))

  invalid_has_path <- '/etc'
  msg_name <- validate_scseq_add_sample(rows = 1:3, invalid_has_path)
  expect_true(grepl('not form a path', msg_name))

  valid_name <- 'sample_name'
  msg_name <- validate_scseq_add_sample(rows= 1:3, valid_name)
  expect_null(msg_name)
})

