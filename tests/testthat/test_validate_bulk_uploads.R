context("validation of uploaded bulk fastq.gz files works")

test_that("files with 0 bytes are invalid", {

    # fastq_dir <- tempdir()
    # fpath <- file.path(fastq_dir, 'blah.fastq.gz')
    # file.create(fpath)

    up_df <- data.frame(size = c(0, 25, 25))
    expect_match(validate_bulk_uploads(up_df), '0 bytes')
})
