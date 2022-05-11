test_that("three unpaired samples is valid", {

    up_df <- data.frame(name = c('A.fastq.gz', 'B.fastq.gz', 'C.fastq.gz'))
    pairs <- reps <- rep(NA, nrow(up_df))
    expect_null(validate_bulk_labels(up_df, reps, pairs, paired = FALSE))
})

test_that("two unpaired samples is not valid", {

    up_df <- data.frame(name = c('A.fastq.gz', 'B.fastq.gz'))
    pairs <- reps <- rep(NA, nrow(up_df))
    expect_match(validate_bulk_labels(up_df, reps, pairs, paired = FALSE), 'At least 3 samples')
})

test_that("three paired samples is valid", {

    up_df <- data.frame(name = c('A1.fastq.gz', 'A2.fastq.gz', 'B1.fastq.gz', 'B2.fastq.gz', 'C1.fastq.gz', 'C2.fastq.gz'))
    reps <- rep(NA, nrow(up_df))
    pairs <- c(1, 1, 2, 2, 3, 3)
    expect_null(validate_bulk_labels(up_df, reps, pairs, paired = TRUE))
})

test_that("two paired samples is not valid", {

    up_df <- data.frame(name = c('A1.fastq.gz', 'A2.fastq.gz', 'B1.fastq.gz', 'B2.fastq.gz'))
    reps <- rep(NA, nrow(up_df))
    pairs <- c(1, 1, 2, 2)
    expect_match(validate_bulk_labels(up_df, reps, pairs, paired = TRUE), 'At least 3 samples')
})

test_that("Paired files must all have an associated pair", {

    up_df <- data.frame(name = c('A1.fastq.gz', 'A2.fastq.gz', 'B1.fastq.gz', 'B2.fastq.gz', 'C1.fastq.gz', 'C2.fastq.gz'))
    reps <- rep(NA, nrow(up_df))
    pairs <- c(1, 1, 2, 2, 3, 3)
    expect_null(validate_bulk_labels(up_df, reps, pairs, paired = TRUE))

    # remove a pair
    pairs[1:2] <- NA
    expect_match(validate_bulk_labels(up_df, reps, pairs, paired = TRUE), 'must belong to a pair')
})
