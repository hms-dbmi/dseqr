test_that("files with 0 bytes are invalid", {

    up_df <- data.frame(size = c(0, 25, 25))
    expect_match(validate_bulk_uploads(up_df), '0 bytes')
})

test_that("fastq files with older illumina sequence identifiers are valid", {

    fastq_dir <- tempdir()
    fq1_path <- file.path(fastq_dir, 'A1.fastq.gz')
    fq2_path <- file.path(fastq_dir, 'A2.fastq.gz')
    writeLines('@HWUSI-EAS100R:6:73:941:1973#0/1', fq1_path)
    writeLines('@HWUSI-EAS100R:6:73:941:1973#0/1', fq2_path)


    up_df <- data.frame(size = c(1, 1), datapath = c(fq1_path, fq2_path))
    expect_null(validate_bulk_uploads(up_df))

    unlink(c(fq1_path, fq2_path))
})


test_that("fastq files with newer illumina sequence identifiers are valid", {

    fastq_dir <- tempdir()
    fq1_path <- file.path(fastq_dir, 'A1.fastq.gz')
    fq2_path <- file.path(fastq_dir, 'A2.fastq.gz')
    writeLines('@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG', fq1_path)
    writeLines('@EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG', fq2_path)


    up_df <- data.frame(size = c(1, 1), datapath = c(fq1_path, fq2_path))
    expect_null(validate_bulk_uploads(up_df))

    unlink(c(fq1_path, fq2_path))
})


test_that("a fastq file with pair identifier '2' needs a fastq file with pair identifier '1'", {

    fastq_dir <- tempdir()
    fq2_path <- file.path(fastq_dir, 'A2.fastq.gz')
    writeLines('@EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG', fq2_path)


    up_df <- data.frame(size = 1, datapath = fq2_path)
    expect_match(validate_bulk_uploads(up_df), 'both files for pair')

    unlink(fq2_path)
})

test_that("there must be an equal number of fastq files with pair identifiers '1' and '2'", {

    fastq_dir <- tempdir()
    fq1_path <- file.path(fastq_dir, 'A1.fastq.gz')
    fq2_path <- file.path(fastq_dir, 'A2.fastq.gz')
    fq3_path <- file.path(fastq_dir, 'B1.fastq.gz')
    fq4_path <- file.path(fastq_dir, 'B2.fastq.gz')

    writeLines('@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG', fq1_path)
    writeLines('@EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG', fq2_path)
    writeLines('@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG', fq3_path)
    writeLines('@EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG', fq4_path)


    up_df <- data.frame(size = c(1, 1, 1, 1), datapath = c(fq1_path, fq2_path, fq3_path, fq4_path))

    # passes with all
    expect_null(validate_bulk_uploads(up_df))

    # fails if remove any
    expect_match(validate_bulk_uploads(up_df[-1, ]), 'both files for pair')
    expect_match(validate_bulk_uploads(up_df[-2, ]), 'both files for pair')
    expect_match(validate_bulk_uploads(up_df[-3, ]), 'both files for pair')
    expect_match(validate_bulk_uploads(up_df[-4, ]), 'both files for pair')

    unlink(c(fq1_path, fq2_path, fq3_path, fq4_path))
})


test_that("fastq files cannot have the same name", {

    fastq_dir <- tempdir()
    fq1_path <- file.path(fastq_dir, 'A1.fastq.gz')
    fq2_path <- file.path(fastq_dir, 'A2.fastq.gz')
    fq3_path <- file.path(fastq_dir, 'B1.fastq.gz')
    fq4_path <- file.path(fastq_dir, 'B2.fastq.gz')

    writeLines('@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG', fq1_path)
    writeLines('@EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG', fq2_path)
    writeLines('@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG', fq3_path)
    writeLines('@EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG', fq4_path)


    up_df <- data.frame(size = c(1, 1, 1, 1),
                        name = c('fq1', 'f2', 'fq3', 'fq4'),
                        datapath = c(fq1_path, fq2_path, fq3_path, fq4_path))

    # passes when all names different
    expect_null(validate_bulk_uploads(up_df))

    # fails if two are the same
    up_df$name[1:2] <- 'fq1'
    expect_match(validate_bulk_uploads(up_df), 'unique names')

    unlink(c(fq1_path, fq2_path, fq3_path, fq4_path))
})

test_that("bulk datasets cannot have invalid names", {

  invalid_start_period <- '.dataset_name'
  msg_name <- validate_bulk_name(invalid_start_period)
  expect_true(grepl('not start with a period', msg_name))

  invalid_path <- '/etc'
  msg_name <- validate_bulk_name(invalid_path)
  expect_true(grepl('not form a path', msg_name))

  valid_name <- 'dataset_name'
  msg_name <- validate_bulk_name(valid_name)
  expect_null(msg_name)
})
