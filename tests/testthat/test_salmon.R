context("test that salmon will run")
library(jsonlite)

test_that("run_salmon works", {

  # setup
  ibd_dir <- system.file('extdata', 'IBD', package="drugseqr")
  test_dir <- file.path(ibd_dir, 'test')
  dir.create(test_dir)
  file.copy(file.path(ibd_dir, 'IBD_1000 - 1285-01-G1_S9_R1_001.fastq.gz'), test_dir)

  # mock pdata to skip GUI
  pdata <- data.frame('File Name' = 'IBD_1000 - 1285-01-G1_S9_R1_001.fastq.gz', stringsAsFactors = FALSE, check.names = FALSE)

  pdata <- run_salmon(test_dir, pdata = pdata)

  # should create quants dir
  qdir <- system.file('extdata', 'IBD', 'test', 'quants', package="drugseqr")
  expect_true(nchar(qdir) > 0 && dir.exists(qdir))

  # should contain 1 directory inside with sample name
  expect_true(dir.exists(file.path(qdir, 'IBD_1000 - 1285-01-G1_S9_R1_001')))

  # teardown
  unlink(test_dir, recursive = TRUE)
})

test_that("run_salmon handles replicate files", {

  # setup
  fnames <- c('IBD_1000 - 1285-01-G7-1_S16_R1_001.fastq.gz', 'IBD_1000 - 1285-01-G7-2_S17_R1_001.fastq.gz')
  ibd_dir <- system.file('extdata', 'IBD', package="drugseqr")
  test_dir <- file.path(ibd_dir, 'test')
  dir.create(test_dir)
  source_paths <- file.path(ibd_dir, fnames)

  file.copy(source_paths, test_dir)

  # mock pdata to skip GUI
  pdata <- data.frame('File Name' = fnames,
                      'Replicate' = c(1, 1),
                      stringsAsFactors = FALSE, check.names = FALSE)

  run_salmon(test_dir, pdata = pdata)

  # should create quants dir
  qdir <- system.file('extdata', 'IBD', 'test', 'quants', package="drugseqr")
  expect_true(nchar(qdir) > 0 && dir.exists(qdir))

  # should contain 1 directory inside with first sample as name
  sdir <- file.path(qdir, 'IBD_1000 - 1285-01-G7-1_S16_R1_001')
  expect_true(dir.exists(sdir))

  # should have used both files for quants
  cmd_info <- read_json(file.path(sdir, 'cmd_info.json'))
  dest_paths <- file.path(test_dir, fnames)

  expect_length(cmd_info$unmatedReads, 2)
  expect_setequal(unlist(cmd_info$unmatedReads), dest_paths)

  # teardown
  unlink(test_dir, recursive = TRUE)
})
