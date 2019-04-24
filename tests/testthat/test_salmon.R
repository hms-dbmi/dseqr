context("test that salmon will run")

test_that("run_salmon works", {

  # setup
  ibd_dir <- system.file('extdata', 'IBD', package="drugseqr")
  test_dir <- file.path(ibd_dir, 'test')
  dir.create(test_dir)
  file.copy(file.path(ibd_dir, 'IBD_1000 - 1285-01-G1_S9_R1_001.fastq.gz'), test_dir)

  run_salmon(test_dir)

  # should create quants dir
  qdir <- system.file('extdata', 'IBD', 'test', 'quants', package="drugseqr")
  expect_true(nchar(qdir) > 0 && dir.exists(qdir))

  # should contain 1 directory inside with sample name
  expect_true(dir.exists(file.path(qdir, 'IBD_1000 - 1285-01-G1_S9_R1_001')))

  # teardown
  unlink(test_dir, recursive = TRUE)
})
