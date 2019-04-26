context("model matrix setup and surrogate variable analysis")
library(Biobase)

test_that("model matrix is correct if control samples appear first", {

  # load eset and add group info
  eset_path <- system.file('extdata', 'IBD', 'eset.rds', package = 'drugseqr')
  eset <- readRDS(eset_path)
  pData(eset)$group <- rep(c('control', 'test'), each=6)

  # run setup
  setup <- diff_setup(eset, svanal=FALSE, rna_seq=TRUE)

  # controls are first 6 samples, tests and next 6 samples
  expect_equivalent(setup$mod[, 'control'], rep(c(1, 0), each=6))
  expect_equivalent(setup$mod[, 'test'], rep(c(0, 1), each=6))

})

test_that("model matrix is correct if test samples appear first", {

  # load eset and add group info
  eset_path <- system.file('extdata', 'IBD', 'eset.rds', package = 'drugseqr')
  eset <- readRDS(eset_path)
  pData(eset)$group <- rep(c('test', 'control'), each=6)

  # run setup
  setup <- diff_setup(eset, svanal=FALSE, rna_seq=TRUE)

  # tests are first 6 samples, control and next 6 samples
  expect_equivalent(setup$mod[, 'control'], rep(c(0, 1), each=6))
  expect_equivalent(setup$mod[, 'test'], rep(c(1, 0), each=6))
})
