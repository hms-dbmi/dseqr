context("model matrix setup and surrogate variable analysis")
library(Biobase)

test_that("diff_setup returns correct model matrix if control samples appear first", {

  # load eset and add group info
  eset_path <- system.file('extdata', 'IBD', 'eset.rds', package = 'drugseqr')
  eset <- readRDS(eset_path)
  pData(eset)$group <- rep(c('control', 'test'), each=5)

  # run setup
  setup <- diff_setup(eset, svanal=FALSE, rna_seq=TRUE)

  # controls are first 5 samples, tests and next 5 samples
  expect_equivalent(setup$mod[, 'control'], rep(c(1, 0), each=5))
  expect_equivalent(setup$mod[, 'test'], rep(c(0, 1), each=5))

})

test_that("diff_setup returns correct model matrix if test samples appear first", {

  # load eset and add group info
  eset_path <- system.file('extdata', 'IBD', 'eset.rds', package = 'drugseqr')
  eset <- readRDS(eset_path)
  pData(eset)$group <- rep(c('test', 'control'), each=5)

  # run setup
  setup <- diff_setup(eset, svanal=FALSE, rna_seq=TRUE)

  # tests are first 5 samples, control and next 5 samples
  expect_equivalent(setup$mod[, 'control'], rep(c(0, 1), each=5))
  expect_equivalent(setup$mod[, 'test'], rep(c(1, 0), each=5))
})


test_that("diff_setup removes 1:many PROBE:SYMBOL before sva", {

  # load eset and add group info
  eset_path <- system.file('extdata', 'IBD', 'eset.rds', package = 'drugseqr')
  eset <- readRDS(eset_path)
  pData(eset)$group <- rep(c('test', 'control'), each=5)

  # esets with uniq and duplicate rows
  eset_uniq <- eset[1:1000, ]
  eset_dups <- eset[rep(1:1000, each=5), ]

  setup_uniq <- diff_setup(eset_uniq, svanal = TRUE, rna_seq = TRUE)
  setup_dups <- diff_setup(eset_dups, svanal = TRUE, rna_seq = TRUE)

  # discovered surrogate variables should be equal
  expect_equal(setup_uniq$modsv, setup_dups$modsv)
})

test_that("diff_setup preserve multiple unique measures for the same gene before sva", {

  # load eset and add group info
  eset_path <- system.file('extdata', 'IBD', 'eset.rds', package = 'drugseqr')
  eset <- readRDS(eset_path)
  pData(eset)$group <- rep(c('test', 'control'), each=5)

  # esets with uniq and duplicate genes but all unique exprs
  eset <- eset[!duplicated(exprs(eset)), ]
  eset_uniq <- eset[1:1000, ]
  eset_dups <- eset_uniq
  fData(eset_dups) <- fData(eset_dups)[rep(1:100, each=10), ]

  setup_uniq <- diff_setup(eset_uniq, svanal = TRUE, rna_seq = TRUE)
  setup_dups <- diff_setup(eset_dups, svanal = TRUE, rna_seq = TRUE)

  # discovered surrogate variables should be equal
  expect_equal(setup_uniq$modsv, setup_dups$modsv)
})
