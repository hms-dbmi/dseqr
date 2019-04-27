context("sanity check diff_expr")
library(Biobase)
library(limma)

test_that("diff_expr run without sva", {

  # load ExpressionSet
  data_dir <- system.file('extdata', 'IBD', package='drugseqr')
  eset <- readRDS(file.path(data_dir, 'eset.rds'))
  group <- factor(pData(eset)$`Binary Affected`, level = c('Normal', 'Abnormal'))

  # run diff_expr without sva
  # mock previous analysis (to skip UI)
  prev_anal <- list(pdata = data.frame(group = c('control', 'test')[group], row.names = sampleNames(eset)))
  expect_error(diff_expr(eset, data_dir, svanal = FALSE, prev_anal = prev_anal), NA)

  # cleanup
  unlink(system.file('extdata', 'IBD', 'diff_expr_symbol.rds', package='drugseqr'))
  unlink(system.file('tests', 'testthat', 'Rplots.pdf', package='drugseqr'))

})

test_that("diff_expr run with sva", {

  # load ExpressionSet
  data_dir <- system.file('extdata', 'IBD', package='drugseqr')
  eset <- readRDS(file.path(data_dir, 'eset.rds'))
  group <- factor(pData(eset)$`Binary Affected`, level = c('Normal', 'Abnormal'))

  # run diff_expr without sva
  # mock previous analysis (to skip UI)
  prev_anal <- list(pdata = data.frame(group = c('control', 'test')[group], row.names = sampleNames(eset)))
  expect_error(diff_expr(eset, data_dir, svanal = TRUE, prev_anal = prev_anal), NA)

  # cleanup
  unlink(system.file('extdata', 'IBD', 'diff_expr_symbol.rds', package='drugseqr'))
  unlink(system.file('tests', 'testthat', 'Rplots.pdf', package='drugseqr'))
})

test_that("diff_expr produces similar results to tximport vignette limma-voom section", {

  # tximport vignette: see https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#limma-voom

  # load DGElist (saved after calcNormFactors line)
  data_dir <- system.file('extdata', 'IBD', package='drugseqr')
  dgel <- readRDS(file.path(data_dir, 'dgel.rds'))

  # load ExpressionSet
  eset <- readRDS(file.path(data_dir, 'eset.rds'))

  # make sure samples are in same order
  expect_equal(row.names(pData(eset)), row.names(dgel$samples))

  # run voom and standard limma pipeline using dgel
  group <- factor(pData(eset)$`Binary Affected`, level = c('Normal', 'Abnormal'))
  design <- model.matrix(~group)
  v <- voom(dgel, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  top_table <- topTable(fit, coef=2, Inf)

  # run diff_expr without sva
  # mock previous analysis (to skip UI)
  prev_anal <- list(pdata = data.frame(group = c('control', 'test')[group], row.names = sampleNames(eset)))
  anal <- diff_expr(eset, data_dir, svanal = FALSE, prev_anal = prev_anal)

  expect_equal(anal$top_table, top_table)

  # run diff_expr with sva
  # put top_table in same order
  sva_anal <- diff_expr(eset, data_dir, prev_anal = prev_anal)
  sva_anal$top_table <- sva_anal$top_table[row.names(top_table), ]

  logFC_cor <- cor(sva_anal$top_table$logFC, top_table$logFC)

  cat('logFC correlation (sva vs not):', logFC_cor)
  expect_equal(logFC_cor, 0.6215673, tolerance=1e-3)

  # cleanup
  unlink(system.file('extdata', 'IBD', 'diff_expr_symbol.rds', package='drugseqr'))
  unlink(system.file('tests', 'testthat', 'Rplots.pdf', package='drugseqr'))
})
