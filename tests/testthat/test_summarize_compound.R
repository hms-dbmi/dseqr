context("collapsing multiple results from same compound")
library(tibble)

query_res <- tibble('Correlation' = c(0.2, 0.3, 0.4),
                    'Compound' = rep('LY1', 3),
                    'Cell Line' = rep('MCF7', 3),
                    'Dose' = rep('1e-2M', 3),
                    'Duration' = rep('6h', 3),
                    'Samples(n)' = rep(3, 3),
                    'Clinical Phase' = factor(c('Phase 1', 'Phase 2', 'Phase 3'), ordered = TRUE),
                    'MOA' = c('dog', 'cat', NA),
                    'Target' = rep(NA, 3))

query_res <- summarize_compound(query_res)

test_that("summarize_compound collapses multiple results into 1", {
  expect_equal(nrow(query_res), 1)
})

test_that("summarize_compound puts all correlations into a list", {
  expect_equal(query_res$Correlation, I(list(c(0.2, 0.3, 0.4))))
})

test_that("summarize_compound creates title column with list", {
  expect_equal(query_res$title, I(list(rep('MCF7_1e-2M_6h_3', 3))))
})

test_that("summarize_compound pastes unique non-NA entries for other columns", {
  expect_equal(query_res$MOA, "dog | cat")
})

test_that("summarize_compound keeps most advanced clinical phase only", {
  expect_equal(as.character(query_res$`Clinical Phase`), "Phase 3")
})

test_that("summarize_compound keeps NA if all are NA", {
  expect_true(is.na(query_res$Target))
})
