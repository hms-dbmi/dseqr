context("drugseqr app")
# This file is for testing the applications in the inst/ directory.

library(shinytest)

test_that("app works", {
  appdir <- system.file(package = "drugseqr", "app")
  expect_pass(testApp(appdir, compareImages = FALSE))
})
