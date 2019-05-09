context("Calculating x position of correlation value")

test_that("calcx return midpoint of width when correlation is 0", {
  xpos <- calcx(0, width=200)
  expect_equal(xpos, 100)
})

test_that("calcx without padding returns exactly devided values", {
  expect_equal(calcx(-1, width=200, pad = 0), 0)
  expect_equal(calcx(-0.5, width=200, pad = 0), 50)
  expect_equal(calcx(0, width=200, pad = 0), 100)
  expect_equal(calcx(0.5, width=200, pad = 0), 150)
  expect_equal(calcx(1, width=200, pad = 0), 200)
})

test_that("calcx with padding avoids endpoints", {
  expect_gt(calcx(-1), calcx(-1, pad = 0))
  expect_lt(calcx(1), calcx(1, pad = 0))
})
