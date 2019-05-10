context("select pairs validation")

test_that("validate_pairs requires two rows to be paired", {

  # single row fails
  rows <- 1
  pairs <- rep(NA, 4)
  reps <- rep(NA, 4)
  expect_message(validate_pairs(pairs, rows, reps), 'at least two rows')
  expect_false(validate_pairs(pairs, rows, reps))

  # using two rows works
  rows2 <- c(1, 2)
  expect_true(validate_pairs(pairs, rows2, reps))
})

test_that("validate_pairs requires two non-replicate rows to be paired", {
  rows <- c(1, 2)
  pairs <- rep(NA, 4)
  reps <- c(1, 1, 2, 2)

  expect_message(validate_pairs(pairs, rows, reps), 'at least two non-replicate rows')
  expect_false(validate_pairs(pairs, rows, reps))

  rows2 <- c(1, 2, 3, 4)
  expect_true(validate_pairs(pairs, rows2, reps))
})

test_that("validate_pairs doesn't allow a replicate not included in a pair", {
  rows <- c(1, 3)
  pairs <- rep(NA, 4)
  reps <- c(1, 1, 2, 2)

  expect_message(validate_pairs(pairs, rows, reps), 'All replicates must be included in the same pair')
  expect_false(validate_pairs(pairs, rows, reps))

  rows2 <- c(1, 2, 3, 4)
  expect_true(validate_pairs(pairs, rows2, reps))
})

test_that("validate_pairs requires exactly two unique replicates, or one replicate and 1 additional sample", {
  rows <- c(1, 2, 3, 4, 5)
  pairs <- rep(NA, 6)
  reps <- c(1, 1, 2, 2, NA, NA)

  expect_message(validate_pairs(pairs, rows, reps), 'exactly two replicate groups')
  expect_false(validate_pairs(pairs, rows, reps))

  rows2 <- c(1, 2, 3, 4)
  expect_true(validate_pairs(pairs, rows2, reps))
})

test_that("validate_pairs doesn't allow overwriting of existing pairs", {
  rows <- c(1, 2, 3, 4)
  reps <- rep(NA, 6)
  pairs <- c(1, 1, NA, NA, 2, 2)

  expect_message(validate_pairs(pairs, rows, reps), 'already belong to a pair')
  expect_false(validate_pairs(pairs, rows, reps))

  rows2 <- c(3, 4)
  expect_true(validate_pairs(pairs, rows2, reps))
})


