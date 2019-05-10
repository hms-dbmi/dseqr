context("select pairs validation")

test_that("validate_reps requires atleast two rows", {

  # single row fails
  rows <- 1
  pairs <- rep(NA, 4)
  reps <- rep(NA, 4)
  expect_message(validate_reps(pairs, rows, reps), 'at least two rows')
  expect_false(validate_reps(pairs, rows, reps))

  # using two rows works
  rows2 <- c(1, 2)
  expect_true(validate_reps(pairs, rows2, reps))
})

test_that("validate_reps won't overwrite existing replicates", {

  rows <- c(1, 2, 3, 4)
  pairs <- rep(NA, 6)
  reps <- c(1, 1, NA, NA, 2, 2)

  expect_message(validate_reps(pairs, rows, reps), 'already belong to a replicate')
  expect_false(validate_reps(pairs, rows, reps))

  rows2 <- c(3, 4)
  expect_true(validate_reps(pairs, rows2, reps))
})

test_that("validate_reps won't allow paired rows to be marked as replicates", {

  rows <- c(1, 2)
  pairs <- c(1, 1, NA, NA, NA)
  reps <- rep(NA, 5)

  expect_message(validate_reps(pairs, rows, reps), 'Replicates must be specified first')
  expect_false(validate_reps(pairs, rows, reps))

  rows2 <- c(3, 4)
  expect_true(validate_reps(pairs, rows2, reps))
})
