test_that("dMcast dimensions are correct", {
  tmp <- data.frame(lapply(letters, function(x) (rep(letters, 1500))))
  colnames(tmp) <- letters
  d <- dMcast(tmp, ~.)
  expect_equal(nrow(d), nrow(tmp))
  a <- dMcast(warpbreaks, ~.)
  expect_equal(nrow(a), nrow(warpbreaks))
  b <- dMcast(warpbreaks, ~., as.factors = TRUE)
  expect_lt(ncol(a), ncol(b))
})

test_that("dMcast air quality tests", {
  melt <- function(data, idColumns) {
    cols <- setdiff(colnames(data), idColumns)
    results <- lapply(cols, function(x) cbind(data[, idColumns], variable = x, value = as.numeric(data[, x])))
    results <- Reduce(rbind, results)
  }
  names(airquality) <- tolower(names(airquality))
  aqm <- melt(airquality, idColumns = c("month", "day"))
  a <- dMcast(aqm, month:day ~ variable, fun.aggregate = "mean", value.var = "value")
  expect_equal(nrow(unique(aqm[, c("month", "day")])), nrow(a))
  aqm[aqm == 30] <- NA
  aqm$day <- as.factor(aqm$day)
  b <- dMcast(aqm, month:day ~ variable, fun.aggregate = "mean", value.var = "value")
  expect_equal(nrow(unique(aqm[, c("month", "day")])), nrow(b))
  aqm[aqm == "temp"] <- NA
  c <- dMcast(aqm, month:day ~ variable, fun.aggregate = "mean", value.var = "value")
  expect_equal(length(unique(aqm$variable)), ncol(c))
  d <- dMcast(aqm, month:day ~ variable, fun.aggregate = "mean", value.var = "value", factor.nas = FALSE)
  expect_equal(ncol(d), 3)
})


test_that("aggregate.Matrix works", {

  # fake sparse matrix
  counts <- matrix(c(1,0,0,0,1,
                     2,0,0,1,1,
                     1,1,1,1,1),
                   nrow = 3,
                   byrow = TRUE,
                   dimnames = list(c('a','a', 'b'), NULL))

  expected <- matrix(c(3,0,0,1,2,
                       1,1,1,1,1),
                     nrow = 2,
                     byrow = TRUE,
                     dimnames = list(c('a', 'b'), NULL))

  res <- aggregate.Matrix(counts, row.names(counts), fun = 'sum')

  expect_equal(as.matrix(res), expected)
})
