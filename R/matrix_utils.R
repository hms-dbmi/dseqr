#' Compute summary statistics of a Matrix
#'
#' Similar to \code{\link[stats]{aggregate}}.  Splits the matrix into groups as
#' specified by groupings, which can be one or more variables. Aggregation
#' function will be applied to all columns in data, or as specified in formula.
#' Warning: groupings will be made dense if it is sparse, though data will not.
#'
#' \code{aggregate.Matrix} uses its own implementations of functions and should
#' be passed a string in the \code{fun} argument.
#'
#' @param x a \code{\link{Matrix}} or matrix-like object
#' @param groupings an object coercible to a group of factors defining the
#'   groups
#' @param form \code{\link[stats]{formula}}
#' @param fun character string specifying the name of aggregation function to be
#'   applied to all columns in data.  Currently "sum", "count", and "mean"
#'   are supported.
#' @param ... arguments to be passed to or from methods.  Currently ignored
#' @return A sparse \code{Matrix}.  The rownames correspond to the values of the
#'   groupings or the interactions of groupings joined by a \code{_}.
#'
#'   There is an attribute \code{crosswalk} that includes the groupings as a
#'   data frame.  This is necessary because it is not possible to include
#'   character or data frame groupings in a sparse Matrix.  If needed, one can
#'   \code{cbind(attr(x,"crosswalk"),x)} to combine the groupings and the
#'   aggregates.
#'
#' @seealso \code{\link[dplyr]{summarise}}
#' @seealso \code{\link[plyr]{summarise}}
#' @seealso \code{\link[stats]{aggregate}}
#' @export
#' @export aggregate.Matrix
#' @examples
#' skus <- Matrix::Matrix(as.matrix(data.frame(
#'   orderNum = sample(1000, 10000, TRUE),
#'   sku = sample(1000, 10000, TRUE),
#'   amount = runif(10000)
#' )), sparse = TRUE)
#' # Calculate sums for each sku
#' a <- aggregate.Matrix(skus[, "amount"], skus[, "sku"], fun = "sum")
#' # Calculate counts for each sku
#' b <- aggregate.Matrix(skus[, "amount"], skus[, "sku"], fun = "count")
#' # Calculate mean for each sku
#' c <- aggregate.Matrix(skus[, "amount"], skus[, "sku"], fun = "mean")
#'
#' m <- Matrix::rsparsematrix(1000000, 100, .001)
#' labels <- as.factor(sample(1e4, 1e6, TRUE))
#' b <- aggregate.Matrix(m, labels)
#'
#' \dontrun{
#' orders <- data.frame(
#'   orderNum = as.factor(sample(1e6, 1e7, TRUE)),
#'   sku = as.factor(sample(1e3, 1e7, TRUE)),
#'   customer = as.factor(sample(1e4, 1e7, TRUE)),
#'   state = sample(letters, 1e7, TRUE), amount = runif(1e7)
#' )
#' system.time(d <- aggregate.Matrix(orders[, "amount", drop = FALSE], orders$orderNum))
#' system.time(e <- aggregate.Matrix(orders[, "amount", drop = FALSE], orders[, c("customer", "state")]))
#' }
aggregate.Matrix <- function(x, groupings = NULL, form = NULL, fun = "sum", ...) {
  if (!methods::is(x, "Matrix")) {
    x <- Matrix::Matrix(as.matrix(x), sparse = TRUE)
  }
  if (fun == "count") {
    x <- x != 0
  }

  groupings2 <- groupings
  if (!methods::is(groupings2, "data.frame")) {
    groupings2 <- as.data.frame(groupings2)
  }

  groupings2 <- data.frame(lapply(groupings2, as.factor))
  groupings2 <- data.frame(interaction(groupings2, sep = "_"))
  colnames(groupings2) <- "A"

  if (is.null(form)) {
    form <- stats::as.formula("~0+.")
  }

  form <- stats::as.formula(form)
  mapping <- dMcast(groupings2, form)
  colnames(mapping) <- substring(colnames(mapping), 2)
  result <- Matrix::t(mapping) %*% x
  if (fun == "mean") {
    result@x <- result@x / (aggregate.Matrix(x, groupings2, fun = "count"))@x
  }
  attr(result, "crosswalk") <- grr::extract(groupings, match(rownames(result), groupings2$A))
  return(result)
}


#' Casts or pivots a long \code{data frame} into a wide sparse matrix
#'
#' Similar in function to \code{\link[reshape2]{dcast}}, but produces a sparse
#' \code{\link{Matrix}} as an output. Sparse matrices are beneficial for this
#' application because such outputs are often very wide and sparse. Conceptually
#' similar to a \code{pivot} operation.
#'
#' Casting formulas are slightly different than those in \code{dcast} and follow
#' the conventions of \code{\link{model.matrix}}. See \code{\link{formula}} for
#' details.  Briefly, the left hand side of the \code{~} will be used as the
#' grouping criteria.  This can either be a single variable, or a group of
#' variables linked using \code{:}.  The right hand side specifies what the
#' columns will be. Unlike \code{dcast}, using the \code{+} operator will append
#' the values for each variable as additional columns.  This is useful for
#' things such as one-hot encoding.  Using \code{:} will combine the columns as
#' interactions.
#'
#' @param data a data frame
#' @param formula casting \code{\link[stats]{formula}}, see details for specifics.
#' @param fun.aggregate name of aggregation function.  Defaults to 'sum'
#' @param value.var name of column that stores values to be aggregated numerics
#' @param as.factors if TRUE, treat all columns as factors, including
#' @param factor.nas if TRUE, treat factors with NAs as new levels.  Otherwise,
#'  rows with NAs will receive zeroes in all columns for that factor
#' @param drop.unused.levels should factors have unused levels dropped? Defaults to TRUE,
#'  in contrast to \code{\link{model.matrix}}
#' @return a sparse \code{Matrix}
#' @seealso \code{\link[reshape]{cast}}
#' @seealso \code{\link[reshape2]{dcast}}
#' @export
#' @examples
#' # Classic air quality example
#' melt <- function(data, idColumns) {
#'   cols <- setdiff(colnames(data), idColumns)
#'   results <- lapply(cols, function(x) cbind(data[, idColumns], variable = x, value = as.numeric(data[, x])))
#'   results <- Reduce(rbind, results)
#' }
#' names(airquality) <- tolower(names(airquality))
#' aqm <- melt(airquality, idColumns = c("month", "day"))
#' dMcast(aqm, month:day ~ variable, fun.aggregate = "mean", value.var = "value")
#' dMcast(aqm, month ~ variable, fun.aggregate = "mean", value.var = "value")
#'
#' # One hot encoding
#' # Preserving numerics
#' dMcast(warpbreaks, ~.)
#' # Pivoting numerics as well
#' dMcast(warpbreaks, ~., as.factors = TRUE)
#'
#' \dontrun{
#' orders <- data.frame(
#'   orderNum = as.factor(sample(1e6, 1e7, TRUE)),
#'   sku = as.factor(sample(1e3, 1e7, TRUE)),
#'   customer = as.factor(sample(1e4, 1e7, TRUE)),
#'   state = sample(letters, 1e7, TRUE),
#'   amount = runif(1e7)
#' )
#' # For simple aggregations resulting in small tables, dcast.data.table (and reshape2) will be faster
#' system.time(a <- dcast.data.table(as.data.table(orders), sku ~ state, sum,
#'   value.var = "amount"
#' )) # .5 seconds
#' system.time(b <- reshape2::dcast(orders, sku ~ state, sum,
#'   value.var = "amount"
#' )) # 2.61 seconds
#' system.time(c <- dMcast(orders, sku ~ state,
#'   value.var = "amount"
#' )) # 8.66 seconds
#'
#' # However, this situation changes as the result set becomes larger
#' system.time(a <- dcast.data.table(as.data.table(orders), customer ~ sku, sum,
#'   value.var = "amount"
#' )) # 4.4 seconds
#' system.time(b <- reshape2::dcast(orders, customer ~ sku, sum,
#'   value.var = "amount"
#' )) # 34.7 seconds
#' system.time(c <- dMcast(orders, customer ~ sku,
#'   value.var = "amount"
#' )) # 14.55 seconds
#'
#' # More complicated:
#' system.time(a <- dcast.data.table(as.data.table(orders), customer ~ sku + state, sum,
#'   value.var = "amount"
#' )) # 16.96 seconds, object size = 2084 Mb
#' system.time(b <- reshape2::dcast(orders, customer ~ sku + state, sum,
#'   value.var = "amount"
#' )) # Does not return
#' system.time(c <- dMcast(orders, customer ~ sku:state,
#'   value.var = "amount"
#' )) # 21.53 seconds, object size = 116.1 Mb
#'
#' system.time(a <- dcast.data.table(as.data.table(orders), orderNum ~ sku, sum,
#'   value.var = "amount"
#' )) # Does not return
#' system.time(c <- dMcast(orders, orderNum ~ sku,
#'   value.var = "amount"
#' )) # 24.83 seconds, object size = 175Mb
#'
#' system.time(c <- dMcast(orders, sku:state ~ customer,
#'   value.var = "amount"
#' )) # 17.97 seconds, object size = 175Mb
#' }
dMcast <- function(data, formula, fun.aggregate = "sum", value.var = NULL, as.factors = FALSE, factor.nas = TRUE, drop.unused.levels = TRUE) {
  values <- 1
  if (!is.null(value.var)) {
    values <- data[, value.var]
  }
  alltms <- stats::terms(formula, data = data)
  response <- rownames(attr(alltms, "factors"))[attr(alltms, "response")]
  tm <- attr(alltms, "term.labels")
  interactionsIndex <- grep(":", tm)
  interactions <- tm[interactionsIndex]
  simple <- setdiff(tm, interactions)
  i2 <- strsplit(interactions, ":")
  newterms <- unlist(lapply(i2, function(x) paste("paste(", paste(x, collapse = ","), ",", "sep='_'", ")")))
  newterms <- c(simple, newterms)
  newformula <- stats::as.formula(paste("~0+", paste(newterms, collapse = "+")))
  allvars <- all.vars(alltms)
  data <- data[, c(allvars), drop = FALSE]
  if (as.factors) {
    data <- data.frame(lapply(data, as.factor))
  }
  characters <- unlist(lapply(data, is.character))
  data[, characters] <- lapply(data[, characters, drop = FALSE], as.factor)
  factors <- unlist(lapply(data, is.factor))
  # Prevents errors with 1 or fewer distinct levels
  data[, factors] <- lapply(data[, factors, drop = FALSE], function(x) {
    if (factor.nas) {
      if (any(is.na(x))) {
        levels(x) <- c(levels(x), "NA")
        x[is.na(x)] <- "NA"
      }
    }
    if (drop.unused.levels) {
      if (nlevels(x) != length(stats::na.omit(unique(x)))) {
        x <- factor(as.character(x))
      }
    }
    y <- stats::contrasts(x, contrasts = FALSE, sparse = TRUE)
    attr(x, "contrasts") <- y
    return(x)
  })
  # Allows NAs to pass
  attr(data, "na.action") <- stats::na.pass
  result <- Matrix::sparse.model.matrix(newformula, data, drop.unused.levels = FALSE, row.names = FALSE)
  brokenNames <- grep("paste(", colnames(result), fixed = TRUE)
  colnames(result)[brokenNames] <- lapply(colnames(result)[brokenNames], function(x) {
    x <- gsub("paste(", replacement = "", x = x, fixed = TRUE)
    x <- gsub(pattern = ", ", replacement = "_", x = x, fixed = TRUE)
    x <- gsub(pattern = '_sep = \"_\")', replacement = "", x = x, fixed = TRUE)
    return(x)
  })

  result <- result * values
  if (isTRUE(response > 0)) {
    responses <- all.vars(stats::terms(stats::as.formula(paste(response, "~0"))))
    result <- aggregate.Matrix(result, data[, responses, drop = FALSE], fun = fun.aggregate)
  }
  return(result)
}
