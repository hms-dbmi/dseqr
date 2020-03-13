#' Read RDS file safely
#'
#' @inheritParams base::readRDS
#' @param .nofile return value when file doesn't exist
#' @param .nullfile return value when file value is \code{NULL}
#'
#' @inherit base::readRDS return
#' @export
#' @keywords internal
#'
readRDS.safe <- function(file, .nofile = NULL, .nullfile = NULL) {
  res <- .nofile
  if (isTruthy(path) && file.exists(path))
    res <- readRDS(path)

  if (is.null(res)) return(.nullfile)
  return(res)
}


#' Disable multiple ids
#'
#' @param ids Character vector of ids to disable
#' @export
#' @keywords internal
disableAll <- function(ids){
  for (id in ids) shinyjs::disable(id)
}


#' Enable multiple ids
#'
#' @param ids Character vector of ids to enable
#' @export
#' @keywords internal
enableAll <- function(ids) {
  for (id in ids) shinyjs::enable(id)
}

#' Check truthiness of multiple objects
#'
#' @param ... objects to check truthiness of
#' @export
#' @keywords internal
isTruthyAll <- function(...) {
  x <- list(...)
  for (xi in x) if (!shiny::isTruthy(xi)) return(FALSE)
  return(TRUE)
}
