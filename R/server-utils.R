#' Read qs file safely
#'
#' @inheritParams qs::qread
#' @param .nofile return value when file doesn't exist
#' @param .nullfile return value when file value is \code{NULL}
#'
#' @inherit qs::qread return
#'
#' @export
#'
qread.safe <- function(file, .nofile = NULL, .nullfile = NULL) {
  res <- .nofile
  if (isTruthy(file) && file_exists(file))
    res <- tryCatch(qs::qread(file), error = function(e) NULL)

  if (is.null(res)) return(.nullfile)
  return(res)
}


#' Disable multiple ids
#'
#' @param ids Character vector of ids to disable
#' @return called for side effects
#'
#' @keywords internal
disableAll <- function(ids, asis = rep(FALSE, length(ids))){
  for (i in seq_along(ids)) shinyjs::disable(ids[i], asis = asis[i])
  shinyjs::runjs("$(\".tooltip\").tooltip(\"hide\");")
}


#' Enable multiple ids
#'
#' @param ids Character vector of ids to enable
#'
#' @return No return value. Called for side effects.
#' @keywords internal
enableAll <- function(ids, asis = rep(FALSE, length(ids))) {
  for (i in seq_along(ids)) shinyjs::enable(ids[i], asis = asis[i])
}

#' Check truthiness of multiple objects
#'
#' @param ... objects to check truthiness of
#'
#' @keywords internal
isTruthyAll <- function(...) {
  x <- list(...)
  for (xi in x) if (!isTruthy(xi)) return(FALSE)
  return(TRUE)
}
