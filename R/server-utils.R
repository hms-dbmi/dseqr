#' Read qs file safely
#'
#' @inheritParams qs::qread
#' @param .nofile return value when file doesn't exist
#' @param .nullfile return value when file value is \code{NULL}
#'
#' @inherit qs::qread return
#'
#' @keywords internal
#'
qread.safe <- function(file, .nofile = NULL, .nullfile = NULL) {
  res <- .nofile
  if (isTruthy(file) && file.exists(file))
    res <- qs::qread(file)

  if (is.null(res)) return(.nullfile)
  return(res)
}


#' Disable multiple ids
#'
#' @param ids Character vector of ids to disable
#'
#' @keywords internal
disableAll <- function(ids){
  for (id in ids) shinyjs::disable(id)
  shinyjs::runjs("$(\".tooltip\").tooltip(\"hide\");")
}


#' Enable multiple ids
#'
#' @param ids Character vector of ids to enable
#'
#' @keywords internal
enableAll <- function(ids) {
  for (id in ids) shinyjs::enable(id)
}

#' Check truthiness of multiple objects
#'
#' @param ... objects to check truthiness of
#'
#' @keywords internal
isTruthyAll <- function(...) {
  x <- list(...)
  for (xi in x) if (!shiny::isTruthy(xi)) return(FALSE)
  return(TRUE)
}
