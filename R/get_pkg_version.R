#' Get version of salmon/kallisto from system command.
#'
#' @param type Either \code{'salmon'} or \code{'kallisto'}.
#'
#' @return Version of salmon/kallisto.
#' @export
#' @keywords internal
#'
get_pkg_version <- function(type) {
  # possibly use older salmon with version appended to executable name
  if (type == 'salmon') {
    version <- system(paste(type, '--version'), intern = TRUE)
    version <- gsub('^salmon ', '', version)

  } else if (type == 'kallisto') {
    version <- system(paste(type, 'version'), intern = TRUE)
    version <- gsub('^kallisto, version ', '', version)

  } else {
    stop('type must be either salmon or kallisto')
  }

  return(version)
}
