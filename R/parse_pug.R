#' Check Pubchem PUG-View for GRAS notice
#'
#' @param pug_view List from reading PUG-View JSON data with \code{\link[rjson]{fromJSON}}
#'
#' @return Boolean indicating if \code{pug_view} has a GRAS notice.
#' @export
#'
#' @examples
check_gras <- function(pug_view) {

  # default
  is_gras <- FALSE

  # get top level headings
  sections <- pug_view$Record$Section
  headings <- sapply(sections, `[[`, 'TOCHeading')

  # check for index of food heading
  foodi <- which(headings == 'Food Additives and Ingredients')

  if (length(foodi)) {
    food_sections <- sections[[foodi]]$Section
    food_headings <- sapply(food_sections, `[[`, 'TOCHeading')
    is_gras <- 'FDA Generally Recognized as Safe - GRAS Notices' %in% food_headings
  }

  return(is_gras)
}
