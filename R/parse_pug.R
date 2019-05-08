#' Check Pubchem PUG-View for GRAS notice
#'
#' @param pug_view List from reading PUG-View JSON data with \code{\link[rjson]{fromJSON}}
#'
#' @return Boolean indicating if \code{pug_view} has a GRAS notice.
#' @export
#'
#' @examples
#'
#' pug_view <- fromJSON(file='/home/alex/Documents/Batcave/zaklab/drugseqr/data-raw/drug_annot/pug_view/views/4091.json')
#'
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

#' Extract DrugBank URL from PUG-View
#'
#' @inheritParams check_gras
#'
#' @return Character vector with DrugBank URL or NA if none exists
#' @export
#'
#' @examples
get_drugbank <- function(pug_view) {

  # default
  db_url <- NA_character_

  # get reference source names
  refs <- pug_view$Record$Reference
  sources <- sapply(refs, `[[`, 'SourceName')

  # get first drugbank index
  dbi <- which(sources == 'DrugBank')[1]

  if (!is.na(dbi)) {
    # get base url
    db_url <- refs[[dbi]]$URL
    db_url <- gsub('^([^#]+/DB[0-9]+)#.+?$', '\\1', db_url)
  }

  return(db_url)
}
