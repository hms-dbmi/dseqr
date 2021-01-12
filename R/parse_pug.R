#' Check Pubchem PUG-View for GRAS notice
#'
#' @param pug_view List from reading PUG-View JSON data with \code{\link[rjson]{fromJSON}}
#'
#' @return Boolean indicating if \code{pug_view} has a GRAS notice.
#' @keywords internal
#'
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
#' @return Character vector with DrugBank ID or NA if none exists
#' @keywords internal
get_drugbank <- function(pug_view) {

  # default
  db_id <- NA_character_

  # get reference source names
  refs <- pug_view$Record$Reference
  sources <- sapply(refs, `[[`, 'SourceName')

  # get first drugbank index
  dbi <- which(sources == 'DrugBank')[1]

  if (!is.na(dbi)) {
    # get base url
    db_url <- refs[[dbi]]$URL
    db_id <- gsub('^[^#]+/(DB[0-9]+)', '\\1', db_url)
  }
  return(db_id)
}


#' Extract Wikipedia URL from PUG-View
#'
#' @inheritParams check_gras
#'
#' @return Character vector with Wikipedia page
#' @keywords internal
get_wikipedia <- function(pug_view) {

  # default
  wiki_name <- NA_character_

  # get reference source names
  refs <- pug_view$Record$Reference
  sources <- sapply(refs, `[[`, 'SourceName')

  # get first drugbank index
  wiki <- which(sources == 'Wikipedia')[1]

  if (!is.na(wiki)) {
    # get base url
    wiki_url <- refs[[wiki]]$URL

    # make sure is wikipedia (wikidata like Pubchem but worse)
    if (grepl('https://en.wikipedia.org/wiki', wiki_url, fixed = TRUE))
      wiki_name <- gsub('^https://en.wikipedia.org/wiki/([^/]+)$', '\\1', wiki_url)
  }
  return(wiki_name)
}
