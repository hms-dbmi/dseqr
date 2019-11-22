#' UI for docs page (add new sections here)
#' @export
#' @keywords internal
docsPageUI <- function(id, tab, active) {
  ns <- NS(id)

  docsSections <- list(
    # Overview section
    docsSection(id = 'js-overview', name = 'Overview',
                content = list(

                  docsSubsection(id = 'js-individual-compiled',
                                 name = 'Individual or compiled',
                                 content = tagList(p(shinipsum::random_text(nwords=500)),
                                                   p(shinipsum::random_text(nwords=500)))
                  ),
                  docsSubsection(id = 'js-data-attrs',
                                 name = 'Data attributes',
                                 content = tagList(p(shinipsum::random_text(nwords=500)),
                                                   p(shinipsum::random_text(nwords=500)))
                  ),
                  docsSubsection(id = 'js-programmatic-api',
                                 name = 'Programmatic API',
                                 content = tagList(p(shinipsum::random_text(nwords=500)),
                                                   p(shinipsum::random_text(nwords=500)))
                  ),
                  docsSubsection(id = 'js-noconflict',
                                 name = 'No conflict',
                                 content = tagList(p(shinipsum::random_text(nwords=500)),
                                                   p(shinipsum::random_text(nwords=500)))
                  )
                ))
  )


  section_info <- extract_section_info(docsSections)


  withTags({
    tabPane(tab, active,
            div(class = "container bs-docs-container",
                div(class = 'row',
                    div(class = 'col-md-8', role="main", docsSections),
                    div(class = "col-md-3 col-md-offset-1 bd-toc", role="complementary", docsSideNav(section_info))
                )
            )
    )
  })
}


#' A primary section in the Docs page
#' @export
#' @keywords internal
docsSection <- function(id, name, content) {

  withTags({
    div(class = 'bs-docs-section',
        h1(id = id, class='page-header',
           tags$a(href = paste0('#', id)),
           name
        ),
        content
    )
  })

}

#' A subsections of a docsSection
#' @export
#' @keywords internal
docsSubsection <- function(id, name, content) {
  withTags({
    tagList(
      h2(id = id,
         tags$a(href = paste0('#', id)),
         name

      ),
      content
    )
  })
}

#' Navigation on right side of Docs page
#' @export
#' @keywords internal
docsSideNav <- function(section_info) {

  section_id_names <- section_info$section_id_names
  subsection_id_names <- section_info$subsection_id_names

  nsections <- length(section_id_names)

  withTags({
    nav(class = "bs-docs-sidebar well-form well-bg",
        ul(class="nav bs-docs-sidenav",
           lapply(seq_len(nsections), function(i) {

             # get id and name for current section
             section_id_name <- section_id_names[[i]]

             # get id and names list of subsections for current section
             subsectioni_id_names <- subsection_id_names[[i]]

             li(
               # link for current section
               a(href=paste0("#", section_id_name[1]), section_id_name[2]),
               ul(class="nav",
                  # link for each subsection of current section
                  lapply(subsectioni_id_names, function(subsectioni_id_name) {
                    li(
                      a(href = paste0("#", subsectioni_id_name[1]), subsectioni_id_name[2])
                    )
                  }
                  )
               )
             )
           })
        ))
  })
}


#' Extracts info needed to construct docsSideNav
#' @export
#' @keywords internal
extract_section_info <- function(docsSections) {
  section_id_names <- lapply(docsSections, function(docsSection) {

    # title is in first child
    # get id and name
    id <- docsSection$children[[1]]$attribs$id
    name <- docsSection$children[[1]]$children[[2]]

    return(c(id, name))
  })

  subsection_id_names <- lapply(docsSections, function(docsSection) {

    # subsections are in the second child
    subsections <- docsSection$children[[2]]

    lapply(subsections, function(subsection) {
      # subsection titles are in the first child
      id <- subsection[[1]]$attribs$id
      name <- subsection[[1]]$children[[2]]

      return(c(id, name))
    })

  })

  return(list(
    section_id_names = section_id_names,
    subsection_id_names = subsection_id_names
  ))
}

