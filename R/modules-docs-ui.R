#' UI for docs page (add new sections here)
#' @export
#' @keywords internal

#' UI for docs page (add new sections here)
#' @export
#' @keywords internal
docsPageUI <- function(id, tab, active) {
  docsSections <- list(
    # Datasets documentation
    docsSection(id = 'ds-docs', name = 'Datasets',
                content = list(
                  docsSubsection(id = 'ds-quantification',
                                 name = 'Quantifying RNA-seq datasets',
                                 content = tagList(
                                   p("To add a new dataset, first make sure the dataset name input is empty. Then type
                                   a name for your new bulk or single-cell dataset and hit enter. You will be prompted to select a folder
                                   containing either fastq.gz files (bulk or single-cell) or CellRanger files (single-cell)."),
                                   p("If a folder with bulk fastq.gz files is selected, you will be prompted to confirm the auto-selected end type.
                                     If the end type is pair-ended, you will also need to indicate pairs. To do this, click the rows corresponding to the
                                     paired files and click Paired. If present, you should also indicate any files that are replicates of the same sample before you run quantification."),
                                   p("If a folder with single-cell fastq.gz files is selected, nothing further is necessary before running quantification. Quantification uses kallisto for both bulk and single-cell RNA-Seq data.
                                     After quantification completes, bulk datasets can be selected in the Datasets tab for differential expression and exploratory analyses whereas single-cell datasets
                                     will show up in the Single Cell tab for exploratory analyses and sample integration.")
                                 )

                  ),
                  docsSubsection(id = 'ds-differential-expression',
                                 name = 'Differential expression',
                                 content = tagList(
                                   p("For quantified bulk RNA-Seq datasets, select the differential expression analysis type in order to run
                                     differential expression analyses between two groups. Any differential expression analyses performed will be accessible
                                     in the Pathways and Drugs tab for further analysis."),
                                   p("")
                                 )
                  ),
                  docsSubsection(id = 'ds-gene-plots',
                                 name = 'Sample-level gene expression',
                                 content = tagList(p(shinipsum::random_text(nwords=100)),
                                                   p(shinipsum::random_text(nwords=100)))
                  ),
                  docsSubsection(id = 'ds-cell-plots',
                                 name = 'Cell-type deconvolution',
                                 content = tagList(p(shinipsum::random_text(nwords=100)),
                                                   p(shinipsum::random_text(nwords=100)))
                  )
                )),

    # Single cell documentation
    docsSection(id = 'sc-docs', name = 'Single Cell',
                content = list(
                  docsSubsection(id = 'sc-label-clusters',
                                 name = 'Labeling cell clusters',
                                 content = tagList(p(shinipsum::random_text(nwords=100)),
                                                   p(shinipsum::random_text(nwords=100)))

                  ),
                  docsSubsection(id = 'sc-integrate-samples',
                                 name = 'Integrate datasets',
                                 content = tagList(p(shinipsum::random_text(nwords=100)),
                                                   p(shinipsum::random_text(nwords=100)))
                  ),
                  docsSubsection(id = 'sc-differential-expression',
                                 name = 'Differential expression analysis',
                                 content = tagList(p(shinipsum::random_text(nwords=100)),
                                                   p(shinipsum::random_text(nwords=100)))
                  )
                )),

    # Pathways documentation
    docsSection(id = 'path-docs', name = 'Pathways',
                content = list(
                  docsSubsection(id = 'path-overview',
                                 name = 'Overview',
                                 content = tagList(p(shinipsum::random_text(nwords=100)),
                                                   p(shinipsum::random_text(nwords=100)))

                  ),
                  docsSubsection(id = 'path-query',
                                 name = 'CMAP02/L1000 query genes',
                                 content = tagList(p(shinipsum::random_text(nwords=100)),
                                                   p(shinipsum::random_text(nwords=100)))
                  )
                )),

    # Drugs documentation
    docsSection(id = 'drugs-docs', name = 'Drugs',
                content = list(
                  docsSubsection(id = 'drugs-overview',
                                 name = 'Overview',
                                 content = tagList(p(shinipsum::random_text(nwords=100)),
                                                   p(shinipsum::random_text(nwords=100)))

                  ),
                  docsSubsection(id = 'drugs-custom-queries',
                                 name = 'Custom queries',
                                 content = tagList(p(shinipsum::random_text(nwords=100)),
                                                   p(shinipsum::random_text(nwords=100)))
                  ),
                  docsSubsection(id = 'drugs-advanced',
                                 name = 'Advanced options',
                                 content = tagList(p(shinipsum::random_text(nwords=100)),
                                                   p(shinipsum::random_text(nwords=100)))
                  )
                ))
  )


  section_info <- extract_section_info(docsSections)


  withTags({
    tabPane(tab, active,
            div(class = "bs-docs-container",
                div(class = 'row',
                    div(class = 'col-md-6 col-md-offset-1', role="main", docsSections),
                    div(class = "col-md-4 col-md-offset-1 bd-toc", role="complementary", docsSideNav(section_info))
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

