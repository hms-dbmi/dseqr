#' UI for navbar
#' @param active the active tab name
#' @export
#' @keywords internal
navbarUI <- function(tabs, active) {

  withTags({
    nav(class = 'navbar navbar-default navbar-static-top',
        div(class = 'container-fluid',
            div(class = 'navbar-header',
                span(class = 'navbar-brand', title = 'drugseqr',
                     span(class = 'brand-icons',
                          i(class = 'glyphicon glyphicon-leaf'),
                          'drugseqr'
                     )
                )
            ),
            ul(class = 'nav navbar-nav', `data-tabsetid` = 'tabset',
               lapply(seq_along(tabs), function(i) {
                 tab <- tabs[i]
                 is.active <- tab == active
                 li(class = ifelse(is.active, 'active', ''),
                    a(href = paste0('#', id_from_tab(tab)), `data-toggle` = 'tab', `data-value` = tab, `aria-expanded` = ifelse(is.active, 'true', 'false'), tab)
                 )
               })
            )
        )
    )
  })
}


#' UI for a tab pane
#'
#' @param tab The name of the tab
#' @param active The name of the active tab
#' @param ... The UI elements to place in the tab
#' @return shiny div tag with UI for tab
#' @export
#' @keywords internal
tabPane <- function(tab, active, ...) {
  active_class <- ifelse(tab == active, 'active', '')
  tags$div(class = paste('tab-pane', active_class), `data-value` = tab, id = id_from_tab(tab), ...)
}


#' Get cluster choices data.frame for selectize dropdown
#'
#' @param clusters Character vector of cluster names
#' @param scseq \code{Seurat} object
#'
#' @return data.frame with columns for rendering selectizeInput cluster choices
#' @export
#' @keywords internal
get_cluster_choices <- function(clusters, scseq) {

  # show the cell numbers/percentages
  ncells <- tabulate(scseq$seurat_clusters)
  pcells <- round(ncells / sum(ncells) * 100)
  pspace <- strrep('&nbsp;&nbsp;', 2 - nchar(pcells))

  # cluster choices are the clusters themselves
  testColor <- get_palette(clusters)
  cluster_choices <- data.frame(name = stringr::str_trunc(clusters, 27),
                                value = clusters,
                                label = clusters,
                                testColor,
                                ncells, pcells, pspace, row.names = NULL)

  return(cluster_choices)
}

#' Get contrast choices data.frame for selectize dropdown
#'
#' @param clusters Character vector of cluster names
#' @param test Name of test contrast
#'
#' @return data.frame with columns for rendering selectizeInput contrast choices
#' @export
#' @keywords internal
get_contrast_choices <- function(clusters, test) {

  # group choices are as compared to other clusters
  ctrls <- clusters[clusters != test]

  colours <- get_palette(clusters)
  names(colours) <- clusters

  contrast_choices <- data.frame(test = stringr::str_trunc(test, 11, ellipsis = '..'),
                                 ctrl = stringr::str_trunc(c('all', ctrls), 11, ellipsis = '..'),
                                 value = c(test, paste0(test, ' vs ', ctrls)),
                                 testColor = colours[test],
                                 ctrlColor = c('white', colours[ctrls]), row.names = NULL)

  return(contrast_choices)

}


#' Convert tab name to formated id
#'
#' used by navbarUI and *PageUI for drugseqr app
#'
#' @param tab The name of the tab (e.g. \code{'Single Cell'})
#' @export
#' @keywords internal
id_from_tab <- function(tab) {
  id <- tolower(tab)
  gsub(' ', '-', id)
}
