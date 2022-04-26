



#' UI for navbar
#'
#' @param tabs Character vector of tab names to display
#' @param active the active tab name
#' @param logout_url url to logout from. NULL results in no Logout link
#' @return shiny.tag with html for navbar
#'
#' @export
#' @examples
#'
#' # UI with only Single-Cell page
#' navbarUI('Single Cell', 'Single Cell')
#'
navbarUI <- function(tabs, active, logout_url = NULL) {

  logout_li <- NULL
  if (!is.null(logout_url)) {
    # target _top break out of iframe for ShinyProxy
    logout_li <- tags$li(class = 'navbar-right', `data-toggle`="collapse",
                         `data-target`=".navbar-collapse.in", a(target="_top", href = logout_url, 'Sign out')
    )

  }

  tags$nav(class = 'navbar navbar-default navbar-static-top',
           div(class = 'container-fluid',
               div(class = 'navbar-header',
                   tags$button(type='button', class='navbar-toggle collapsed', `data-toggle`='collapse', `data-target`='#bs-navbar', `aria-expanded`='false',
                               span(class = 'sr-only', 'Toggle navigation'),
                               span(class = 'icon-bar'),
                               span(class = 'icon-bar'),
                               span(class = 'icon-bar')
                   ),
                   span(class = 'navbar-brand', title = 'dseqr',
                        tags$a(class = 'brand-icons',
                               href="/",
                               tags$img(src="favicon.png"),
                               span('seqr')
                        )
                   )
               ),
               div(id = 'bs-navbar', class = 'collapse navbar-collapse',
                   tags$ul(class = 'nav navbar-nav shiny-tab-input shiny-bound-input', `data-tabsetid` = 'tabset', id = 'tab',
                           lapply(seq_along(tabs), function(i) {
                             tab <- tabs[i]
                             is.active <- tab == active
                             tags$li(class = ifelse(is.active, 'active', ''), `data-toggle`="collapse", `data-target`=".navbar-collapse.in",
                                     a(href = paste0('#', id_from_tab(tab)), `data-toggle` = 'tab', `data-value` = tab, `aria-expanded` = ifelse(is.active, 'true', 'false'), tab)
                             )
                           }),
                           # github linkout section
                           tags$li(class = 'navbar-right', `data-toggle`="collapse", `data-target`=".navbar-collapse.in",
                                   a(href = 'https://github.com/hms-dbmi/dseqr', target="_blank", icon('github'), style = 'padding-bottom: 0px; font-size:17px;',)
                           ),
                           # docs section
                           tags$li(class = 'navbar-right', `data-toggle`="collapse", `data-target`=".navbar-collapse.in", id = 'docs-link',
                                   a(href = "https://docs.dseqr.com/docs/single-cell/add-dataset/", `target` = '_blank', 'Docs')
                           ),
                           # logout section
                           logout_li
                   )
               )
           )
  )
}


#' UI for secondary navbar
#'
#' @param hide should delete and add be disabled (for demo)
#' @return shiny.tag with html for secondary navbar
#'
#' @export
#' @examples
#'
#' navbar2UI(hide = FALSE)
#'
navbar2UI <- function(hide) {
  class <- 'action-button shiny-bound-input btn-intro btn navbar-btn'
  add <- ifelse(hide, 'disabled', '')

  ui <- tags$nav(
    class = 'navbar navbar-default secondary-navbar navbar-expand',
    tags$div(
      class = 'container-fluid',
      tags$div(
        tags$ul(
          class = 'nav navbar-non-responsive',
          tags$div(class = "secondary-navbar-btn-group pull-right",


                   # Tour Dropdown for Single Cell Tab
                   tags$div(class='btn-group',
                            tags$button(
                              id = 'tour_dropdown',
                              class = paste0(class, ' dropdown-toggle'),
                              `data-toggle` = 'dropdown',
                              type = 'button',
                              haspopup = 'true',
                              `aria-expanded`='false',
                              icon('info', 'fa-fw'),
                              'Tour',
                              tags$span(class='caret')),
                            tags$ul(class="dropdown-menu",
                                    tags$li(class = 'dropdown-header', 'Tour Comparison Between:'),
                                    tags$li(tags$a(
                                      id = 'tour_sc_clusters',
                                      role='button',
                                      class = 'action-button shiny-bound-input',
                                      'Clusters')
                                    ),
                                    tags$li(tags$a(
                                      id = 'tour_sc_samples',
                                      role='button',
                                      class = 'action-button shiny-bound-input',
                                      'Samples')
                                    )
                            )
                   ),

                   # Tour Button
                   tags$li(id = 'start_tour_container',
                           style = 'display: none',
                           tags$button(class = 'btn',
                                       id = 'start_tour',
                                       class=class,
                                       icon('info', 'fa-fw'),
                                       'Tour')
                   ),

          ),

          # Projects
          tags$button(class = 'btn',
                      id = 'select_project',
                      class=class,
                      tags$i(class= 'far fa-folder-open fa-fw'),
                      tags$span(class = 'hidden-xs', 'Projects')),

          tags$i(class = 'fa fa-angle-right', style = 'padding: 0 5px; color: white;'),

          # Dataset Management dropdown
          tags$div(class='btn-group',
                   tags$button(
                     id = 'datasets_dropdown',
                     class = paste0(class, ' dropdown-toggle'),
                     `data-toggle` = 'dropdown',
                     type = 'button',
                     haspopup = 'true',
                     `aria-expanded`='false',
                     tags$i(class= 'far fa-file fa-fw'),
                     tags$span(class = 'hidden-xs', 'Dataset Management'),
                     tags$span(class='caret')),
                   tags$ul(class="dropdown-menu",
                           tags$li(class=add,
                                   tags$a(id = 'add_dataset',
                                          role='button',
                                          class = 'action-button shiny-bound-input',
                                          icon('plus', 'fa-fw'),
                                          'Import')
                           ),
                           tags$li(class=add,
                                   tags$a(id = 'integrate_dataset',
                                          role='button',
                                          class = 'action-button shiny-bound-input',
                                          tags$i(class= 'far fa-object-group fa-fw'),
                                          'Integrate')
                           ),
                           tags$li(role = 'separator', class='divider'),
                           tags$li(class=add,
                                   tags$a(id = 'export_dataset',
                                          role='button',
                                          class = 'action-button shiny-bound-input',
                                          tags$i(class = 'fab fa-r-project'),
                                          'Export')
                           ),
                           tags$li(class=add,
                                   tags$a(id = 'remove_dataset',
                                          role='button',
                                          class = 'action-button shiny-bound-input',
                                          tags$i(class= 'far fa-trash-alt fa-fw'),
                                          'Delete')
                           ),
                   )
          )
        )
      )
    )
  )


  return(ui)
}


#' UI for a tab pane
#'
#' @param tab The name of the tab
#' @param active The name of the active tab
#' @param ... The UI elements to place in the tab
#' @return shiny div tag with UI for tab
#'
#' @keywords internal
tabPane <- function(tab, active, ...) {
  active_class <- ifelse(tab == active, 'active', '')
  tags$div(class = paste('tab-pane', active_class), `data-value` = tab, id = id_from_tab(tab), ...)
}

#' Convert tab name to formated id
#'
#' used by navbarUI and *PageUI for dseqr app
#'
#' @param tab The name of the tab (e.g. \code{'Single Cell'})
#'
#' @keywords internal
id_from_tab <- function(tab) {
  id <- tolower(tab)
  gsub(' ', '-', id)
}


#' Full width button group with validation
#'
#' @param ... actionButtons
#' @param label Character vector. Label above button group.
#' @param container_id id of container. Used to toggle 'has-error' class with shinyjs::toggleClass.
#' @param help_id id of help block. Used to show help message in change shinyjs::html
#' @param class classes to add to btn-group divs.
#'
#' @keywords internal
justifiedButtonGroup <- function(..., label, container_id = NULL, help_id = NULL, class = '') {
  btns <- list(...)

  if (length(class) == 1) class <- rep(class, length(btns))

  tags$div(class = 'form-group selectize-fh', id = container_id,
           tags$label(class = 'control-label',  label),
           tags$div(class = 'btn-group btn-group-justified', role = 'group',
                    lapply(seq_along(btns), function(i) {
                      if (is.null(btns[[i]])) return(NULL)

                      tags$div(
                        class = paste('btn-group', class[i]),
                        role = 'group',
                        id = paste0(container_id, '-', i),
                        btns[[i]]
                      )
                    })
           ),
           span(id = help_id, class = 'help-block')
  )
}


# deprecated but keep because nice design
dsLabelRowsUI <- function(id) {
  ns <- NS(id)
  div(class = 'btn-group btn-group-justified',
      div(class = 'btn-group',
          tags$button(type = 'button',
                      class = 'btn btn-default dropdown-toggle',
                      `data-toggle`='dropdown',
                      `aria-haspopup`='true',
                      `aria-expanded`='false',
                      span('Label Selected Rows')
          ),
          tags$ul(class = 'dropdown-menu', style = 'width: 100%;',
                  tags$li(class="dropdown-header", 'Files'),
                  dropdownMenuButton(ns('pair'), 'Paired'),
                  dropdownMenuButton(ns('rep'), 'Replicates'),
                  tags$li(role = 'separator', class='divider'),
                  tags$li(class="dropdown-header", 'Group'),
                  dropdownMenuButton(ns('test'), 'Test'),
                  dropdownMenuButton(ns('ctrl'), 'Control'),
                  tags$li(role = 'separator', class='divider'),
                  dropdownMenuButton(ns('reset'), 'Reset Labels')
          )
      )
  )
}


#' Dropdown menu button for dsLabelRowsUI
#'
#' @keywords internal
#' @return \code{shiny.tag}
dropdownMenuButton <- function(id, label) {
  tags$li(
    tags$a(
      id = id, role = 'button', class = 'action-button shiny-bound-input', label
    )
  )
}
