#' UI for navbar
#' @param active the active tab name
#'
#' @export
navbarUI <- function(tabs, active) {

  withTags({
    tags$nav(class = 'navbar navbar-default navbar-static-top',
             div(class = 'container-fluid',
                 div(class = 'navbar-header',
                     tags$button(type='button', class='navbar-toggle collapsed', `data-toggle`='collapse', `data-target`='#bs-navbar', `aria-expanded`='false',
                                 span(class = 'sr-only', 'Toggle navigation'),
                                 span(class = 'icon-bar'),
                                 span(class = 'icon-bar'),
                                 span(class = 'icon-bar')
                     ),
                     span(class = 'navbar-brand', title = 'drugseqr',
                          span(class = 'brand-icons',
                               tags$i(class = 'glyphicon glyphicon-leaf'),
                               'drugseqr'
                          )
                     )
                 ),
                 div(id = 'bs-navbar', class = 'collapse navbar-collapse',
                     tags$ul(class = 'nav navbar-nav shiny-tab-input shiny-bound-input', `data-tabsetid` = 'tabset', id = 'tabs',
                             lapply(seq_along(tabs), function(i) {
                               tab <- tabs[i]
                               is.active <- tab == active
                               tags$li(class = ifelse(is.active, 'active', ''), `data-toggle`="collapse", `data-target`=".navbar-collapse.in",
                                       a(href = paste0('#', id_from_tab(tab)), `data-toggle` = 'tab', `data-value` = tab, `aria-expanded` = ifelse(is.active, 'true', 'false'), tab)
                               )
                             }),
                             # TODO: uncomment once repo is public
                             # github linkout section
                             # li(class = 'navbar-right',
                             #    a(href = 'https://github.com/hms-dbmi/drugseqr', icon('github'), style = 'padding-bottom: 0px; font-size:17px;',)
                             # ),
                             # docs section
                             tags$li(class = 'navbar-right',
                                     a(href = paste0('#', id_from_tab('Docs')), `data-toggle` = 'tab', `data-value` = 'Docs', `aria-expanded` = 'false', 'Docs')
                             )
                     )
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
#'
#' @keywords internal
tabPane <- function(tab, active, ...) {
  active_class <- ifelse(tab == active, 'active', '')
  tags$div(class = paste('tab-pane', active_class), `data-value` = tab, id = id_from_tab(tab), ...)
}

#' Convert tab name to formated id
#'
#' used by navbarUI and *PageUI for drugseqr app
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
#'
#' @keywords internal
justifiedButtonGroup <- function(..., label, container_id = NULL, help_id = NULL) {
  tags$div(class = 'form-group selectize-fh', id = container_id,
           tags$label(class = 'control-label',  label),
           tags$div(class = 'btn-group btn-group-justified', role = 'group',
                    lapply(list(...), function(btn) {
                      tags$div(class = 'btn-group', role = 'group', btn)
                    })
           ),
           span(id = help_id, class = 'help-block')
  )
}


# deprecated but keep because nice design
dsLabelRowsUI <- function(id) {
  ns <- NS(id)
  withTags({
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
  })
}


#' Dropdown menu button for dsLabelRowsUI
#'
#' @keywords internal
dropdownMenuButton <- function(id, label) {
  tags$li(
    tags$a(
      id = id, role = 'button', class = 'action-button shiny-bound-input', label
    )
  )
}


