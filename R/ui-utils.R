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


#' selectizeInput with validation
#' @inheritParams shiny::selectizeInput
#' @param container_id id of container. Used to toggle 'has-error' class with shinyjs::toggleClass.
#' @param help_id id of help block. Used to show help message in change shinyjs::html
#' @export
#' @keywords internal
selectizeInputWithValidation <- function(id, label, options = NULL, container_id = NULL, help_id = NULL) {
  options <- ifelse(is.null(options), '{}', jsonlite::toJSON(options, auto_unbox = TRUE))

  tags$div(class = 'form-group selectize-fh', id = container_id,
           tags$label(class = 'control-label', `for` = id, label),
           tags$div(
             tags$select(id = id, style = 'display: none'),
             tags$script(type = 'application/json', `data-for` = id, HTML(options))
           ),
           tags$span(class = 'help-block', id = help_id)
  )
}

#' textInput with validation
#' @inheritParams shiny::textInput
#' @param container_id id of container. Used to toggle 'has-error' class with shinyjs::toggleClass.
#' @param help_id id of help block. Used to show help message in change shinyjs::html
#' @export
#' @keywords internal
textInputWithValidation <- function(id, label, container_id = NULL, help_id = NULL) {
  tags$div(class = 'form-group selectize-fh', id = container_id,
           tags$label(class = 'control-label', `for` = id, label),
           tags$input(id = id, type = 'text', class = 'form-control shiny-bound-input', value = '', placeholder = ''),
           tags$span(class='help-block', id = help_id)
  )
}

#' textInput with buttons
#' @inheritParams shiny::textInput
#' @param ... actionButtons
#' @export
#' @keywords internal
textInputWithButtons <- function(id, label, ...) {
  tags$div(class = 'form-group selectize-fh',
           tags$label(class = 'control-label', `for` = id, label),
           tags$div(class = 'input-group',
                    tags$input(id = id, type = 'text', class = 'form-control shiny-bound-input', value = '', placeholder = ''),
                    tags$span(class = 'input-group-btn', ...)
           )
  )
}

#' selectizeInput with buttons
#' @inheritParams shiny::textInput
#' @param ... selectizeInput
#' @export
#' @keywords internal
selectizeInputWithButtons <- function(id, label, options = NULL, ...) {

  options <- ifelse(is.null(options), '{}', jsonlite::toJSON(options, auto_unbox = TRUE))

  tags$div(class = 'form-group selectize-fh',
           tags$label(class = 'control-label', `for` = id, label),
           tags$div(class = 'input-group',
                    tags$div(
                      tags$select(id = id, style = 'display: none'),
                      tags$script(type = 'application/json', `data-for` = id, HTML(options))
                    ),
                    tags$div(class = 'input-group-btn',
                             # the buttons
                             ...
                    )
           )
  )
}


#' Full width button group with validation
#'
#' @param ... actionButtons
#' @param label Character vector. Label above button group.
#' @param container_id id of container. Used to toggle 'has-error' class with shinyjs::toggleClass.
#' @param help_id id of help block. Used to show help message in change shinyjs::html
#' @export
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
            button(type = 'button',
                   class = 'btn btn-default dropdown-toggle',
                   `data-toggle`='dropdown',
                   `aria-haspopup`='true',
                   `aria-expanded`='false',
                   span('Label Selected Rows')
            ),
            ul(class = 'dropdown-menu', style = 'width: 100%;',
               li(class="dropdown-header", 'Files'),
               dropdownMenuButton(ns('pair'), 'Paired'),
               dropdownMenuButton(ns('rep'), 'Replicates'),
               li(role = 'separator', class='divider'),
               li(class="dropdown-header", 'Group'),
               dropdownMenuButton(ns('test'), 'Test'),
               dropdownMenuButton(ns('ctrl'), 'Control'),
               li(role = 'separator', class='divider'),
               dropdownMenuButton(ns('reset'), 'Reset Labels')
            )
        )
    )
  })
}


#' Dropdown menu button for dsLabelRowsUI
#' @export
#' @keywords internal
dropdownMenuButton <- function(id, label) {
  tags$li(
    tags$a(
      id = id, role = 'button', class = 'action-button shiny-bound-input', label
    )
  )
}

