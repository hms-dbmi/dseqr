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
               # docs section
               # li(class = 'navbar-right',
                  # a(href = paste0('#', id_from_tab('Docs')), `data-toggle` = 'tab', `data-value` = 'Docs', `aria-expanded` = 'false', 'Docs')
               # )
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
textInputWithButtons <- function(id, label, ..., container_id = NULL, help_id = NULL, tooltip = TRUE) {

  buttons <- list(...)
  buttons <- buttons[!sapply(buttons, is.null)]

  # generate tooltips
  button_tooltips <- NULL
  if (tooltip) {
    button_tooltips <- tags$div(
      lapply(buttons, function(btn)
        shinyBS::bsTooltip(id = btn$attribs$id, title = btn$attribs$title, options = list(container = 'body')))
    )
  }

  tags$div(class = 'form-group selectize-fh', id = container_id, class = 'validate-wrapper',
           tags$label(class = 'control-label', `for` = id, label),
           tags$div(class = 'input-group',
                    tags$input(id = id, type = 'text', class = 'form-control shiny-bound-input', value = '', placeholder = ''),
                    tags$span(class = 'input-group-btn',
                              lapply(buttons, function(btn) {
                                if (tooltip) btn$attribs$title <- NULL
                                return(btn)
                              }))
           ),
           tags$span(class = 'help-block', id = help_id),
           button_tooltips
  )
}

#' selectizeInput with validation
#' @inheritParams shiny::selectizeInput
#' @param container_id id of container. Used to toggle 'has-error' class with shinyjs::toggleClass.
#' @param help_id id of help block. Used to show help message in change shinyjs::html
#' @export
#' @keywords internal
selectizeInputWithValidation <- function(id, label, options = NULL, container_id = NULL, help_id = NULL, label_title = NULL, placement = 'top') {
  options <- ifelse(is.null(options), '{}', jsonlite::toJSON(options, auto_unbox = TRUE))

  label_tooltip <- NULL
  if (!is.null(label_title)) {
    label_id <- paste0(id, '-label-info')
    label_tooltip <- shinyBS::bsTooltip(label_id, title = label_title, placement = placement, options = list(container = 'body'))
    label <- tags$span(label, span(class='hover-info', span(id = label_id, icon('info', 'fa-fw'))))
  }

  tags$div(class = 'form-group selectize-fh', id = container_id,
           tags$label(class = 'control-label', `for` = id, label, title = label_title),
           tags$div(
             tags$select(id = id, style = 'display: none'),
             tags$script(type = 'application/json', `data-for` = id, HTML(options))
           ),
           tags$span(class = 'help-block', id = help_id),
           label_tooltip
  )
}


#' selectizeInput with buttons and validation
#' @inheritParams shiny::textInput
#' @param ... selectizeInput
#' @export
#' @keywords internal
selectizeInputWithButtons <- function(id,
                                      label,
                                      ...,
                                      options = NULL,
                                      container_id = NULL,
                                      help_id = NULL,
                                      label_title = NULL,
                                      tooltip = TRUE,
                                      placement = NULL,
                                      hide_btns = FALSE) {

  mult <- isTRUE(options$multiple)
  if(mult) {
    select_tag <- tags$select(id = id, style = 'display: none', multiple = TRUE)
  } else {
    select_tag <- tags$select(id = id, style = 'display: none')
  }

  buttons <- list(...)
  buttons <- buttons[!sapply(buttons, is.null)]

  if (length(buttons) == 0)
    return(selectizeInputWithValidation(id, label, options, container_id, help_id, label_title, placement = placement))

  # generate tooltips
  button_tooltips <- NULL
  if (tooltip) {
    button_tooltips <- tags$div(
      lapply(buttons, function(btn) {
        placement <- ifelse(is.null(placement),
                            ifelse(any(grepl('dropdown', unlist(btn$attribs))), 'right', 'bottom'),
                            placement)
        shinyBS::bsTooltip(id = btn$attribs$id, title = btn$attribs$title, placement = placement, options = list(container = 'body'))
      })
    )
  }


  # add info icon to label with tooltip
  label_tooltip <- NULL
  if (!is.null(label_title)) {
    label_id <- paste0(id, '-label-info')
    label_tooltip <- shinyBS::bsTooltip(label_id, title = label_title, placement = 'top', options = list(container = 'body'))
    label <- tags$span(label, span(class='hover-info', span(id = label_id, icon('info', 'fa-fw'))))
  }

  # if hide_btns allows full-width selectize
  ig_class  <- ifelse(hide_btns, '', 'input-group full-height-btn')
  fhs_class <- ifelse(hide_btns, '', 'full-height-selectize')

  options <- ifelse(is.null(options), '{}', jsonlite::toJSON(options, auto_unbox = TRUE))

  tags$div(class = 'form-group selectize-fh', id = container_id,
           tags$label(class = 'control-label', `for` = id, label),
           tags$div(class = ig_class, id = paste0(id, '-input-group'),
                    tags$div(class = fhs_class, id = paste0(id, '-full-height-selectize'),
                             select_tag,
                             tags$script(type = 'application/json', `data-for` = id, HTML(options))
                    ),
                    lapply(buttons, function(btn) {
                      is_dropdown <- any(grepl('dropdown', unlist(btn$attribs)))


                      if (hide_btns & is_dropdown)
                        btn$children[[1]]$attribs$style <- paste0('display: none;', btn$children[[1]]$attribs$style)

                      if (hide_btns & !is_dropdown)
                        btn$attribs$style <- paste0('display: none;', btn$attribs$style)

                      if (!is_dropdown)
                        btn <- tags$div(class = 'input-group-btn', id = paste0(btn$attribs$id, '-parent'), style = btn$attribs$`parent-style`, btn)

                      # remove title since using tooltips
                      if (tooltip) btn$attribs$title <- NULL

                      return(btn)
                    })
           ),
           tags$span(class = 'help-block', id = help_id),
           button_tooltips,
           label_tooltip
  )
}


#' Show/hide buttons in selectizeInputWithButtons
#'
#' @param selectize_id id of selectizeInputWithButtons element
#' @param button_ids character vector of ids for buttons in selectizeInputWithButtons
#' @param condition will show if TRUE and hide if FALSE
#' @export
toggleSelectizeButtons <- function(selectize_id, button_ids, condition) {
  # allows to take full-width/rounded corners
  toggleClass(paste0(selectize_id, '-input-group'),
              class = 'input-group full-height-btn',
              condition = condition)

  toggleClass(paste0(selectize_id, '-full-height-selectize'),
              class = 'full-height-selectize',
              condition = condition)

  #hide buttons
  for (button_id in button_ids)
    toggle(button_id, condition = condition)
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

inputLabelWithInfo <- function(id, label, tooltip = TRUE) {
  tagList(
    tags$span(id = id, label, span(class='hover-info', icon('info', 'fa-fw'))),
    shinyBS::bsTooltip(id = id, )

  )
}

