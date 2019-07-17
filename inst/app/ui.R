
dsPageUI <- function(id, tab, active) {
  ns <- NS(id)
  withTags({
    tabPane(tab, active,
            div(class = 'row',
                div(class = 'col-sm-6',
                    dsFormInput(ns('form'))
                )
            ),
            hr(),
            div(id = ns('quant_table_container'), style = 'display: none;',
                dsTable(ns('quant'))
            ),
            div(id = ns('anal_table_container'), style = 'display: none;',
                dsTable(ns('anal'))
            )
    )
  })
}

dsFormInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(class = "well-form well-bg",
        dsDatasetInput(ns('selected_dataset')),
        div(id = ns('quant_dataset_panel'), style = 'display: none;',
            dsFormQuantInput(ns('quant_form'))
        ),
        div(id = ns('anal_dataset_panel'), style = 'display: none;',
            dsFormAnalInput(ns('anal_form'))
        )
    )
  })
}

dsFormQuantInput <- function(id) {
  ns <- NS(id)

  tagList(
    dsEndTypeInput(ns('end_type')),
    justifiedButtonGroup(
      id = ns('quant_labels'),
      label = 'Label selected rows:',
      help_block = span(id = ns('error_msg'), class = 'help-block'),
      actionButton(ns('pair'), 'Paired'),
      actionButton(ns('rep'), 'Replicate'),
      actionButton(ns('reset'), 'Reset')
    ),
    actionButton(ns('run_quant'), 'Run Quantification', width = '100%', class = 'btn-primary')
  )
}

dsFormAnalInput <- function(id) {
  ns <- NS(id)

  tagList(
    selectizeInputWithValidation(
      ns('anal_name'),
      'Analysis name:',
      options = list(create = TRUE, placeholder = 'Type name to create new analysis'),
      container_id = ns('anal_name_container'),
      help_id = ns('anal_name_help')
    ),
    div(id = ns('anal_buttons_panel'),
      justifiedButtonGroup(
        id = ns('anal_labels'),
        label = 'Label selected rows:',
        help_block = span(id = ns('labels_help'), class = 'help-block'),
        actionButton(ns('test'), 'Test'),
        actionButton(ns('ctrl'), 'Control'),
        actionButton(ns('reset'), 'Reset')
      ),
      actionButton(ns('run_anal'), 'Run Analysis', width = '100%', class = 'btn-primary')
    )
  )
}

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



dsDatasetInput <- function(id) {
  ns <- NS(id)

  withTags({
    selectizeInputWithButtons(ns('dataset_name'), 'Dataset name:', options = list(create = TRUE, placeholder = 'Type name to add quant dataset'),
                              button(id = ns('dataset_dir'), type = 'button', class="shinyDirectories btn btn-default action-button shiny-bound-input disabled",
                                     `data-title` = 'Folder with fastq.gz files',
                                     title = 'Select folder with fastq.gz files',
                                     i(class = 'far fa-folder')

                              )
    )
  })
}

dsEndTypeInput <- function(id) {
  ns <- NS(id)

  selectizeInput(ns('end_type'),
                 'Confirm end-type:',
                 choices = NULL, width = '100%')
}

dsTable <- function(id) {
  ns <- NS(id)
  DT::dataTableOutput(ns("pdata"))
}

# general utility components
textInputWithButtons <- function(id, label, ...) {
  tags$div(class = 'form-group selectize-fh',
           tags$label(class = 'control-label', `for` = id, label),
           tags$div(class = 'input-group',
                    tags$input(id = id, type = 'text', class = 'form-control shiny-bound-input', value = '', placeholder = ''),
                    tags$span(class = 'input-group-btn', ...)
           )
  )
}

textInput <- function(id, label, container_id = NULL, help_id = NULL) {
  tags$div(class = 'form-group selectize-fh', id = container_id,
           tags$label(class = 'control-label', `for` = id, label),
           tags$input(id = id, type = 'text', class = 'form-control shiny-bound-input', value = '', placeholder = ''),
           tags$span(class='help-block', id = help_id)
  )
}

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

dropdownMenuButton <- function(id, label) {
  tags$li(
    tags$a(
      id = id, role = 'button', class = 'action-button shiny-bound-input', label
    )
  )
}

justifiedButtonGroup <- function(..., label, id = NULL, help_block = NULL) {
  tags$div(class = 'form-group selectize-fh', id = id,
           tags$label(class = 'control-label',  label),
           tags$div(class = 'btn-group btn-group-justified', role = 'group',
                    lapply(list(...), function(btn) {
                      tags$div(class = 'btn-group', role = 'group', btn)
                    })
           ),
           help_block
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


#' UI for Drugs page
#' @export
#' @keywords internal
drugsPageUI <- function(id, tab, active) {
  ns <- NS(id)

  withTags({
    tabPane(tab, active,
            div(class = 'row',
                div(class = 'col-sm-6',
                    drugsFormInput(ns('form'))
                )
            ),
            hr(),
            drugsTableOutput(ns('table'))
    )
  })
}

#' Output table for Drugs Page
#' @export
#' @keywords internal
drugsTableOutput <- function(id) {
  ns <- NS(id)

  withTags({
    div(class = 'dt-container',
        DT::dataTableOutput(ns("query_table"))
    )
  })
}

#' Input form for Drugs page
#' @export
#' @keywords internal
drugsFormInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(class = "well-form well-bg",
        querySignatureInput(ns('signature')),
        selectedDrugStudyInput(ns('drug_study')),
        advancedOptionsInput(ns('advanced'))

    )
  })
}

#' Query signature input for Drugs page
#' @export
#' @keywords internal
querySignatureInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(class = 'form-group selectize-fh',
        label(class = 'control-label', `for` = ns('query'), 'Select query signature:'),
        div(
          select(id = ns('query'), style = 'display: none'),
          script(type = 'application/json', `data-for` = ns('query'), HTML('{"optgroupField": "dataset_name"}'))
        )
    )
  })
}

#' advanced options input for drugs page
#' @export
#' @keywords internal
advancedOptionsInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(id = ns('advanced-panel'), class = 'hidden-form', style = 'display: none;',
        selectizeInput(ns('cells'), 'Select cell lines:', choices = NULL, multiple = TRUE, options = list(placeholder = "showing all"), width = '100%'),
        shinyWidgets::radioGroupButtons(ns('sort_by'), 'Sort based on correlation:', choices = c('minimum' = 'min_cor', 'average' = 'avg_cor'), justified = TRUE)
    )
  })
}

#' Select drugs study (CMAP or L1000) for drugs page
#' @export
#' @keywords internal
selectedDrugStudyInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(class = 'form-group selectize-fh',
        label(class = 'control-label', `for` = ns('study'), 'Select drug study:'),
        div(class = 'input-group',
            div(
              select(id = ns('study'), style = 'display: none'),
              script(type = 'application/json', `data-for` = ns('study'), HTML('{}'))
            ),
            div(class = 'input-group-btn',
                shinyBS::bsButton(ns('clinical'), label = '', icon = icon('pills'), style = 'default', onclick = 'toggleClinicalTitle(this)', title = 'only show compounds with a clinical phase'),
                shinyBS::bsButton(ns('advanced'), label = '', icon = icon('cogs'), style = 'default', title = 'toggle advanced options')
            )
        )
    )
  })
}


#-----
tabs <- c('Datasets', 'Single Cell', 'Drugs')
active <- 'Datasets'

bootstrapPage(
  useShinyjs(),
  includeScript(path = 'www/contrasts.js'),
  includeScript(path = 'www/toggleClinicalTitle.js'),
  includeCSS(path = 'www/custom.css'),
  includeCSS(path = 'www/drugs.css'),
  navbarUI(tabs, active),
  fluidPage(
    # make sure css/js loaded from packages where not using functions (not using default)
    tags$div(style = 'display: none;', selectizeInput('blah1', label = NULL, choices = '')),
    tags$div(style = 'display: none;', shinyFiles::shinyDirButton("blah2", title='', label='', icon=icon('plus'))),

    tags$div(class = "tab-content", `data-tabsetid` = "tabset", id = "tabs",
             # single cell tab
             scPageUI("sc", tab = 'Single Cell', active),
             drugsPageUI("drug", tab = 'Drugs', active),
             dsPageUI('datasets', tab = 'Datasets', active)
    )
  )
)
