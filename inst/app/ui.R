# UI for Drugs page
drugPageUI <- function(id, tab, active) {
  ns <- NS(id)
  active_class <- ifelse(tab == active, 'active', '')

  withTags({
    div(class = paste('tab-pane', active_class), `data-value` = tab, id = id_from_tab(tab),
        div(class = 'row',
            div(class = 'col-sm-6',
                drugFormInput(ns('form'))
            )
        ),
        hr(),
        drugTableOutput(ns('table'))
    )
  })
}

# Output table for Drugs Page
drugTableOutput <- function(id) {
  ns <- NS(id)

  withTags({
    div(class = 'dt-container',
        DT::dataTableOutput(ns("query_res"))
    )
  })
}

#' Input form for Drugs page
drugFormInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(class = "well-form well-bg",
        querySignatureInput(ns('signature')),
        selectedDrugStudyInput(ns('drug_study')),
        advancedOptionsInput(ns('advanced'))

    )
  })
}

querySignatureInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(class = 'form-group selectize-fh',
        label(class = 'control-label', `for` = ns('query'), 'Select query signature:'),
        div(class = 'input-group',
            div(
              select(id = ns('query'), style = 'display: none'),
              script(type = 'application/json', `data-for` = ns('query'), HTML('{}'))
            ),
            div(class = 'input-group-btn',
                shinyBS::bsButton(ns('run_query'), label = '', icon = icon('search'), title = 'Run query')
            )
        )
    )
  })
}

advancedOptionsInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(id = ns('advanced-panel'), class = 'hidden-form', style = 'display: none;',
        selectizeInput(ns('cells'), 'Select cell lines:', choices = NULL, multiple = TRUE, options = list(placeholder = "showing all"), width = '100%')
    )
  })
}

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


tabs <- c('Datasets', 'Single Cell', 'Bulk', 'Pathways', 'Drugs')
active <- 'Drugs'

bootstrapPage(
  useShinyjs(),
  includeScript(path = 'www/contrasts.js'),
  includeScript(path = 'www/toggleClinicalTitle.js'),
  includeCSS(path = 'www/custom.css'),
  includeCSS(path = 'www/drugs.css'),
  navbarUI(tabs, active),
  fluidPage(
    # make sure selectize loaded (not using default)
    tags$div(style = 'display: none', selectizeInput('blah1', label = NULL, choices = '')),

    # THE TABS ----
    tags$div(class = "tab-content", `data-tabsetid` = "tabset", id = "tabs",
             # single cell tab
             scPageUI("sc", tab = 'Single Cell', active),
             bulkPageUI("bulk", tab = 'Bulk', active),
             drugPageUI("drug", tab = 'Drugs', active)
    )
  )
)
