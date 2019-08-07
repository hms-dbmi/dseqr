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
        tags$div(id = ns('sc_clusters_container'), style = '',
                 scSampleComparisonInput(ns)
        ),
        selectedDrugStudyInput(ns('drug_study')),
        advancedOptionsInput(ns('advanced'))

    )
  })
}

#' Input for Single Cell sample comparison
#'
#' Used in both Drugs and Pathways tab
#'
#' @export
#' @keywords internal
scSampleComparisonInput <- function(ns) {

  button <- actionButton(ns('run_comparison'), '',
                         icon = icon('chevron-right', 'far fa-fw'),
                         title = 'Compare test to control cells')

  selectizeInputMultWithButton(ns('selected_clusters'), label = 'Compare samples for:', button = button)

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
          script(type = 'application/json', `data-for` = ns('query'), HTML('{"optgroupField": "type"}'))
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
        selectizeInput(ns('cells'),
                       'Select cell lines:',
                       choices = NULL,
                       multiple = TRUE,
                       options = list(placeholder = "showing all",
                                      render = I('{option: cellOptions}'),
                                      optgroupField = 'primary_site'),
                       width = '100%'),
        shinyWidgets::radioGroupButtons(ns('sort_by'), 'Sort based on correlation:', choices = c('minimum' = 'min_cor', 'average' = 'avg_cor'), justified = TRUE),
        shiny::sliderInput(ns('min_signatures'), 'Minimum number of signatures:',
                           min = 1,
                           max = 10,
                           value = 1,
                           ticks = TRUE,
                           width = '100%')
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

