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
        tags$div(id = ns('sc_clusters_container'), style = 'display: none;',
                 scSampleComparisonInput(ns)
        ),
        customQueryFormInput(ns('custom-query')),
        selectedDrugStudyInput(ns('drug_study')),
        advancedOptionsInput(ns('advanced'))

    )
  })
}

#' Input form for custom query on Drugs page
#' @export
#' @keywords internal
customQueryFormInput <- function(id) {
  ns <- NS(id)

  options = list(delimiter = ' ', create = I("function(input, callback){return {value: input,label: input};}"))

  tags$div(id = ns('custom_query_container'), class = 'hidden-form', style = 'display: none;',
           selectizeInput(ns('dn_genes'),
                          label = 'Genes to downregulate:',
                          choices = NULL,
                          multiple = TRUE,
                          width = '100%',
                          options = options),
           selectizeInput(ns('up_genes'),
                          label = 'Genes to upregulate:',
                          choices = NULL,
                          multiple = TRUE,
                          width = '100%',
                          options = options),
           textInputWithButtons(id = ns('custom_name'),
                                label = 'Name for custom query:',
                                actionButton(ns('submit_custom'), '', icon('plus', 'fa-fw')),
                                container_id = ns('validate'),
                                help_id = ns('error_msg'))
  )

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

  selectizeInputWithButtons(id = ns('query'), label = 'Select query signature:',
                            shiny::actionButton(ns('show_custom'), '', icon('object-group', 'fa-fw'), title = 'Toggle custom signature'),
                            options = list(optgroupField = 'type'))
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

  selectizeInputWithButtons(id = ns('study'), label = 'Select perturbation study:',
                            shiny::actionButton(ns('clinical'), label = '', icon = icon('pills', 'fa-fw'), onclick = 'toggleClinicalTitle(this)', title = 'only show compounds with a clinical phase'),
                            shiny::actionButton(ns('advanced'), label = '', icon = icon('cogs', 'fa-fw'), title = 'toggle advanced options'))

}

