#' UI for Drugs page
#' @export
#' @keywords internal
drugsPageUI <- function(id, tab, active) {
  ns <- NS(id)

  withTags({
    tabPane(tab, active,
            rightClickMenu(),
            div(class = 'row',
                div(class = 'col-sm-6',
                    drugsFormInput(ns('form'))
                )
            ),
            hr(),
            div(class = 'hide-btn-container',
                drugsGenesPlotlyOutput(ns('genes')),
                downloadButton(ns('dl_drugs'), class = 'hide-btn', label = NULL),
                drugsTableOutput(ns('table')),
                shinyBS::bsTooltip(ns('dl_drugs'),
                                   'Download full query results',
                                   options = list(
                                     container = 'body',
                                     template = '<div class="tooltip ggplot" role="tooltip"><div class="tooltip-arrow"></div><div class="tooltip-inner"></div></div>'
                                   )
                ),
                shinyBS::bsTooltip(ns('show_genes'),
                                   'Show genes plot',
                                   options = list(
                                     container = 'body',
                                     template = '<div class="tooltip ggplot" role="tooltip"><div class="tooltip-arrow"></div><div class="tooltip-inner"></div></div>'
                                   )
                )
            )

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
        selectedAnalInput(ns('anal')),
        customQueryFormInput(ns('custom-query')),
        selectedDrugStudyInput(ns('drug_study')),
        div(class = 'hidden-forms',
            advancedOptionsInput(ns('advanced')),
            selectedPertSignatureInput(ns('genes'))
        )

    )
  })
}

#' Analysis input for Drugs and Pathways page
#' @export
#' @keywords internal
selectedAnalInput <- function(id) {
  ns <- NS(id)

  tagList(
    selectizeInputWithButtons(id = ns('query'), label = 'Select dataset or query signature:',
                              shiny::actionButton(ns('show_custom'), '', icon('object-group', 'fa-fw'), title = 'Toggle custom signature'),
                              options = list(optgroupField = 'type')),
    tags$div(id = ns('sc_clusters_container'), style = 'display: none;',
             scSampleComparisonInput(ns)
    ),
    tags$div(id = ns('bulk_groups_container'), style = 'display: none;',
             bulkAnalInput(ns('drugs'), with_dl = FALSE)
    )
  )

}


#' UI for query/drug genes plotly
#' @export
#' @keywords internal
drugsGenesPlotlyOutput <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(style='display: none;', id = ns('container'),
             tags$div(class = 'scroll-plot',
                      plotly::plotlyOutput(ns('plotly'), height = '550px')
             ),
             hr()
    )
  )
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


#' UI for seperate drugsPertInput for CMAP/L1000
#' @export
#' @keywords internal
selectedPertSignatureInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(id = ns('pert_container'),  class = 'hidden-form bottom', style = 'display: none;',
        drugsPertInput(ns('pert'))
    )
  })
}

#' UI to select perturbation for Drugs page
#' @export
#' @keywords internal
drugsPertInput <- function(id) {

  ns <- NS(id)
  selectizeInputWithValidation(ns('pert'),
                               label = 'Select perturbation for plot:',
                               label_title = 'Perturbation signature (correlation)')
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

  selectizeInputWithButtons(ns('selected_clusters'),
                            label = 'Compare samples for:',
                            button,
                            options = list(multiple = TRUE))

}

#' advanced options input for drugs page
#' @export
#' @keywords internal
advancedOptionsInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(id = ns('advanced-panel'), class = 'hidden-form bottom', style = 'display: none;',
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
                            shiny::actionButton(ns('direction'), label = '', icon = icon('arrows-alt-v', 'fa-fw'), title = 'change direction of correlation', `parent-style` = 'display: none;'),
                            shiny::actionButton(ns('clinical'), label = '', icon = icon('pills', 'fa-fw'), onclick = 'toggleClinicalTitle(this)', title = 'only show compounds with a clinical phase'),
                            shiny::actionButton(ns('advanced'), label = '', icon = icon('cogs', 'fa-fw'), title = 'toggle advanced options'),
                            shiny::actionButton(ns('show_genes'), label = '', icon = icon('chart-line', 'fa-fw'), title = 'show plot of query genes'))

}


#' Custom right click menu for selecting correlation point as query
#' @export
#' @keywords internal
rightClickMenu <- function() {
  withTags({
    ul(
      class = 'custom-menu', id = 'cor-menu',
      li(id='menu-title', 'Signature title'),
      li('data-action' = 'load', id = 'cor-signature', 'Load as query')
    )
  })
}


tabs <- c('Bulk Data', 'Single Cell', 'Pathways', 'Drugs')
active <- 'Drugs'

bootstrapPage(
  useShinyjs(),
  # scrollspy for docs tab
  extendShinyjs(text = "shinyjs.init = function() {$('body').scrollspy({ target: '.bs-docs-sidenav', offset: 60 });}"),
  includeScript(path = 'www/renderSelectize.js'),
  includeScript(path = 'www/toggleClinicalTitle.js'),
  includeScript(path = 'www/contextMenu.js'),
  includeScript(path = 'www/anchor-polyfill.js'),
  includeCSS(path = 'www/custom.css'),
  includeCSS(path = 'www/bs-docs.css'),
  includeCSS(path = 'www/drugs.css'),
  includeCSS(path = 'www/pathways.css'),
  navbarUI(tabs, active),
  fluidPage(
    # make sure css/js loaded from packages where not using functions (not using default)
    tags$div(style = 'display: none;', selectizeInput('blah1', label = NULL, choices = '')),
    tags$div(style = 'display: none;', shinyFiles::shinyDirButton("blah2", title='', label='', icon=icon('plus'))),

    tags$div(class = "tab-content", `data-tabsetid` = "tabset", id = "tabs",
             scPageUI("sc", tab = 'Single Cell', active),
             drugsPageUI("drug", tab = 'Drugs', active),
             dsPageUI('bulk', tab = 'Bulk Data', active),
             pathPageUI('pathways', tab = 'Pathways', active),
             docsPageUI('docs', tab = 'docs', active)

    )
  )
)



