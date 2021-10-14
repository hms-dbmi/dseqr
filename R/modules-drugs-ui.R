#' UI for Drugs page
#'
#' @inheritParams scPageUI
#'
#' @return shiny.tag with html for drugs tab
#'
#' @export
drugsPageUI <- function(id, tab, active) {
  ns <- NS(id)

  withTags({
    tabPane(tab, active,
            rightClickMenu(),
            div(class = 'row',
                div(class = 'col-sm-12 col-md-5',
                    drugsFormInput(ns('form'))
                ),
                div(class = 'col-sm-12 col-md-7 mobile-margin-md',
                    drugsGenesPlotlyOutput(ns('genes'))
                )
            ),
            hr(),
            div(class = 'hide-btn-container',
                drugsTableOutput(ns('table')),
                shinyBS::bsTooltip(ns('show_genes'),
                                   'Show genes plot',
                                   options = list(
                                     container = 'body',
                                     template = '<div class="tooltip plot" role="tooltip"><div class="tooltip-arrow"></div><div class="tooltip-inner"></div></div>'
                                   )
                )
            )
    )

  })
}

#' Input form for Drugs page
#'
#' @keywords internal
drugsFormInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(class = "well-form well-bg",
        selectedAnalInput(ns('drugs')),
        customQueryFormInput(ns('custom-query')),
        div(id = ns('drug_study_container'), style = 'display: none',
            selectedDrugStudyInput(ns('drug_study'))
        ),
        div(id = ns('pert_signature_container'), style = 'display: none',
            selectedPertSignatureInput(ns('genes'))
        ),
        div(class = 'hidden-forms',
            advancedOptionsInput(ns('advanced'))
        )
    )

  })
}


#' UI for query/drug genes plotly
#'
#' @keywords internal
drugsGenesPlotlyOutput <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = 'scroll-plot', id = ns('container'), style = 'display: none',
             plotly::plotlyOutput(ns('plotly'), width = 'auto')
    )
  )
}


#' Output table for Drugs Page
#'
#' @keywords internal
drugsTableOutput <- function(id) {
  ns <- NS(id)

  tagList(
    hidden(downloadButton(ns('dl_drugs'), class = 'hide-btn', label = NULL)),
    tags$div(class = 'dt-container',
             DT::dataTableOutput(ns("query_table"))
    ),
    shinyBS::bsTooltip(ns('dl_drugs'),
                       'Download full query results',
                       options = list(
                         container = 'body',
                         template = '<div class="tooltip plot" role="tooltip"><div class="tooltip-arrow"></div><div class="tooltip-inner"></div></div>'
                       )
    )
  )
}


#' UI for seperate drugsPertInput for CMAP/L1000
#'
#' @keywords internal
selectedPertSignatureInput <- function(id) {
  ns <- NS(id)

  div(id = 'drugs-intro-pert',
      shinypanel::selectizeInputWithValidation(
        ns('pert'),
        label = 'Select perturbation for plot:',
        label_title = 'Perturbation signature (correlation)')

  )
}



#' Input form for custom query on Drugs page
#'
#' @keywords internal
customQueryFormInput <- function(id) {
  ns <- NS(id)

  options = list(delimiter = ' ', create = I("function(input, callback){return {value: input,label: input};}"))

  tags$div(id = ns('custom_query_container'), class = 'hidden-form', style = 'display: none;',
           shinypanel::textInputWithButtons(
             inputId = ns('custom_name'),
             label = 'Name for custom query:',
             actionButton(ns('click_custom'), '', icon('upload', 'fa-fw'), title = 'Upload csv with gene names and effect size'),
             container_id = ns('validate'),
             help_id = ns('error_msg')),

           # hidden upload button
           div(style = 'display: none',
               fileInput(ns('up_custom'), '', accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
           )
  )

}


#' advanced options input for drugs page
#'
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
        shinyWidgets::radioGroupButtons(ns('sort_by'), 'Sort based on correlation:', choices = c('average' = 'avg_cor', 'minimum' = 'min_cor'), justified = TRUE),
        shiny::sliderInput(ns('min_signatures'), 'Minimum number of signatures:',
                           min = 1,
                           max = 10,
                           value = 3,
                           ticks = TRUE,
                           width = '100%')
    )
  })
}


#' Select drugs study (CMAP or L1000) for drugs page
#'
#' @keywords internal
selectedDrugStudyInput <- function(id) {
  ns <- NS(id)

  div(id='drugs-intro-pert-study',

      shinypanel::selectizeInputWithButtons(
        inputId = ns('study'), label = 'Select perturbation study:',
        shiny::actionButton(ns('direction'), label = '', icon = icon('arrows-alt-v', 'fa-fw'), title = 'change direction of correlation', `parent-style` = 'display: none;'),
        shiny::actionButton(ns('clinical'), label = '', icon = icon('pills', 'fa-fw'), title = 'only show compounds with a clinical phase'),
        shiny::actionButton(ns('advanced'), label = '', icon = icon('cogs', 'fa-fw'), title = 'toggle advanced options'))
  )


}


#' Custom right click menu for selecting correlation point as query
#'
#' @keywords internal
rightClickMenu <- function() {
  tags$ul(
    class = 'custom-menu', id = 'cor-menu',
    tags$li(id='menu-title', 'Signature title'),
    tags$li('data-action' = 'load', id = 'cor-signature', 'Load as query')
  )
}


#' Analysis input for Drugs page
#'
#' @keywords internal
selectedAnalInput <- function(id, label = 'Select a dataset or query signature:') {
  ns <- NS(id)

  tagList(
    div(id='drugs-intro-query',
        shinypanel::selectizeInputWithButtons(
          inputId = ns('query'), label = label,
          shiny::actionButton(ns('show_custom'), '', icon = tags$i(class ='far fa-fw fa-file-alt'), title = 'Toggle custom signature'),
          options = list(optgroupField = 'type'))
    ),
    tags$div(id='drugs-intro-comparison',
             tags$div(id = ns('sc_clusters_container'), style = 'display: none;',
                      scSampleGroupsInput(ns('sample_groups')),
                      scSampleClustersInput(ns('sample_clusters'))
             ),
             tags$div(id = ns('bulk_groups_container'), style = 'display: none;',
                      bulkAnalInput(ns('bulk'), with_dl = FALSE)
             )

    )
  )

}
