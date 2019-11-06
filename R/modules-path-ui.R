#' UI for pathways page
#' @export
#' @keywords internal
pathPageUI <- function(id, tab, active) {
  ns <- NS(id)
  withTags({
    tabPane(tab, active,
            div(class = 'row',
                div(class = 'col-sm-6',
                    pathFormInput(ns('form'))
                )
            ),
            hr(),
            div(id = ns('pert-legend'), class = 'pert-label',
                label(class='control-label', 'Perturbation direction:'),
                div(span(class = 'input-swatch opposite'), 'opposite'),
                div(span(class = 'input-swatch same'), 'same')
            ),
            div(class = 'scroll-plot',
                plotly::plotlyOutput(ns('path_plot'), height = '550px')

            )
    )
  })
}



#' Input form for pathways page
#' @export
#' @keywords internal
pathFormInput <- function(id) {
  ns <- NS(id)

  tags$div(class = "well-form well-bg",
           selectizeInputWithValidation(ns('anal'), 'Select an analysis:', options = list(optgroupField = 'type')),
           tags$div(id = ns('sc_clusters_container'), style = 'display: none;',
                    scSampleComparisonInput(ns)
           ),
           selectizeInputWithButtons(ns('pathway'),
                                     label = tags$span('Select a pathway:', span(class='hover-info', icon('info', 'fa-fw'))),
                                     actionButton(ns('kegg'), '', icon = icon('external-link-alt', 'fa-fw'), title = 'Go to KEGG'),
                                     options = list(optgroupField = 'direction_label', searchField = c('text', 'optgroup')),
                                     label_title = 'Pathway (FDR)'),
           tags$div(id = ns('l1000_pert_container'), style = 'display: none;',
                    pathPertInput(ns('l1000'), type = 'L1000')
           ),
           tags$div(id = ns('cmap2_pert_container'), style = 'display: none;',
                    pathPertInput(ns('cmap2'), type = 'CMAP02')
           )
           # selectizeInput(ns('custom_path_genes'), 'Custom gene set:', choices = NULL, multiple = TRUE, options = list(render = I('{option: pathGene, item: pathGene}')), width = '100%')

  )
}


#' Input for single cell clusters in pathFormInput
#' @export
#' @keywords internal
scSampleComparisonInput <- function(ns) {

  button <- actionButton(ns('run_comparison'), '',
                         icon = icon('chevron-right', 'far fa-fw'),
                         title = 'Compare test to control cells')

  selectizeInputWithButtons(ns('selected_clusters'),
                            label = tags$span('Compare samples for:', span(class='hover-info', icon('info', 'fa-fw'))),
                            button,
                            options = list(multiple = TRUE),
                            label_title = 'Cluster (n test :: n ctrl)')
}

#' Input for perturbation selection for Pathways page
#' @export
#' @keywords internal
pathPertInput <- function(id, type) {
  ns <- NS(id)

  if (type == 'CMAP02') {
    res <- selectizeInputWithValidation(ns('pert'),
                                        label = tags$span('Select perturbation signature:', span(class='hover-info', icon('info', 'fa-fw'))),
                                        label_title = 'Signature (correlation)')

  } else if (type == 'L1000') {
    res <- selectizeInputWithButtons(ns('pert'),
                                     label = tags$span('Select perturbation signature:', span(class='hover-info', icon('info', 'fa-fw'))),
                                     actionButton(ns('l1000_type'), '', icon = icon('pills', 'fa-fw'), title = 'Toggle between drug and genetic perts'),
                                     label_title = 'Signature (correlation)'
    )
  }
  return(res)
}
