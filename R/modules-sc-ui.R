
#' UI for Single Cell Exploration page
#' @export
#' @keywords internal
scPageUI <- function(id, tab, active) {
  ns <- NS(id)
  active_class <- ifelse(tab == active, 'active', '')
  withTags({
    div(class = paste('tab-pane', active_class), `data-value` = tab, id = id_from_tab(tab),
        div(class = 'row',
            div(class = 'col-sm-6',
                scFormInput(ns('form'))
            ),
            div(class = 'col-sm-6',
                scClusterPlotOutput(ns('cluster_plot'))
            )
        ),
        hr(),
        # row for cluster comparison
        div(class = 'row', id = ns('cluster_comparison_row'),
            div(class = "col-sm-6 col-lg-6 col-lg-push-6",
                scMarkerPlotOutput(ns('marker_plot_cluster'))
            ),
            div(class = "col-sm-6 col-lg-6 col-lg-pull-6",
                scBioGpsPlotOutput(ns('biogps_plot'))
            )
        ),
        # row for sample (test vs ctrl) comparison
        div(class = 'row', id = ns('sample_comparison_row'), style = 'display: none;',
            div(class = "col-sm-6 col-lg-6 col-lg-push-6",
                scMarkerPlotOutput(ns('marker_plot_ctrl'))
            ),
            div(class = "col-sm-6 col-lg-6 col-lg-pull-6",
                scMarkerPlotOutput(ns('marker_plot_test'))
            )
        )
    )
  })
}


#' Input form for Single Cell Exploration page
#' @export
#' @keywords internal
scFormInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(class = "well-form well-bg",
        selectedAnalInput(ns('anal')),
        integrationFormInput(ns('integration')),
        comparisonTypeToggle(ns('comparison')),
        # inputs for comparing clusters
        div(id = ns('cluster_comparison_inputs'),
            clusterComparisonInput(ns('cluster')),
            selectedGeneInput(ns('gene_clusters'))
        ),
        # inputs for comparing samples
        div(id = ns('sample_comparison_inputs'), style = 'display: none',
            sampleComparisonInput(ns('sample')),
            selectedGeneInput(ns('gene_samples'), sample_comparison = TRUE)
        )
    )
  })
}



#' Selection form/button for sample comparisons (test vs control)
#' @export
#' @keywords internal
sampleComparisonInput <- function(id) {
  ns <- NS(id)

  button <- actionButton(ns('run_comparison'), '',
                         icon = icon('chevron-right', 'far fa-fw'),
                         title = 'Compare test to control cells')

  selectizeInputWithButtons(ns('selected_clusters'),
                            label = tags$span('Compare samples for:', span(class='hover-info', icon('info', 'fa-fw'))),
                            button,
                            options = list(multiple = TRUE),
                            label_title = 'Cluster (n test :: n ctrl)')
}


#' Input form to control/test/all groups for integrated datasets
#' @export
#' @keywords internal
comparisonTypeToggle <- function(id) {
  ns <- NS(id)

  withTags({
    div(style = 'display: none;', id = ns('comparison_type_container'), class = 'selectize-fh form-group',
        shinyWidgets::radioGroupButtons(ns('comparison_type'), "Perform comparisons between:",
                                        choices = c('clusters', 'samples'),
                                        selected = 'clusters', justified = TRUE)

    )

  })

}

#' Input form/associated buttons for selecting single cell analysis
#' @export
#' @keywords internal
selectedAnalInput <- function(id) {
  ns <- NS(id)

  selectizeInputWithButtons(ns('selected_anal'), 'Select a dataset:',
                            showLabelTransferButton(ns('label-transfer')),
                            showIntegrationButton(ns('integration')),
                            plotStylesButton(ns('styles'))
  )
}

#' Button with to toggle display of label transfer inputs
#' @export
#' @keywords internal
showLabelTransferButton <- function(id) {
  ns <- NS(id)

  actionButton(ns('show_label_transfer'), '',
               icon = icon('tag', 'fa-fw'),
               title = 'Toggle label transfer', class = 'squashed-btn')
}



#' Button with sliders for adjusting plot jitter and point size
#' @export
#' @keywords internal
plotStylesButton <- function(id) {
  ns <- NS(id)
  dropdownButton(
    sliderInput(ns('point_size'), 'Point size:',
                width = '100%', ticks = FALSE,
                min = 0.5, max = 4, value = 2.5, step = 0.5),
    sliderInput(ns('point_jitter'), 'Point jitter:',
                width = '100%', ticks = FALSE,
                min = 0, max = 3, value = 0, step = 0.5),
    circle = FALSE, right = TRUE, icon = icon('cog', 'fa-fw'),
    title = 'Show plot style controls'
  )
}



#' Button with to toggle display of integrationFormInput
#' @export
#' @keywords internal
showIntegrationButton <- function(id) {
  ns <- NS(id)

  actionButton(ns('show_integration'), '',
               icon = icon('object-group', 'far fa-fw'),
               title = 'Toggle dataset integration', class = 'squashed-btn')
}


#' Input form for integrating single cell datasets
#' @export
#' @keywords internal
integrationFormInput <- function(id) {
  ns <- NS(id)
  withTags({
    tagList(
      div(class = 'hidden-forms',
          div(id = ns('label-transfer-form'), class = 'hidden-form', style = 'display: none;',
              selectizeInputWithButtons(ns('ref_name'), 'Transfer labels from:',
                                        actionButton(ns('submit_transfer'), '', icon('chevron-right', 'fa-fw'), title = 'Run label transfer'),
                                        actionButton(ns('overwrite_annot'), '', icon = icon('tag', 'fa-fw'), title = 'Overwrite previous labels'),
                                        options = list(optgroupField = 'type', render = I('{option: transferLabelOption}'))
              )
          ),

          div(id = ns('integration-form'), class = 'hidden-form', style = 'display: none;',
              selectizeInput(ns('test_integration'), 'Integration test datasets:', multiple = TRUE, choices = '', width = '100%'),
              selectizeInput(ns('ctrl_integration'), 'Integration control datasets:', multiple = TRUE, choices = '', width = '100%'),
              selectizeInput(ns('exclude_clusters'), 'Integration excluded clusters:', multiple = TRUE, choices = '', width = '100%', options = list(optgroupField = 'anal')),
              textInputWithButtons(ns('integration_name'),
                                   'Name for new integrated analysis:',
                                   actionButton(ns('submit_integration'), '', icon = icon('plus', 'fa-fw'), title = 'Integrate datasets'),
                                   help_id = ns('error_msg'))
          )
      )
    )
  })
}



#' Input form and buttons to select a cluster or contrast and rename a cluster
#' @export
#' @keywords internal
clusterComparisonInput <- function(id) {
  ns <- NS(id)

  withTags({
    tagList(
      div(id = ns('selected_cluster_panel'),
          div(id = ns('select_panel'),
              div(class = 'form-group selectize-fh',
                  label(class = 'control-label', `for` = ns('selected_cluster'), 'Show marker genes for:', span(class = 'hover-info', icon('info', 'fa-fw')), title = 'Cluster (n cells :: % of total)'),
                  div(class = 'input-group',
                      div(
                        select(id = ns('selected_cluster'), style = 'display: none'),
                        script(type = 'application/json', `data-for` = ns('selected_cluster'), HTML('{}'))
                      ),
                      div(class = 'input-group-btn',

                          actionButton(ns('show_contrasts'), '',
                                       icon = icon('chevron-right', 'fa-fw'),
                                       title = 'Toggle single group comparisons'),
                          actionButton(ns('show_rename'), '',
                                       icon = icon('tag', 'fa-fw'),
                                       title = 'Toggle rename cluster')
                      )
                  )
              )
          ),
          div(id = ns('rename_panel'), style = 'display: none',
              div(class = 'form-group selectize-fh',
                  label(class = 'control-label', `for` = ns('new_cluster_name'), 'New cluster name:'),
                  div(class = 'input-group',
                      input(id = ns('new_cluster_name'), type = 'text', class = 'form-control shiny-bound-input', value = '', placeholder = ''),
                      span(class = 'input-group-btn',
                           actionButton(ns('rename_cluster'), '',
                                        icon = icon('plus', 'fa-fw'),
                                        title = 'Rename cluster')
                      )
                  )
              )
          )
      )
    )
  })
}


#' Input form to select gene for scBioGpsPlotOutput and scMarkerPlotOutput
#' @export
#' @keywords internal
selectedGeneInput <- function(id, sample_comparison = FALSE) {
  ns <- NS(id)

  exclude_ambient_button <- NULL
  if (sample_comparison)
    exclude_ambient_button <- actionButton(ns('exclude_ambient'), '',
                                           icon = icon('ban', 'fa-fw'),
                                           title = 'Toggle excluding ambient genes', class = 'squashed-btn')


  selectizeInputWithButtons(id = ns('selected_gene'),
                            label = tags$span('Show expression for:', span(class='hover-info', icon('info', 'fa-fw'))),
                            exclude_ambient_button,
                            downloadButton(ns('download'), label = NULL, icon = icon('download', 'fa-fw'), title = 'Download results'),
                            actionButton(ns('genecards'), label = NULL, icon = icon('external-link-alt', 'fa-fw'), title = 'Go to GeneCards'),
                            label_title = 'Gene (% exp test :: % exp ctrl)'
  )
}


#' Output plot of single cell clusters
#' @export
#' @keywords internal
scClusterPlotOutput <- function(id) {
  ns <- NS(id)
  plotOutput(ns('cluster_plot'), height = '500px')
}

#' Output plot of single cell markers
#' @export
#' @keywords internal
scMarkerPlotOutput <- function(id) {
  ns <- NS(id)
  plotOutput(ns('marker_plot'), height = '500px')
}

#' Output plot of biogps data for a gene
#' @export
#' @keywords internal
scBioGpsPlotOutput <- function(id) {
  ns <- NS(id)
  plotOutput(ns('biogps_plot'), height = '500px')
}
