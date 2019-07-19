


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

#' @export
#' @keywords internal
bulkPageUI <- function(id, tab, active) {
  ns <- NS(id)
  active_class <- ifelse(tab == active, 'active', '')
  withTags({
    div(class = paste('tab-pane', active_class), `data-value` = tab, id = id_from_tab(tab),
        div(class = 'row'),
        hr(),
        div(class = 'row')
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

  withTags({
    div(class = 'form-group selectize-fh',
        label(class = 'control-label', `for` = ns('selected_clusters'), 'Compare samples for:'),
        div(class = 'input-group full-height-btn',
            div(class = 'full-height-selectize',
                select(id = ns('selected_clusters'), style = 'display: none', multiple = TRUE),
                script(type = 'application/json', `data-for` = ns('selected_clusters'), HTML('{}'))
            ),
            div(class = 'input-group-btn',
                actionButton(ns('run_comparison'), '',
                             icon = icon('chevron-right', 'far fa-fw'),
                             title = 'Compare test to control cells')
            )
        )
    )
  })
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

  withTags({
    div(class = 'form-group selectize-fh',
        label(class = 'control-label', `for` = ns('selected_anal'), 'Select a dataset:'),
        div(class = 'input-group',
            div(
              select(id = ns('selected_anal'), style = 'display: none'),
              script(type = 'application/json', `data-for` = ns('selected_anal'), HTML('{}'))
            ),
            div(class = 'input-group-btn',
                showIntegrationButton(ns('integration')),
                plotStylesButton(ns('styles'))
            )
        )
    )
  })
}

#' Button with sliders for adjusting plot jitter and point size
#' @export
#' @keywords internal
plotStylesButton <- function(id) {
  ns <- NS(id)
  shinyWidgets::dropdownButton(
    sliderInput(ns('point_size'), 'Point size:',
                width = '100%', ticks = FALSE,
                min = 0.5, max = 4, value = 2.5, step = 0.5),
    sliderInput(ns('point_jitter'), 'Point jitter:',
                width = '100%', ticks = FALSE,
                min = 0, max = 3, value = 0, step = 0.5),
    circle = FALSE, right = TRUE, icon = icon('cog', 'fa-fw')
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
    div(id = ns('integration-form'), class = 'hidden-form', style = 'display: none;',
        selectizeInput(ns('test_integration'), 'Test datasets:', multiple = TRUE, choices = '', width = '100%'),
        selectizeInput(ns('ctrl_integration'), 'Control datasets:', multiple = TRUE, choices = '', width = '100%'),
        div(class = 'form-group selectize-fh',
            label(class = 'control-label', `for` = ns('integration_name'), 'Name for new integrated analysis:'),
            div(class = 'validate-wrapper', id = ns('validate'),

                div(class = 'input-group',
                    input(id = ns('integration_name'), type = 'text', class = 'form-control shiny-bound-input', value = '', placeholder = ''),
                    span(class = 'input-group-btn',
                         actionButton(ns('submit_integration'), '',
                                      icon = icon('plus', 'fa-fw'),
                                      title = 'Integrate datasets')
                    )
                ),
                span(id = ns('error_msg'), class = 'help-block')
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
                  label(class = 'control-label', `for` = ns('selected_cluster'), 'Show marker genes for:'),
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


  withTags({
    div(class = 'form-group selectize-fh',
        label(class = 'control-label', `for` = ns('selected_gene'), 'Show expression for:'),
        div(class = 'input-group',
            div(
              select(id = ns('selected_gene'), style = 'display: none;'),
              script(type = 'application/json', `data-for` = ns('selected_gene'), HTML('{}'))
            ),
            div(class = 'input-group-btn',
                exclude_ambient_button,
                actionButton(ns('genecards'), label = NULL, icon = icon('external-link-alt', 'fa-fw'))
            )
        )
    )
  })
}


#' Output plot of single cell clusters
#' @export
#' @keywords internal
scClusterPlotOutput <- function(id) {
  ns <- NS(id)
  plotOutput(ns('cluster_plot'))
}

#' Output plot of single cell markers
#' @export
#' @keywords internal
scMarkerPlotOutput <- function(id) {
  ns <- NS(id)
  plotOutput(ns('marker_plot'))
}

#' Output plot of biogps data for a gene
#' @export
#' @keywords internal
scBioGpsPlotOutput <- function(id) {
  ns <- NS(id)
  plotOutput(ns('biogps_plot'))
}