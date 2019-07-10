#' UI for navbar
#' @param active the active tab name
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
            )
        )
    )
  })
}

#' UI for Single Cell Exploration page
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

BulkPageUI <- function(id, tab, active) {
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


DrugPageUI <- function(id, tab, active) {

  ns <- NS(id)

  active_class <- ifelse(tab == active, 'active', '')


  withTags({
    div(class = paste('tab-pane', active_class), `data-value` = tab, id = id_from_tab(tab),
        div(class = 'row',
            div(class = 'col-sm-6',
                div(class = "well-form well-bg",
                    div(class = 'form-group selectize-fh',
                        label(class = 'control-label', `for` = ns('selected_anal'), 'Select study:'),
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
                    ),
                      div(id = ns('advanced-panel'), class = 'hidden-form', style = 'display: none;',
                        selectizeInput(ns('cells'), 'Select cell lines:', choices = NULL, multiple = TRUE, options = list(placeholder = "showing all"), width = '100%'),
    # checkboxGroupInput('options', NULL, c('Plot Histogram', 'Blah')),
                        plotOutput(outputId = ns("histPlot"), width = 800)
                        )

                )
            )
        ),
        hr(),
        div(class = 'dt-container',
            DT::dataTableOutput(ns("query_res"))
        )
    )
  })
}

#' Input form for Single Cell Exploration page
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
showIntegrationButton <- function(id) {
  ns <- NS(id)

  actionButton(ns('show_integration'), '',
               icon = icon('object-group', 'far fa-fw'),
               title = 'Toggle dataset integration', class = 'squashed-btn')
}


#' Input form for integrating single cell datasets
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

                div(style = 'display: inline-block',
                    uiOutput(ns("genecards"))

                )
            )
        )
    )
  })
}


#' Output plot of single cell clusters
scClusterPlotOutput <- function(id) {
  ns <- NS(id)
  plotOutput(ns('cluster_plot'))
}

#' Output plot of single cell markers
scMarkerPlotOutput <- function(id) {
  ns <- NS(id)
  plotOutput(ns('marker_plot'))
}

#' Output plot of biogps data for a gene
scBioGpsPlotOutput <- function(id) {
  ns <- NS(id)
  plotOutput(ns('biogps_plot'))
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
             BulkPageUI("bulk", tab = 'Bulk', active),
             DrugPageUI("drug", tab = 'Drugs', active)
    )
  )
)
