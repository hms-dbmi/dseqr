scPageInput <- function(id) {
  ns <- NS(id)
  withTags({
    div(class = 'tab-pane active', `data-value` = 'Single-Cell', id = 'single-cell-tab',
        div(class = 'row',
            div(class = 'col-sm-6',
                scFormInput(ns('form'))
            ),
            div(class = 'col-sm-6',
                scClusterPlotUI(ns('cluster_plot'))
            )
        ),
        hr(),
        div(class = 'row',
            div(class = "col-sm-6 col-lg-6 col-lg-push-6",
                scMarkerPlotUI(ns('marker_plot'))
            ),
            div(class = "col-sm-6 col-lg-6 col-lg-pull-6",
                scBioGpsPlotUI(ns('biogps_plot'))
            )
        )
    )
  })
}

scFormInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(class = "well-form well-bg",
        selectedAnalInput(ns('anal')),
        integrationFormInput(ns('integration')),
        br(),
        selectedClusterInput(ns('cluster')),
        br(),
        selectedGeneInput(ns('gene')),
        br(),
        selectedGroupsInput(ns('groups'))
    )
  })
}

scClusterPlotUI <- function(id) {
  ns <- NS(id)
  plotOutput(ns('cluster_plot'))
}

scMarkerPlotUI <- function(id) {
  ns <- NS(id)
  plotOutput(ns('marker_plot'))
}

scBioGpsPlotUI <- function(id) {
  ns <- NS(id)
  plotOutput(ns('biogps_plot'))
}

integrationFormInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(id = ns('integration-form'), class = 'hidden-form', style = 'display: none;',
        textInput(ns('integration_name'), 'Name for new integrated analysis:', width = '100%'),
        selectizeInput(ns('test_integration'), 'Test datasets:', multiple = TRUE, choices = '', width = '100%'),
        selectizeInput(ns('ctrl_integration'), 'Control datasets:', multiple = TRUE, choices = '', width = '100%'),

        div(class = 'validate-wrapper',
            actionButton(ns('submit_integration'), 'Integrate Datasets',
                         icon = shiny::icon('object-group', 'fa-fw'),
                         title = 'Integrate datasets',
                         class = 'btn-block btn-default'),
            span(class = 'help-block')
        )
    )
  })
}


selectedAnalInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(class = 'form-group selectize-fh',
        label(class = 'control-label', `for` = ns('selected_anal'), 'Select a dataset:'),
        div(class = 'input-group',
            div(
              select(id = ns('selected_anal')),
              script(type = 'application/json', `data-for` = ns('selected_anal'), HTML('{}'))
            ),
            div(class = 'input-group-btn',
                showIntegrationInput(ns('integration')),
                plotStylesInput(ns('styles'))
            )
        )
    )
  })
}

plotStylesInput <- function(id) {
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


showIntegrationInput <- function(id) {
  ns <- NS(id)

  actionButton(ns('show_integration'), '',
               icon = icon('object-group', 'far fa-fw'),
               title = 'Toggle dataset integration', class = 'squashed-btn')
}


selectedClusterInput <- function(id) {
  ns <- NS(id)

  withTags({
    tagList(

      div(id = ns('select_panel'),
          div(class = 'form-group selectize-fh',
              label(class = 'control-label', `for` = ns('selected_cluster'), 'Show marker genes for:'),
              div(class = 'input-group',
                  div(
                    select(id = ns('selected_cluster')),
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
  })
}


selectedGeneInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(class = 'form-group selectize-fh',
        label(class = 'control-label', `for` = ns('selected_gene'), 'Show expression for:'),
        div(class = 'input-group',
            div(
              select(id = ns('selected_gene')),
              script(type = 'application/json', `data-for` = ns('selected_gene'), HTML('{}'))
            ),
            div(class = 'input-group-btn',
                uiOutput(ns("genecards"))

            )
        )
    )
  })
}

selectedGroupsInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(style = 'display: block;',
        shinyWidgets::radioGroupButtons(ns('selected_group'), "Show cells for group:",
                                        choices = c('test', 'all', 'ctrl'),
                                        selected = 'all', justified = TRUE)

    )

  })

}

## ui.R ##
bootstrapPage(
  useShinyjs(),
  includeScript(path = 'www/contrasts.js'),
  includeCSS(path = 'www/custom.css'),
  htmlTemplate("navbar.html"),
  fluidPage(
    # make sure selectize loaded (not using default)
    tags$div(style = 'display: none', selectizeInput('blah1', label = NULL, choices = '')),

    # THE TABS ----
    tags$div(class = "tab-content", `data-tabsetid` = "1823", id = "tabs",

             # single cell tab
             scPageInput("sc")


    )
  )
)
