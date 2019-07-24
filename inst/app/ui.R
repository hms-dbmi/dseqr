
tabs <- c('Datasets', 'Single Cell', 'Pathways', 'Drugs')
active <- 'Single Cell'


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

selectedGeneInput <- function(id, sample_comparison = FALSE) {
  ns <- NS(id)

  exclude_ambient_button <- NULL
  if (sample_comparison)
    exclude_ambient_button <- actionButton(ns('exclude_ambient'), '',
                                           icon = icon('ban', 'fa-fw'),
                                           title = 'Toggle excluding ambient genes', class = 'squashed-btn')


  selectizeInputWithButtons(id = ns('selected_gene'),
                            label = 'Show expression for:',
                            exclude_ambient_button,
                            downloadButton(ns('download'), label = NULL, icon = icon('download', 'fa-fw'), title = 'Download'),
                            actionButton(ns('genecards'), label = NULL, icon = icon('external-link-alt', 'fa-fw'), title = 'Go to GeneCards')
  )
}


selectizeInputWithButtons <- function(id, label, ..., options = NULL) {

  mult <- isTRUE(options$multiple)
  if(mult) {
    select_tag <- tags$select(id = id, style = 'display: none', multiple = TRUE)
  } else {
    select_tag <- tags$select(id = id, style = 'display: none')
  }

  buttons <- list(...)
  buttons <- buttons[!sapply(buttons, is.null)]

  options <- ifelse(is.null(options), '{}', jsonlite::toJSON(options, auto_unbox = TRUE))

  tags$div(class = 'form-group selectize-fh',
           tags$label(class = 'control-label', `for` = id, label),
           tags$div(class = 'input-group full-height-btn',
                    tags$div(class = 'full-height-selectize',
                             select_tag,
                             tags$script(type = 'application/json', `data-for` = id, HTML(options))
                    ),
                    lapply(buttons, function(btn) {
                      tags$div(class = 'input-group-btn',
                               # the buttons
                               btn
                      )
                    })
           )
  )
}

bootstrapPage(
  useShinyjs(),
  includeScript(path = 'www/contrasts.js'),
  includeScript(path = 'www/cellOptions.js'),
  includeScript(path = 'www/pathOptions.js'),
  includeScript(path = 'www/toggleClinicalTitle.js'),
  includeCSS(path = 'www/custom.css'),
  includeCSS(path = 'www/drugs.css'),
  includeCSS(path = 'www/pathways.css'),
  navbarUI(tabs, active),
  fluidPage(
    # make sure css/js loaded from packages where not using functions (not using default)
    tags$div(style = 'display: none;', selectizeInput('blah1', label = NULL, choices = '')),
    tags$div(style = 'display: none;', shinyFiles::shinyDirButton("blah2", title='', label='', icon=icon('plus'))),

    tags$div(class = "tab-content", `data-tabsetid` = "tabset", id = "tabs",
             # single cell tab
             scPageUI("sc", tab = 'Single Cell', active),
             drugsPageUI("drug", tab = 'Drugs', active),
             dsPageUI('datasets', tab = 'Datasets', active),
             pathPageUI('pathways', tab = 'Pathways', active)
    )
  )
)

