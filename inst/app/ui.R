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
            div(class = 'scroll-plot',
                plotlyOutput(ns('path_plot'), height = '550px')

            )
    )
  })
}


selectizeInputWithButtons <- function(id, label, ..., options = NULL) {

  options <- ifelse(is.null(options), '{}', jsonlite::toJSON(options, auto_unbox = TRUE))

  tags$div(class = 'form-group selectize-fh',
           tags$label(class = 'control-label', `for` = id, label),
           tags$div(class = 'input-group full-height-btn',
                    tags$div(class = 'full-height-selectize',
                      tags$select(id = id, style = 'display: none'),
                      tags$script(type = 'application/json', `data-for` = id, HTML(options))
                    ),
                    tags$div(class = 'input-group-btn',
                             # the buttons
                             ...
                    )
           )
  )
}


#' Input form for pathways page
#' @export
#' @keywords internal
pathFormInput <- function(id) {
  ns <- NS(id)

    tags$div(class = "well-form well-bg",
        selectizeInputWithValidation(ns('anal'), 'Select an analysis:', options = list(optgroupField= 'dataset_name')),
        selectizeInputWithButtons(ns('pathway'), 'Select a pathway:', actionButton(ns('kegg'), '', icon = icon('external-link-alt', 'fa-fw'), title = 'Go to KEGG'))
    )
}


tabs <- c('Datasets', 'Single Cell', 'Pathways', 'Drugs')
active <- 'Pathways'

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
