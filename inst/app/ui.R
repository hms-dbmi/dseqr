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
                ),
                div(class = 'col-sm-6 mobile-margin',
                    shiny::plotOutput(ns('heatmap'), height = '601px')
                )
            ),
            hr(),
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
           selectedAnalInput(ns('path'), label = 'Select a dataset:', with_custom = FALSE),
           selectizeInputWithButtons(ns('pathway'),
                                     label = 'Select a GO term:',
                                     label_title = 'GO term | Direction | p<10-5',
                                     actionButton(ns('direction'), '', icon = icon('arrows-alt-v', 'fa-fw'), title = 'Toggle sort direction'))
  )
}


tabs <- c('Bulk Data', 'Single Cell', 'Pathways', 'Drugs')
active <- 'Bulk Data'

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
             bulkPageUI('bulk', tab = 'Bulk Data', active),
             scPageUI("sc", tab = 'Single Cell', active),
             pathPageUI('pathways', tab = 'Pathways', active),
             drugsPageUI("drug", tab = 'Drugs', active),
             docsPageUI('docs', tab = 'docs', active)

    )
  )
)



