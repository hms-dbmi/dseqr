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
