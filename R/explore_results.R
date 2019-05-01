
explore_results <- function(query_res) {

  # setup ----

  #  user interface ----

  ui <- miniUI::miniPage(
    shinyjs::useShinyjs(),
    shiny::tags$head(
      shiny::tags$style("#query_res {white-space: nowrap;}") # table text on 1 line
    ),
    # title bar
    miniUI::gadgetTitleBar("Explore Results"),
    miniUI::miniContentPanel(
      shiny::fillCol(flex = c(NA, NA, 1),
                     shiny::selectizeInput('study',
                                           'Select study:',
                                           choices = c('CMAP02', 'L1000')),
                     shiny::hr(),
                     DT::dataTableOutput("query_res")
      )
    )
  )

  # server ----

  server <- function(input, output, session) {

    # show query data
    output$query_res <- DT::renderDataTable({

      DT::datatable(
        query_res,
        class = 'cell-border',
        rownames = FALSE,
        escape = FALSE, # to allow HTML in table
        options = list(
          scrollY = TRUE,
          paging = FALSE,
          bInfo = 0
        )
      )
    })


    # click 'Done' ----

    shiny::observeEvent(input$done, {
      shiny::stopApp()
    })

  }

  shiny::runGadget(shiny::shinyApp(ui, server), viewer = shiny::browserViewer())
}
