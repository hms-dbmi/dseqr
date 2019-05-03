

#' Select contrast for an experiment.
#'
#' Used to select a contrast (control and test group samples) for an experiment.
#'
#' @param eset \code{ExpressionSet}
#'
#' @return \code{eset} with selected samples retained and column \code{group} with values \code{'control'} and \code{'test'}.
#' @export
#'
#' @examples
#'
#' eset_path <- system.file('extdata', 'IBD', 'eset.rds', package='drugseqr')
#' eset <- readRDS(eset_path)
#'
#' select_contrast(eset)
#'
#'
select_contrast <- function(eset) {

  # setup ----

  background <- '#e9305d url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAAPklEQVQoU43Myw0AIAgEUbdAq7VADCQaPyww55dBKyQiHZkzBIwQLqQzCk9E4Ytc6KEPMnTBCG2YIYMVpHAC84EnVbOkv3wAAAAASUVORK5CYII=) repeat'

  # stores group names and selected rows
  # used for validating selections and retuned when complete
  sels <- list()

  # Group column is populated by user selections
  pdata <- Biobase::pData(eset)
  pdata <- tibble::add_column(pdata, Group = NA, Title = row.names(pdata), .before = 1)


  #  user interface ----

  ui <- miniUI::miniPage(
    shinyjs::useShinyjs(),
    shiny::tags$head(
      shiny::tags$style("#pdata {white-space: nowrap;}"), # table text on 1 line
      shiny::tags$style(".dt-fake-height {height: 1px;}"), # to make 100% height div work
      shiny::tags$style("td.dt-nopad {padding: 0px !important; height: 100%;}"), # td for bg color group column
      shiny::tags$style("td.dt-nopad div {height: 100%; width: 100%; text-align: center;}"), # div inside td for bg color group column
      shiny::tags$style("td.dt-nopad span {display: inline-block; padding: 8px 10px; color: white;}") # span inside div inside td for bg color group column
    ),
    # title bar
    miniUI::gadgetTitleBar(shiny::textOutput('title', inline=TRUE), left = miniUI::miniTitleBarButton("reset", "Reset")),
    miniUI::miniContentPanel(
      shiny::fillCol(flex = c(.15, 1),
                     shiny::hr(),
                     DT::dataTableOutput("pdata")
      )
    ),
    miniUI::miniButtonBlock(
      shiny::actionButton("add", "Add Control Group", class = "btn-primary")
    )
  )

  # server ----

  server <- function(input, output, session) {

    # make reactive state value to keep track of ctrl vs test group
    state <- shiny::reactiveValues(ctrl = 1, pdata = 0)

    # reactive label based on ctrl vs test state
    label_r <- shiny::eventReactive(state$ctrl, {
      if (state$ctrl == 1) {
        return('Add Control Rows')

      } else if (state$ctrl == 0) {
        return('Add Test Rows')

      } else if (state$ctrl == 3) {
        return('Click Reset or Done')
      }
    })

    # title is dynamic label
    output$title <- shiny::renderText((
      gsub('^Add', 'Select', label_r())
    ))

    # label updates should update add button
    shiny::observe({
      shiny::updateActionButton(session, "add", label = label_r())
    })


    # show phenotype data
    output$pdata <- DT::renderDataTable({

      DT::datatable(
        isolate(pdata_r()),
        class = 'cell-border dt-fake-height',
        rownames = FALSE,
        escape = FALSE, # to allow HTML in table
        options = list(
          columnDefs = list(list(className = 'dt-nopad', targets = 0)),
          scrollY = FALSE,
          scrollX = TRUE,
          paging = FALSE,
          bInfo = 0
        )
      )
    })

    # pdata reactive so that will update group column
    pdata_r <- shiny::eventReactive(state$pdata, {
      return(pdata)
    })

    # proxy used to replace data
    pdata_proxy <- DT::dataTableProxy("pdata")
    shiny::observe({
      DT::replaceData(pdata_proxy, pdata_r(), rownames = FALSE)
    })


    # click 'Add' ----

    shiny::observeEvent(input$add, {

      # get rows
      rows  <- input$pdata_rows_selected

      # check for incomplete/wrong input
      if (length(sels) == 2) {
        message("Contrast is fully specified. Click 'Reset' to clear selections or 'Done' to exit.")

      } else if (length(rows) == 0) {
        message("Select rows.")

      } else if (Position(function(x) setequal(x, rows), sels, nomatch = 0) > 0) {
        message("Selection in use for control group.")

      } else if (state$ctrl == 1) {
        # add ctrl group data to sels
        sels <<- c(sels, list(ctrl = rows))

        # update pdata Group column
        pdata[rows, 'Group'] <<- paste0('<div style="background: ', background, ';"><span>control</span></div>')

        # update states to trigger updates
        state$ctrl <- 0
        state$pdata <- state$pdata + 1

      } else {
          # add test group data to sels
          sels <<- c(sels, list(test = rows))

          # add group to table
          pdata[rows, 'Group'] <<- '<div style="background-color: #e6194b;"><span>test</span></div>'
          state$pdata <- state$pdata + 1
          state$ctrl <- 3
          shinyjs::disable("add")
      }
    })

    # click 'Done' ----

    shiny::observeEvent(input$done, {
      if (state$ctrl == 0) {
        message("Need to add test group.")
      } else if (length(sels) == 0) {
        message("No contrast selected.")
      } else {
        group <- rep(NA, ncol(eset))
        group[sels$test] <- 'test'
        group[sels$ctrl] <- 'control'
        Biobase::pData(eset)$group <- group

        # retain selected samples only
        shiny::stopApp(eset[, !is.na(group)])
      }
    })


    # click 'Reset' ----

    shiny::observeEvent(input$reset, {
      sels <<- list()
      pdata$Group <<- NA

      # remove groups from table and reset control
      state$ctrl <- 1
      state$pdata <- 0
      shinyjs::enable("add")
    })
  }

  shiny::runGadget(shiny::shinyApp(ui, server), viewer = shiny::paneViewer())
}
