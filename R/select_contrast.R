#' Select contrast for an experiment.
#'
#' Used to select a contrast (control and test group samples) for an experiment.
#'
#' @param eset \code{ExpressionSet}
#'
#' @return Named list of length two. In order, list items contain rows corresponding to control and test group. Names are specified group names.
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

  # stores group names and selected rows
  # used for validating selections and retuned when complete
  previous <- list()

  # Group column is populated by user selections
  pdata <- Biobase::pData(eset)
  pdata <- tibble::add_column(pdata, Group = NA, Title = row.names(pdata), .before = 1)


  # js to update table ----

  addGroup <- "shinyjs.addGroup = function(params){

     // DT first table is header, second is data
     var pdata = $('#pdata table:eq(1)')

     // setup default params (used for Reset)
     var nrow = $('tr', pdata).length


     var allRows = [];
     for (var i = 1; i <= nrow; i++) {
       allRows.push(i);
     }

     var defaultParams = {
      rows: allRows,
      group: '',
      color: 'transparent'
    };

    params = shinyjs.getParams(params, defaultParams);

    // add group name and color to selected rows
    for (var i = 0; i < params.rows.length; i++) {

      $('tr:eq('+ params.rows[i] +') td:eq(0)', pdata)
        .text(params.group)
        .css('background-color', params.color);
    }
  }"


  #  user interface ----

  ui <- miniUI::miniPage(
    shinyjs::useShinyjs(),
    shinyjs::extendShinyjs(text = addGroup, functions = c('addGroup')),
    shiny::tags$head(
      shiny::tags$style("#pdata {white-space: nowrap;}") # table text on 1 line
    ),
    # title bar
    miniUI::gadgetTitleBar("Select Contrast", left = miniUI::miniTitleBarButton("reset", "Reset")),
    miniUI::miniContentPanel(
      shiny::fillCol(flex = c(NA, NA, 1),
                     shiny::textInput(
                       "group",
                       "Control group name:",
                       width = "200px"
                     ),
                     shiny::hr(),
                     DT::dataTableOutput("pdata")
      )
    ),
    miniUI::miniButtonBlock(
      shiny::actionButton("add", "Add Group")
    )
  )

  # server ----

  server <- function(input, output, session) {

    # show phenotype data
    output$pdata <- DT::renderDataTable({

      dt <- DT::datatable(
        pdata,
        class = 'cell-border',
        rownames = FALSE,
        options = list(
          scrollY = FALSE,
          paging = FALSE,
          bInfo = 0
        )
      )
    })

    # proxy is used to deselect rows after adding a group
    pdata_proxy <- DT::dataTableProxy("pdata")


    # make reactive state value to keep track of ctrl vs test group
    state <- shiny::reactiveValues(ctrl = 1)


    # click 'Add' ----

    shiny::observeEvent(input$add, {

      # get group name
      rows  <- input$pdata_rows_selected
      group <- input$group
      group_data <- list()
      group_data[[group]] <- rows

      print(group_data)
      print(previous)


      # check for incomplete/wrong input
      if (length(previous) == 2) {
        message("Contrast is fully specified. Click 'Reset' to clear selections or 'Done' to exit.")

      } else if (group == "") {
        message("Enter group name.")

      } else if (length(rows) == 0) {
        message("Select rows.")

      } else if (make.names(group) != group) {
        message("Group name invalid.")

      } else if (group %in% names(previous) &&
                 !setequal(previous[[group]], rows)) {
        message("Group name in use with different samples.")

      } else if (Position(function(x) setequal(x, rows), previous, nomatch = 0) > 0 &&
                 !group %in% names(previous)) {
        message("Selection in use for control group.")

      } else if (state$ctrl == 1) {
        # add ctrl group data to previous
        previous <<- c(previous, group_data)

        # update group name prompt
        shiny::updateTextInput(session, "group", label = "Test group name:", value = "")

        # update pdata Group column
        pdata[rows, 'Group'] <<- group
        shinyjs::js$addGroup(rows, group, "#A6CEE3")
        DT::selectRows(pdata_proxy, NULL)

        # update states to trigger updates
        state$ctrl <- 0

      } else {

        if (group %in% names(previous)) {
          message("Group name in use for control group.")

        } else {
          # add test group data to previous
          previous <<- c(previous, group_data)

          # update inputs
          shiny::updateTextInput(session, "group", label = "Control group name:", value = "")

          # add group to table
          shinyjs::js$addGroup(rows, group, "#1F78B4")
          DT::selectRows(pdata_proxy, NULL)
          state$ctrl <- 1
        }
      }
    })

    # click 'Done' ----

    shiny::observeEvent(input$done, {
      if (state$ctrl == 0) {
        message("Need to add test group.")
      } else if (length(previous) == 0) {
        message("No contrast selected.")
      } else {
        # TODO: implement return value
        stopApp(previous)
      }
    })


    # click 'Reset' ----

    shiny::observeEvent(input$reset, {
      previous <<- list()

      # remove groups from table and reset control
      shinyjs::js$addGroup()
      state$ctrl <- 1
      shiny::updateTextInput(session, "group", label = "Control group name:", value = "")
      # TODO: handle reset
    })
  }

  shiny::runGadget(shiny::shinyApp(ui, server), viewer = shiny::paneViewer())
}
