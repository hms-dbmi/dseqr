#' Select contrast for expression set
#'
#' @param eset \code{ExpressionSet}
#'
#' @return
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

  # ------------------- Setup

  # objects we will update
  previous <- list()

  pdata <- Biobase::pData(eset)
  pdata <- tibble::add_column(pdata, Group = NA, Title = row.names(pdata), .before = 1)

  contrast <- c()



  # ------------------- user interface

  addGroup <- "shinyjs.addGroup = function(params){
     var defaultParams = {
      rows : null,
      group : '',
      color : 'white'
    };
    params = shinyjs.getParams(params, defaultParams);

    for (var i = 0; i < params.rows.length; i++) {
      var row = params.rows[i];

      // DT first table is header, second is data
      $('#pdata table:eq(1) tr:eq('+row+') td:eq(0)')
        .text(params.group)
        .css('background-color', params.color);
    }
  }"

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



  # ------------------------- server


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

    pdata_proxy <- DT::dataTableProxy("pdata")


    # make reactive state value to keep track of ctrl vs test group
    state <- shiny::reactiveValues(ctrl = 1, pdata = 0)



    # ------------------- click 'Add'

    shiny::observeEvent(input$add, {

      # get group name
      rows  <- input$pdata_rows_selected
      group <- input$group


      # check for incomplete/wrong input
      if (group == "") {
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
        message("Selection in use with different group name.")


      } else if (state$ctrl == 1) {
        # add ctrl group data to contrast
          contrast <<- group

        # update group name prompt
        shiny::updateTextInput(session, "group", label = "Test group name:", value = "")

        # update pdata Group column
        pdata[rows, 'Group'] <<- group
        shinyjs::js$addGroup(rows, group, "#A6CEE3")
        DT::selectRows(pdata_proxy, NULL)

        # update states to trigger updates
        state$ctrl <- 0
        state$pdata <- state$pdata + 1

      } else {

        if (group %in% contrast) {
          message("Group name in use for control group.")

        } else {
          # add test group to contrast
          contrast <<- c(contrast, group)

          # update inputs
          shiny::updateTextInput(session, "group", label = "Control group name:", value = "")

          # update pdata Group column
          pdata[rows, 'Group'] <<- group

          #update states
          shinyjs::js$addGroup(rows, group, "#1F78B4")
          DT::selectRows(pdata_proxy, NULL)
          state$ctrl <- 1
          state$pdata <- state$pdata + 1
        }
      }
    })





    # ------------------- click 'Done'


    shiny::observeEvent(input$done, {
      if (state$ctrl == 0) {
        message("Need to add test group.")
      } else if (nrow(contrasts) == 0) {
        message("No contrasts selected.")
        stopApp(NULL)
      } else {
        # TODO: implement return value
        stopApp(NULL)
      }
    })


    # ------------------- click 'Reset'


    shiny::observeEvent(input$reset, {
      # TODO: handle reset
    })

  }

  shiny::runGadget(shiny::shinyApp(ui, server), viewer = shiny::paneViewer())
}
