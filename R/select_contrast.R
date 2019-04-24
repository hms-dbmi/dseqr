select_contrast <- function(eset) {

  # ------------------- Setup

  # striped background for test group
  group_colors <- c("#A6CEE3", "#1F78B4")
  background <- '#51bd64 url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAAPklEQVQoU43Myw0AIAgEUbdAq7VADCQaPyww55dBKyQiHZkzBIwQLqQzCk9E4Ytc6KEPMnTBCG2YIYMVpHAC84EnVbOkv3wAAAAASUVORK5CYII=) repeat'

  # objects we will update
  previous <- list()

  pdata <- Biobase::pData(eset)
  pdata <- tibble::add_column(pdata, Group = c(rep('A', 6), rep('B', 6)), Title = row.names(pdata), .before = 1)

  contrasts <- data.frame(Control = character(0),
                          Test = character(0), stringsAsFactors = FALSE)



  # ------------------- user interface

  ui <- miniUI::miniPage(
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

      # color control/test group for easy identification
      groups <- unique(pdata$Group)
      dt <- DT::formatStyle(dt, 'Group', target = 'cell',
                            backgroundColor = DT::styleEqual(c(NA, groups),
                                                             c('#FFFFFF',
                                                               group_colors[seq_along(groups)])))

      return(dt)
    })


    # make reactive state value to keep track of ctrl vs test group
    state <- shiny::reactiveValues(ctrl = 1, contrast = 0)



    # ------------------- click 'Add'

    shiny::observeEvent(input$add, {

      # construct group data
      rows  <- input$pdata_rows_selected
      group <- input$group
      group_data <- list()
      group_data[[group]] <- rows


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
        # add ctrl group data to previous and contrasts
        if (!group %in% names(previous))
          previous <<- c(previous, group_data)
        contrasts[nrow(contrasts) + 1, ] <<- c(group, NA)

        # update inputs
        shiny::updateTextInput(session, "group",
                        label = "Test group name:", value = "")
        shiny::updateSelectInput(session, "prev",
                          choices = c("", names(previous)))

        # update states
        state$ctrl <- 0
        state$contrast <- state$contrast + 1


      } else {

        if (group == contrasts[nrow(contrasts), "Control"]) {
          message("Group name in use for control group.")

        } else {
          # add test group data to previous and contrasts
          if (!group %in% names(previous))
            previous <<- c(previous, group_data)
          contrasts[nrow(contrasts), "Test"] <<- group

          # update inputs
          shiny::updateTextInput(session, "group",
                                 label = "Control group name:", value = "")
          shiny::updateSelectInput(session, "prev",
                                   choices = c("", names(previous)))

          #update states
          state$ctrl <- 1
          state$contrast <- state$contrast + 1
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
