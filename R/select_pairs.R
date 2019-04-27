

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
#' data_dir <- system.file('extdata', 'IBD', package='drugseqr')
#' pdata_path <- file.path(data_dir, 'Phenotypes.csv')
#'
#' select_pairs(data_dir, pdata_path)
#'
#'
select_pairs <- function(data_dir, pdata_path) {

  # setup ----

  background <- '#e9305d url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAAPklEQVQoU43Myw0AIAgEUbdAq7VADCQaPyww55dBKyQiHZkzBIwQLqQzCk9E4Ytc6KEPMnTBCG2YIYMVpHAC84EnVbOkv3wAAAAASUVORK5CYII=) repeat'

  # load pdata and determine row to file correspondence
  # needs to be data.frame for ExpressionSet construction
  pdata <- tryCatch(data.table::fread(pdata_path, fill=TRUE, data.table = FALSE),
                    error = function(err) {err$message <- "Couldn't read pdata"; stop(err)})

  # match pdata and file names
  fastqs <- list.files(data_dir, '.fastq.gz')
  pdata <- match_pdata(pdata, fastqs)

  pdata <- tibble::add_column(pdata, Pair = NA, Replicate = NA, .before = 1)

  # auto-detect if paired
  fastq_id1s <- get_fastq_id1s(file.path(data_dir, fastqs))
  paired <- detect_paired(fastq_id1s)

  # select and mark auto-detected pair type
  end_types <- c('single-ended', 'pair-ended')
  if (paired) end_types <- end_types[c(2, 1)]
  end_types[1] <- paste(end_types[1], '(detected)')

  message("For paired-end experiments:
            - select and add paired rows (usually contain _1 and _2 in filename)\n")

  message("For paired-end and single-ended experiments:
            - confirm that file names (in column 'File Name') are assigned to the correct sample rows
            - confirm that the experiment was correctly identified as single-ended or paired-end
            - select samples to treat as a single library (if any - e.g. same sample sequenced in replicate)\n")


  # things we will update/return to user
  pairs <- list()
  dups <- list()


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
      shiny::fillCol(flex = c(NA, NA, 1),
                     shiny::selectizeInput('end_type',
                                           'Confirm end-type:',
                                           choices = end_types),
                     shiny::hr(),
                     DT::dataTableOutput("pdata")
      )
    ),
    miniUI::miniButtonBlock(
      shiny::actionButton("pairs", "Pair Samples"),
      shiny::actionButton("dups", "Mark Replicates")
    )
  )

  # server ----

  server <- function(input, output, session) {

    # make reactive state value to keep track of ctrl vs test group
    state <- shiny::reactiveValues(ctrl = 1, pdata = 0)



    # show phenotype data
    output$pdata <- DT::renderDataTable({

      DT::datatable(
        isolate(pdata_r()),
        class = 'cell-border dt-fake-height',
        rownames = FALSE,
        escape = FALSE, # to allow HTML in table
        options = list(
          columnDefs = list(list(className = 'dt-nopad', targets = c(0, 1))),
          scrollY = FALSE,
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
        stopApp(eset[, !is.na(group)])
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
