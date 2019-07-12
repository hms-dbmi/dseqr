datasetsPage <- function(input, output, session) {
  dsForm <- callModule(datasetsForm, 'form')

  callModule(dsPairsTable, 'pairs', fastq_dir = dsForm$fastq_dir)
}


datasetsForm <- function(input, output, session) {

  dsAnalysis <- callModule(dsAnalysis, 'selected_anal')

  end_type <- callModule(dsEndType, 'end_type',
                         fastq_dir = dsAnalysis$fastq_dir)

  return(list(
    fastq_dir = dsAnalysis$fastq_dir
  ))

}

dsAnalysis <- function(input, output, session) {

  indices_dir <- '~/Documents/Batcave/zaklab/drugseqr/data-raw/indices/kallisto'


  # get directory with fastqs
  bulk_dir <- c('bulk-data' = '~/Documents/Batcave/zaklab/drugseqr/data-raw/bulk/example-data')
  shinyFiles::shinyDirChoose(input, "anal_dir", roots = bulk_dir)
  fastq_dir <- reactiveVal('~/Documents/Batcave/zaklab/drugseqr/data-raw/bulk/example-data/IBD')

  # fastq_dir <- reactive({
  #   dir <- parseDirPath(bulk_dir, input$anal_dir)
  #   req(dir)
  #   dir
  # })

  anal_name <- reactive(input$anal_name)

  return(list(
    anal_name = anal_name,
    fastq_dir = fastq_dir
  ))

}


dsEndType <- function(input, output, session, fastq_dir) {

  # get fastq files in directory
  fastq_files <- reactive({
    fastq_dir <- fastq_dir()
    req(fastq_dir)

    list.files(fastq_dir, '.fastq.gz$')
  })

  # auto detected if paired
  detected_paired <- reactive({
    fastqs <- fastq_files()
    fastq_dir <- fastq_dir()
    req(fastqs, fastq_dir)

    # auto-detect if paired
    fastq_id1s <- get_fastq_id1s(file.path(fastq_dir, fastqs))
    detect_paired(fastq_id1s)
  })



  observe({
    # label the end type choices with auto detected
    end_types <- c('single-ended' = 'single-ended', 'pair-ended' = 'pair-ended')
    if (detected_paired()) end_types <- end_types[c(2, 1)]
    names(end_types)[1] <- paste(names(end_types)[1], '(detected)')

    updateSelectizeInput(session, 'end_type', choices = end_types)
  })

  return(end_type = reactive(input$end_type))
}


dsPairsTable <- function(input, output, session, fastq_dir) {

  # things user will update and return
  state <- reactiveValues(pdata = 0)
  pdata_r <- reactiveVal()
  pairs_r <- reactiveVal()
  reps_r <- reactiveVal()

  # setup ----

  # label button clicks
  labels <- callModule(dsLabelRows, 'label_rows')

  group_colors <- c("#C7E9C0", "#C6DBEF", "#FCBBA1", "#FDD0A2", "#BCBDDC", "#D9D9D9", "#F6E8C3", "#DC143C",
                    "#A1D99B", "#9ECAE1", "#FC9272", "#FDAE6B", "#9E9AC8", "#BDBDBD", "#DFC27D", "#FFFFFF",
                    "#C3B091", "#007FFF", "#00FFFF", "#7FFFD4", "#228B22", "#808000", "#7FFF00", "#BFFF00",
                    "#FFD700", "#DAA520", "#FF7F50", "#FA8072","#FC0FC0", "#CC8899", "#E0B0FF", "#B57EDC", "#843179")

  ncolors <- length(group_colors)
  background <- 'url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAAPklEQVQoU43Myw0AIAgEUbdAq7VADCQaPyww55dBKyQiHZkzBIwQLqQzCk9E4Ytc6KEPMnTBCG2YIYMVpHAC84EnVbOkv3wAAAAASUVORK5CYII=) repeat'


  # reset everything when new fastq_dir
  observeEvent(fastq_dir(), {
    fastq_dir <- fastq_dir()
    req(fastq_dir)

    fastqs <- list.files(fastq_dir, '.fastq.gz$')
    pdata <- tibble::tibble('File Name' = fastqs)
    pdata <- tibble::add_column(pdata, Group = NA, Pair = NA, Replicate = NA, .before = 1)

    pdata_r(pdata)

    clear <- rep(NA, nrow(pdata))
    pairs_r(clear)
    reps_r(clear)
  })

  # redraw table when new pdata (otherwise update data using proxy)
  output$pdata <- DT::renderDataTable({

    DT::datatable(
      pdata_r(),
      class = 'cell-border dt-fake-height',
      rownames = FALSE,
      escape = FALSE, # to allow HTML in table
      options = list(
        columnDefs = list(list(className = 'dt-nopad', targets = c(0, 1, 2))),
        scrollY = FALSE,
        paging = FALSE,
        bInfo = 0
      )
    )
  })


  # proxy pdata to update Pair/Replicate column
  proxy_pata <- shiny::eventReactive(state$pdata, {

    pdata <- pdata_r()
    reps <- reps_r()
    pairs <- pairs_r()

    # update pdata Replicate column
    rep_nums <- sort(unique(setdiff(reps, NA)))
    for (rep_num in rep_nums) {
      color <- group_colors[ncolors - rep_num]
      rows <- which(reps == rep_num)
      pdata[rows, 'Replicate'] <- paste('<div style="background-color:', color, ';"></div>')
    }

    # update pdata Pair column
    if (TRUE) {
      pair_nums <- sort(unique(setdiff(pairs, NA)))
      for (pair_num in pair_nums) {
        color <- group_colors[pair_num]
        rows <- which(pairs == pair_num)
        pdata[rows, 'Pair'] <- paste('<div style="background:', color, background, ';"></div>')
      }
    } else {
      pdata[1:nrow(pdata), 'Pair'] <- NA
    }

    return(pdata)
  })

  # proxy used to replace data
  proxy <- DT::dataTableProxy("pdata")
  shiny::observe({

    DT::replaceData(proxy, proxy_pata(), rownames = FALSE)
  })


  # click 'Pair Samples' ----

  shiny::observeEvent(labels$pair(), {
    req(labels$pair())
    reps <- reps_r()
    pairs <- pairs_r()

    # get rows
    rows  <- input$pdata_rows_selected

    # check for incomplete/wrong input
    if (validate_pairs(pairs, rows, reps)) {

      # add rows as a pair
      pair_num <- length(unique(setdiff(pairs, NA))) + 1
      pairs[rows] <- pair_num
      pairs_r(pairs)

      # update states to trigger updates
      state$pdata <- state$pdata + 1
    }
  })

  # click 'Mark Replicates' ----

  shiny::observeEvent(labels$rep(), {
    req(labels$rep())
    reps <- reps_r()
    pairs <- pairs_r()

    # get rows
    rows  <- input$pdata_rows_selected
    if (validate_reps(pairs, rows, reps)) {
      # add rows as replicates
      rep_num <- length(unique(setdiff(reps, NA))) + 1
      reps[rows] <- rep_num
      reps_r(reps)

      # update states to trigger updates
      state$pdata <- state$pdata + 1
    }
  })


  # click 'Reset' ----

  shiny::observeEvent(input$reset, {
    pdata <- pdata()
    pairs <<- rep(NA, nrow(pdata))
    reps <<- rep(NA, nrow(pdata))
    pdata$Pair <<- NA
    pdata$Replicate <<- NA

    # remove groups from table and reset control
    state$pdata <- state$pdata + 1
  })
}



dsLabelRows <- function(input, output, session, fastq_dir) {

  return(list(
    test = reactive(input$test),
    ctrl = reactive(input$ctrl),
    reset = reactive(input$reset),
    rep = reactive(input$rep),
    pair = reactive(input$pair)
  ))
}



server <- function(input, output, session) {

  # get arguments from calling function
  # defaults for server
  data_dir <- getShinyOption('data_dir', '/srv/shiny-server/drugseqr/scseq/sjia')

  # for testing don't seem to be able to pass arguments as options
  if (isTRUE(getOption('shiny.testmode'))) {

    # reset data for testing
    data_dir <- 'tests/data/test'
    static_dir <- 'tests/data/static'
    unlink(data_dir, recursive = TRUE)
    dir.create(data_dir)
    file.copy(list.files(static_dir, full.names = TRUE), data_dir, recursive = TRUE)
  }


  # single cell analysis and options
  scPage <- callModule(scPage, 'sc',
                       data_dir = data_dir)


  drugsPage <- callModule(drugsPage, 'drug')

  datasetsPage <- callModule(datasetsPage, 'datasets')


}
