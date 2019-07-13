dsPage <- function(input, output, session, data_dir) {


  new_dataset <- reactiveVal()

  dsForm <- callModule(dsForm, 'form', data_dir,
                       new_dataset = new_dataset)


  observe({
    toggle('new_table_container', condition = dsForm$show_new())
    toggle('prev_table_container', condition = dsForm$show_prev())
  })

  new_pdata <- callModule(dsTable, 'new',
                          fastq_dir = dsForm$fastq_dir,
                          labels = dsForm$new_labels,
                          paired = dsForm$paired)

  # dsPrevTable <- callModule(dsTable, 'prev',
  #                           fastq_dir = dsForm$fastq_dir,
  #                           labels = dsForm$new_inputs$labels,
  #                           paired = dsForm$new_inputs$paired)

  # run quantification for new dataset
  observeEvent(dsForm$run_quant(), {
    indices_dir <- '/home/alex/Documents/Batcave/zaklab/drugseqr/data-raw/indices/kallisto'

    # setup
    pdata <- new_pdata()
    paired <- dsForm$paired()
    fastq_dir <- dsForm$fastq_dir()

    shinyjs::disable(selector = 'input')

    # quantification
    run_kallisto_bulk(indices_dir = indices_dir,
                      data_dir = fastq_dir,
                      pdata = pdata,
                      paired = paired)

    shinyjs::enable(selector = 'input')

    # save in file to indicate that has been quantified
    quant_path <- file.path(data_dir, 'quantified.rds')

    dataset_name <- dsForm$dataset_name()
    names(fastq_dir) <- dataset_name

    quant <- readRDS(quant_path)
    saveRDS(c(quant, fastq_dir), quant_path)

    new_dataset(dataset_name)
  })

  # run differential expression analysis for previous dataset


}


dsForm <- function(input, output, session, data_dir, new_dataset) {

  dataset <- callModule(dsSelectedDataset, 'selected_dataset',
                        data_dir = data_dir,
                        new_dataset = new_dataset)


  # show new, previous or neither
  show_new <- reactive({
    dataset$dataset_name() != '' && dataset$is.create()
  })
  show_prev <- reactive({
    dataset$dataset_name() != '' && !dataset$is.create()
  })


  observe({
    toggle('new_dataset_panel', condition = show_new())
    toggle('prev_dataset_panel', condition = show_prev())
  })

  new_inputs <- callModule(dsNewInputs, 'new_dataset',
                           fastq_dir = dataset$fastq_dir)


  return(list(
    fastq_dir = dataset$fastq_dir,
    paired = new_inputs$paired,
    new_labels = new_inputs$labels,
    run_quant = new_inputs$run_quant,
    dataset_name = dataset$dataset_name,
    show_new = show_new,
    show_prev = show_prev
  ))

}

dsNewInputs <- function(input, output, session, fastq_dir) {

  paired <- callModule(dsEndType, 'end_type',
                       fastq_dir = fastq_dir)

  labels <- callModule(dsLabelNewRows, 'label_rows',
                       paired = paired)

  quantModal <- modalDialog(
    withTags({
      dl(
        dt('End type'),
        dd('Experiment is selected as single ended'),
        hr(),
        dt('Replicates'),
        dd('If any - e.g. same sample sequenced in replicate. These will be treated as a single library.')
      )
    }),
    title = 'Double check:',
    size = 's',
    footer = tagList(
      modalButton("Cancel"),
      actionButton("confirm", "Quantify", class = 'pull-left btn-danger')
    )
  )

  # Show modal when button is clicked.
  observeEvent(input$run_quant, {
    showModal(quantModal)
  })

  run_quant <- reactive({
    req(input$confirm)
    input$confirm
  })

  return(list(
    paired = paired,
    labels = labels,
    run_quant = run_quant
  ))
}

dsLabelNewRows <- function(input, output, session, paired) {

  observe(
    shinyjs::toggleClass("pair", 'disabled', condition = !paired())
  )

  return(list(
    reset = reactive(input$reset),
    rep = reactive(input$rep),
    pair = reactive(input$pair)
  ))
}

dsSelectedDataset <- function(input, output, session, data_dir, new_dataset) {

  quant_path <- file.path(data_dir, 'quantified.rds')
  if (!file.exists(quant_path)) saveRDS('', quant_path)

  # dataset names/dirs stored as names/values in quantified.rds
  prev_paths <- reactiveVal(readRDS(quant_path))
  prev_datasets <- reactive(names(prev_paths()))
  is.create <- reactiveVal(FALSE)

  # get directory with fastqs
  bulk_dir <- c('bulk' = data_dir)
  shinyFiles::shinyDirChoose(input, "dataset_dir", roots = bulk_dir)

  # update prev_paths/datasets if new one quantified
  observe({
    req(new_dataset())
    prev_paths(readRDS(quant_path))
  })

  observe({
    req(prev_datasets())
    updateSelectizeInput(session, 'dataset_name', choices = prev_datasets())
  })


  # is the dataset a new one?
  observe({
    req(input$dataset_name)
    # for creating new
    is.create(!input$dataset_name %in% prev_datasets())
  })

  # add disabled class if not creating
  observe({
    toggleClass('dataset_dir', 'disabled', condition = !is.create())
  })


  # open selector if creating
  observeEvent(is.create(), {
    req(is.create())
    if (is.create()) {
      shinyjs::click('dataset_dir')
    }
  })


  fastq_dir <- reactive({
    req(input$dataset_name)

    if (is.create()) {
      req(input$dataset_dir)
      dir <- shinyFiles::parseDirPath(bulk_dir, input$dataset_dir)
    }
    else {
      ds <- prev_datasets()
      req(ds)
      dir <- prev_paths()[which(ds == input$dataset_name)]
    }

    normalizePath(dir)
  })


  dataset_name <- reactive(input$dataset_name)

  return(list(
    dataset_name = dataset_name,
    fastq_dir = fastq_dir,
    is.create = is.create
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

  return(paired = reactive(input$end_type == 'pair-ended'))
}


dsTable <- function(input, output, session, fastq_dir, labels, paired) {

  # things user will update and return
  pdata_r <- reactiveVal()
  pairs_r <- reactiveVal()
  reps_r <- reactiveVal()


  # colors
  background <- 'url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAAPklEQVQoU43Myw0AIAgEUbdAq7VADCQaPyww55dBKyQiHZkzBIwQLqQzCk9E4Ytc6KEPMnTBCG2YIYMVpHAC84EnVbOkv3wAAAAASUVORK5CYII=) repeat'
  group_colors <- c("#C7E9C0", "#C6DBEF", "#FCBBA1", "#FDD0A2", "#BCBDDC", "#D9D9D9", "#F6E8C3", "#DC143C",
                    "#A1D99B", "#9ECAE1", "#FC9272", "#FDAE6B", "#9E9AC8", "#BDBDBD", "#DFC27D", "#FFFFFF",
                    "#C3B091", "#007FFF", "#00FFFF", "#7FFFD4", "#228B22", "#808000", "#7FFF00", "#BFFF00",
                    "#FFD700", "#DAA520", "#FF7F50", "#FA8072","#FC0FC0", "#CC8899", "#E0B0FF", "#B57EDC", "#843179")

  ncolors <- length(group_colors)


  # reset everything when new fastq_dir
  observeEvent(fastq_dir(), {
    fastq_dir <- fastq_dir()
    req(fastq_dir)

    pdata_path <- file.path(fastq_dir, 'pdata.rds')

    # initial creation of saved pdata
    if (!file.exists(pdata_path)) {
      fastqs <- list.files(fastq_dir, '.fastq.gz$')
      pdata <- tibble::tibble('File Name' = fastqs)
      pdata <- tibble::add_column(pdata, Pair = NA, Replicate = NA, .before = 1)
      saveRDS(pdata, pdata_path)
    }

    pdata <- readRDS(pdata_path)

    pdata_r(pdata)
    pairs_r(pdata$Pair)
    reps_r(pdata$Replicate)
  })

  # redraw table when new pdata (otherwise update data using proxy)
  output$pdata <- DT::renderDataTable({

    DT::datatable(
      pdata_r(),
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


  # pdata that gets returned with reps and pairs
  returned_pdata <- reactive({
    pdata <- pdata_r()
    pdata$Replicate <- reps_r()
    pdata$Pair <- pairs_r()

    return(pdata)
  })


  # save returned pdata so that don't lose work
  observe({
    pdata <- returned_pdata()
    req(pdata)

    pdata_path <- file.path(isolate(fastq_dir()), 'pdata.rds')
    saveRDS(pdata, pdata_path)
  })


  # pdata to update Pair/Replicate column in proxy (uses html)
  html_pata <- reactive({

    # things that trigger update
    pdata <- pdata_r()
    req(pdata)
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
    if (paired()) {
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
    DT::replaceData(proxy, html_pata(), rownames = FALSE)
  })


  # click 'Paired'
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
    }
  })

  # click 'Replicate'
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
    }
  })


  # click 'Reset'
  shiny::observeEvent(labels$reset(), {
    pdata <- pdata_r()
    clear <- rep(NA, nrow(pdata))
    reps_r(clear)
    pairs_r(clear)
  })

  return(returned_pdata)
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


  sc_dir <- file.path(data_dir, 'single-cell')
  bulk_dir <- file.path(data_dir, 'bulk')

  # single cell analysis and options
  scPage <- callModule(scPage, 'sc',
                       data_dir = sc_dir)


  drugsPage <- callModule(drugsPage, 'drug')

  datasetsPage <- callModule(dsPage, 'datasets',
                             data_dir = bulk_dir)


}
