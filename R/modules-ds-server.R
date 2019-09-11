#' Logic Datasets page
#' @export
#' @keywords internal
dsPage <- function(input, output, session, data_dir) {

  new_anal <- reactiveVal()
  new_dataset <- reactiveVal()
  msg_quant <- reactiveVal()
  msg_anal <- reactiveVal()

  dsForm <- callModule(dsForm, 'form', data_dir,
                       new_dataset = new_dataset,
                       msg_quant = msg_quant,
                       msg_anal = msg_anal,
                       new_anal = new_anal)

  callModule(dsMDSplotly, 'mds_plotly',
             data_dir = data_dir,
             dataset_dir = dsForm$dataset_dir,
             anal_name = dsForm$anal_name,
             new_anal = new_anal)


  observe({
    toggle('quant_table_container', condition = dsForm$show_quant())
    toggle('anal_table_container', condition = dsForm$show_anal())
  })

  dsQuantTable <- callModule(dsQuantTable, 'quant',
                             fastq_dir = dsForm$fastq_dir,
                             labels = dsForm$quant_labels,
                             paired = dsForm$paired)

  dsAnalTable <- callModule(dsAnalTable, 'anal',
                            fastq_dir = dsForm$fastq_dir,
                            labels = dsForm$anal_labels,
                            data_dir = data_dir,
                            dataset_dir = dsForm$dataset_dir,
                            anal_name = dsForm$anal_name)

  observe({
    msg_quant(dsQuantTable$valid_msg())
  })


  observe({
    pdata <- dsAnalTable$pdata()
    valid_msg <- validate_pdata(pdata)
    msg_anal(valid_msg)
  })

  # run quantification for quant dataset
  observeEvent(dsForm$run_quant(), {
    indices_dir <- system.file('indices/kallisto', package = 'drugseqr.data', mustWork = TRUE)

    # setup
    pdata <- dsQuantTable$pdata()
    paired <- dsForm$paired()
    dataset_dir <- dsForm$dataset_dir()
    dataset_name <- dsForm$dataset_name()

    # disable inputs
    shinyjs::disable(selector = 'input')

    # Create a Progress object
    progress <- Progress$new(session, min=0, max = nrow(pdata)+1)
    progress$set(message = "Quantifying files", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())

    # Create a callback function to update progress.
    updateProgress <- function(amount = NULL, detail = NULL) {
      progress$inc(amount = amount, detail = detail)
    }

    # quantification
    fastq_dir <- file.path(data_dir, 'bulk', dataset_dir)
    run_kallisto_bulk(indices_dir = indices_dir,
                      data_dir = fastq_dir,
                      pdata = pdata,
                      paired = paired,
                      updateProgress = updateProgress)

    # generate eset and save
    progress$set(message = 'Annotating dataset')
    eset <- load_seq(fastq_dir)
    new_dataset(dataset_name)

    # save to bulk datasets to indicate that has been quantified
    save_bulk_dataset(dataset_name, dataset_dir, data_dir)

    # re-enable inputs
    shinyjs::enable(selector = 'input')
    progress$inc(1)
  })


  observeEvent(dsForm$run_anal(), {
    # visual that running
    disable(selector = 'input')

    progress <- Progress$new(session, min=0, max = 3)
    on.exit(progress$close())

    # get what need
    eset <- dsAnalTable$eset()
    pdata <- dsAnalTable$pdata()
    fastq_dir <- dsForm$fastq_dir()

    dataset_name <- dsForm$dataset_name()
    dataset_dir <- dsForm$dataset_dir()
    anal_name <- dsForm$anal_name()

    req(eset, pdata, fastq_dir, anal_name, dataset_dir)

    # setup for non-interactive differential expression
    colnames(pdata) <- tolower(colnames(pdata))
    pdata <- pdata[!is.na(pdata$group), ]
    pdata <- data.frame(pdata, row.names = pdata$title)

    # run differential expression
    progress$set(message = "Differential expression", value = 1)
    anal <- diff_expr(eset, data_dir = fastq_dir, anal_name = anal_name, prev_anal = list(pdata = pdata))

    # run pathway analysis
    progress$set(message = "Pathway analysis", value = 2)
    path_anal <- diff_path(eset, prev_anal = anal, data_dir = fastq_dir, anal_name = anal_name, rna_seq = TRUE, NI = 24)

    # add to analysed bulk anals
    save_bulk_anals(dataset_name = dataset_name,
                    dataset_dir = dataset_dir,
                    anal_name = anal_name,
                    data_dir = data_dir)


    # visual that done
    progress$inc(1)
    enable(selector = 'input')
    new_anal(anal_name)
  })


  return(list(
    new_anal = new_anal,
    data_dir = dsForm$fastq_dir
  ))


}

#' Logic for Dataset MDS plotly
#' @export
#' @keywords internal
dsMDSplotly <- function(input, output, session, data_dir, dataset_dir, anal_name, new_anal) {

  # path to saved analysis
  anal_path <- reactive({
    anal_name <- anal_name()
    req(anal_name)
    dataset_dir <- dataset_dir()

    group_file <- paste0('diff_expr_symbol_', anal_name, '.rds')
    file.path(data_dir, 'bulk', dataset_dir, group_file)
  })

  anal <- reactive({
    anal_path <- anal_path()
    req(file.exists(anal_path))
    readRDS(anal_path)
  })


  # MDS plot
  output$plotly <- plotly::renderPlotly({
    # get mds scaling
    mds <- anal()$mds

    # plot
    plotMDS(scaling = mds$scaling, scaling_sva = mds$scaling_sva)
  })


}

#' Logic for Datasets form
#' @export
#' @keywords internal
dsForm <- function(input, output, session, data_dir, new_dataset, msg_quant, msg_anal, new_anal) {

  dataset <- callModule(dsDataset, 'selected_dataset',
                        data_dir = data_dir,
                        new_dataset = new_dataset)


  # show quant, anals or neither
  show_quant <- reactive({
    dataset$dataset_name() != '' && dataset$is.create()
  })

  show_anal <- reactive({
    dataset$dataset_name() != '' && !dataset$is.create()
  })


  observe({
    toggle('quant_dataset_panel', condition = show_quant())
    toggle('anal_dataset_panel', condition = show_anal())
  })

  quant <- callModule(dsFormQuant, 'quant_form',
                      fastq_dir = dataset$fastq_dir,
                      error_msg = msg_quant)

  anal <- callModule(dsFormAnal, 'anal_form',
                     error_msg = msg_anal,
                     data_dir = data_dir,
                     dataset_name = dataset$dataset_name,
                     new_anal = new_anal)


  return(list(
    fastq_dir = dataset$fastq_dir,
    paired = quant$paired,
    quant_labels = quant$labels,
    anal_labels = anal$labels,
    run_quant = quant$run_quant,
    run_anal = anal$run_anal,
    dataset_name = dataset$dataset_name,
    dataset_dir = dataset$dataset_dir,
    anal_name = anal$anal_name,
    show_quant = show_quant,
    show_anal = show_anal
  ))

}

#' Logic for selected dataset part of dsFrom
#' @export
#' @keywords internal
dsDataset <- function(input, output, session, data_dir, new_dataset) {


  # get directory with fastqs
  roots <- c('bulk' = file.path(data_dir, 'bulk'))
  shinyFiles::shinyDirChoose(input, "dataset_dir", roots = roots)


  datasets <- reactive({
    new_dataset()
    load_bulk_datasets(data_dir)
  })

  dataset_name <- reactive(input$dataset_name)

  # is the dataset a quantified one?
  is.create <- reactive({
    dataset_name <- dataset_name()
    datasets <- datasets()
    req(dataset_name)

    !dataset_name %in% datasets$dataset_name
  })

  observe({
    req(datasets())
    updateSelectizeInput(session, 'dataset_name', choices = rbind(rep(NA, 4), datasets()), server = TRUE)
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

  dataset_dir <- reactive({
    dataset_name <- dataset_name()
    req(dataset_name)

    if (is.create()) {
      req(input$dataset_dir)
      dir <- shinyFiles::parseDirPath(roots, input$dataset_dir)
      dir <- basename(as.character(dir))
    }
    else {
      datasets <- datasets()
      req(datasets)
      dir <- datasets[datasets$dataset_name == dataset_name, 'dataset_dir']
    }
    return(dir)
  })


  fastq_dir <- reactive({
    dataset_dir <- dataset_dir()
    req(dataset_dir)
    file.path(data_dir, 'bulk', dataset_dir)
  })


  return(list(
    fastq_dir = fastq_dir,
    dataset_name = dataset_name,
    dataset_dir = dataset_dir,
    is.create = is.create
  ))

}


#' Logic for dataset quantification part of dsForm
#' @export
#' @keywords internal
dsFormQuant <- function(input, output, session, fastq_dir, error_msg) {

  paired <- callModule(dsEndType, 'end_type',
                       fastq_dir = fastq_dir)

  observe(shinyjs::toggleClass("pair", 'disabled', condition = !paired()))


  reset <- reactive(input$reset)
  rep <- reactive(input$rep)
  pair <- reactive(input$pair)


  observe({
    error_msg <- error_msg()
    toggleClass('quant_labels', 'has-error', condition = !is.null(error_msg))
    html('error_msg', html = error_msg)

  })

  quantModal <- function() {
    modalDialog(
      withTags({
        dl(
          dt('End type:'),
          dd('Experiment is selected as single ended'),
          hr(),
          dt('Replicates:'),
          dd('If any - e.g. same sample sequenced in replicate. These will be treated as a single library.')
        )
      }),
      title = 'Double check:',
      size = 's',
      footer = tagList(
        modalButton("Cancel"),
        actionButton(session$ns("confirm"), "Quantify", class = 'pull-left btn-warning')
      )
    )
  }

  # Show modal when button is clicked.
  observeEvent(input$run_quant, {
    showModal(quantModal())
  })

  run_quant <- reactive({
    req(input$confirm)
    removeModal()
    input$confirm
  })

  return(list(
    paired = paired,
    labels = list(
      reset = reset,
      pair = pair,
      rep = rep
    ),
    run_quant = run_quant
  ))
}

#' Logic for end type selection is dsFormQuant
#' @export
#' @keywords internal
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

#' Logic for differential expression analysis part of dsForm
#' @export
#' @keywords internal
dsFormAnal <- function(input, output, session, error_msg, dataset_name, data_dir, new_anal) {


  run_anal <- reactiveVal()
  labels <- list(
    test = reactive(input$test),
    ctrl = reactive(input$ctrl),
    reset = reactive(input$reset)
  )

  anal_name <- reactive(input$anal_name)
  has_anal_name <- reactive(anal_name() != '')



  # analyses (can be multiple) from dataset
  dataset_anals <- reactive({
    # reload if new analysis
    new_anal()
    dataset_name <- dataset_name()
    req(dataset_name)

    anals <- load_bulk_anals(data_dir)
    c('', anals[anals$dataset_name == dataset_name, 'anal_name'])
  })

  observe({
    updateSelectizeInput(session, 'anal_name', choices = dataset_anals())
  })

  is_prev_anal <- reactive({
    anal_name() %in% setdiff(dataset_anals(), '')
  })

  observe({
    toggle('anal_buttons_panel', condition = !is_prev_anal())
  })


  # clear anal name error if type
  observe({
    if (has_anal_name()) {
      removeClass('anal_name_container' ,'has-error')
      html('anal_name_help', '')
    }
  })

  # clear labels error on click
  observe({
    input$test
    input$ctrl
    input$reset
    removeClass('anal_labels', 'has-error')
    html('labels_help', '')
  })

  observeEvent(input$run_anal, {
    # check for analysis name
    if (!has_anal_name()) {
      addClass('anal_name_container', 'has-error')
      html('anal_name_help', 'Requires analysis name.')
      return(NULL)
    }

    # check that pdata is correct
    error_msg <- error_msg()
    if (!is.null(error_msg)) {
      addClass('anal_labels', 'has-error')
      html('labels_help', error_msg)
      return(NULL)
    }

    run_anal(input$run_anal)
  })

  download_content <- reactive({
    anal_name <- anal_name()
    dataset_name <- dataset_name()

    req(anal_name, dataset_name)

    # path to analysis result
    anal_file <- paste0('diff_expr_symbol_', anal_name, '.rds')
    anal_path <- file.path(data_dir, 'bulk', dataset_name, anal_file)
    tt <- readRDS(anal_path)$top_table

    tt[order(tt$P.Value), ]
  })

  output$download <- downloadHandler(
    filename = function() {
      date <- paste0(Sys.Date(), '.csv')
      paste('bulk', dataset_name(), anal_name(), date , sep='_')
    },
    content = function(con) {
      write.csv(download_content(), con)
    }
  )


  return(list(
    labels = labels,
    run_anal = run_anal,
    anal_name = anal_name
  ))
}


#' Logic for dataset quantification table
#' @export
#' @keywords internal
dsQuantTable <- function(input, output, session, fastq_dir, labels, paired) {

  # things user will update and return
  pdata_r <- reactiveVal()
  pairs_r <- reactiveVal()
  reps_r <- reactiveVal()
  valid_msg <- reactiveVal()
  is_rendered <- reactiveVal(FALSE)


  # colors
  background <- 'url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAAPklEQVQoU43Myw0AIAgEUbdAq7VADCQaPyww55dBKyQiHZkzBIwQLqQzCk9E4Ytc6KEPMnTBCG2YIYMVpHAC84EnVbOkv3wAAAAASUVORK5CYII=) repeat'
  group_colors <- c("#C7E9C0", "#C6DBEF", "#FCBBA1", "#FDD0A2", "#BCBDDC", "#D9D9D9", "#F6E8C3", "#DC143C",
                    "#A1D99B", "#9ECAE1", "#FC9272", "#FDAE6B", "#9E9AC8", "#BDBDBD", "#DFC27D", "#FFFFFF",
                    "#C3B091", "#007FFF", "#00FFFF", "#7FFFD4", "#228B22", "#808000", "#7FFF00", "#BFFF00",
                    "#FFD700", "#DAA520", "#FF7F50", "#FA8072","#FC0FC0", "#CC8899", "#E0B0FF", "#B57EDC", "#843179")

  ncolors <- length(group_colors)


  # reset everything when quant fastq_dir
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
    pairs_r(pdata$Pair)
    reps_r(pdata$Replicate)

    pdata$Pair <- pdata$Replicate <- NA
    pdata_r(pdata)
  })

  # redraw table when quant pdata (otherwise update data using proxy)
  output$pdata <- DT::renderDataTable({
    dummy_pdata <- tibble::tibble(Pair = NA, Replicate = NA, 'File Name' = NA)
    is_rendered(TRUE)

    DT::datatable(
      dummy_pdata,
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
    req(is_rendered())
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
    msg <- validate_pairs(pairs, rows, reps)
    valid_msg(msg)

    if (is.null(msg)) {

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
    msg <- validate_reps(pairs, rows, reps)
    valid_msg(msg)

    if (is.null(msg)) {
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

  return(list(
    pdata = returned_pdata,
    valid_msg = valid_msg
  ))
}

#' Logic for differential expression analysis table
#' @export
#' @keywords internal
dsAnalTable <- function(input, output, session, fastq_dir, labels, data_dir, dataset_dir, anal_name) {

  background <- '#337ab7 url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAAPklEQVQoU43Myw0AIAgEUbdAq7VADCQaPyww55dBKyQiHZkzBIwQLqQzCk9E4Ytc6KEPMnTBCG2YIYMVpHAC84EnVbOkv3wAAAAASUVORK5CYII=) repeat'

  # things user will update and return
  pdata_r <- reactiveVal()
  group_r <- reactiveVal()

  eset <- reactive({
    fastq_dir <- fastq_dir()
    eset_path <- file.path(fastq_dir, 'eset.rds')

    req(file.exists(eset_path))
    readRDS(eset_path)
  })

  # pdata that gets returned with group column
  returned_pdata <- reactive({
    pdata <- pdata_r()
    req(pdata)

    pdata$Group <- group_r()
    return(pdata)
  })

  # path to saved analysis
  anal_path <- reactive({
    anal_name <- anal_name()
    dataset_dir <- dataset_dir()
    if (anal_name == '') return(NULL)

    group_file <- paste0('diff_expr_symbol_', anal_name, '.rds')
    file.path(data_dir, 'bulk', dataset_dir, group_file)
  })


  # reset when new eset or analysis name
  observe({
    anal_name()
    eset <- eset()
    anal_path <- anal_path()

    req(eset)

    pdata <- Biobase::pData(eset) %>%
      dplyr::select(-lib.size, -norm.factors) %>%
      tibble::add_column(Group = NA, Title = colnames(eset), .before = 1)


    group <- rep(NA, nrow(pdata))

    # load annotations from previous analysis if available
    if (!is.null(anal_path) && file.exists(anal_path)) {

      anal <- readRDS(anal_path)
      in.anal <- which(row.names(pdata) %in% row.names(anal$pdata))
      group[in.anal] <- anal$pdata$group

      is.test <- which(group == 'test')
      is.ctrl <- which(group == 'ctrl')
      pdata[is.test, 'Group'] <- '<div style="background-color: #e6194b;"><span>test</span></div>'
      pdata[is.ctrl, 'Group'] <- paste0('<div style="background: ', background, ';"><span>control</span></div>')

    }

    pdata_r(pdata)
    group_r(group)
  })


  # redraw table when new pdata (otherwise update data using proxy)
  output$pdata <- DT::renderDataTable({

    DT::datatable(
      pdata_r(),
      class = 'cell-border dt-fake-height',
      rownames = FALSE,
      escape = FALSE, # to allow HTML in table
      extensions = "FixedColumns",
      options = list(
        columnDefs = list(list(className = 'dt-nopad', targets = 0)),
        scrollX = TRUE,
        fixedColumns = TRUE,
        paging = FALSE,
        bInfo = 0
      )
    )
  })


  html_pdata <- reactive({

    # things that trigger update
    pdata <- pdata_r()
    group <- group_r()
    req(pdata)

    is.test <- which(group == 'test')
    is.ctrl <- which(group == 'ctrl')
    pdata[is.test, 'Group'] <- '<div style="background-color: #e6194b;"><span>test</span></div>'
    pdata[is.ctrl, 'Group'] <- paste0('<div style="background: ', background, ';"><span>control</span></div>')


    return(pdata)
  })


  # proxy used to replace data
  proxy <- DT::dataTableProxy("pdata")
  shiny::observe({
    req(html_pdata())
    DT::replaceData(proxy, html_pdata(), rownames = FALSE)
  })



  # click 'Test'
  shiny::observeEvent(labels$test(), {
    test <- labels$test()
    group <- group_r()
    req(test)

    group[input$pdata_rows_selected] <- 'test'
    group_r(group)
  })

  # click 'Control'
  shiny::observeEvent(labels$ctrl(), {
    ctrl <- labels$ctrl()
    group <- group_r()
    req(ctrl)

    group[input$pdata_rows_selected] <- 'ctrl'
    group_r(group)

  })


  # click 'Reset'
  shiny::observeEvent(labels$reset(), {
    req(labels$reset())
    group <- group_r()
    group_r(rep(NA, length(group)))

  })


  return(list(
    pdata = returned_pdata,
    eset = eset
  ))

}
