dsPage <- function(input, output, session, data_dir) {


  new_anal <- reactiveVal()
  quant_dataset <- reactiveVal()
  msg_quant <- reactiveVal()
  msg_anal <- reactiveVal()
  is_prev_anal <- reactiveVal()

  dsForm <- callModule(dsForm, 'form', data_dir,
                       quant_dataset = quant_dataset,
                       msg_quant = msg_quant,
                       msg_anal = msg_anal,
                       is_prev_anal = is_prev_anal)


  observe({
    toggle('quant_table_container', condition = dsForm$show_quant())
    toggle('anal_table_container', condition = dsForm$show_anal())
  })

  dsQuantTable <- callModule(dsQuantDatasetTable, 'quant',
                             fastq_dir = dsForm$fastq_dir,
                             labels = dsForm$quant_labels,
                             paired = dsForm$paired)

  dsAnalTable <- callModule(dsAnalDatasetTable, 'anal',
                            fastq_dir = dsForm$fastq_dir,
                            labels = dsForm$anal_labels,
                            data_dir = data_dir,
                            dataset_name = dsForm$dataset_name,
                            anal_name = dsForm$anal_name)

  observe({
    msg_quant(dsQuantTable$valid_msg())
  })

  observe({
    is_prev_anal(dsAnalTable$is_prev_anal())
  })

  observe({
    pdata <- dsAnalTable$pdata()
    valid_msg <- validate_pdata(pdata)
    msg_anal(valid_msg)
  })

  # run quantification for quant dataset
  observeEvent(dsForm$run_quant(), {
    #TODO make not hardcoded
    indices_dir <- '/home/alex/Documents/Batcave/zaklab/drugseqr/data-raw/indices/kallisto'


    # setup
    pdata <- dsQuantTable$pdata()
    paired <- dsForm$paired()
    dataset_dir <- dsForm$dataset_dir()
    dataset_name <- dsForm$dataset_name()

    # disable inputs
    shinyjs::disable(selector = 'input')

    # Create a Progress object
    progress <- Progress$quant(session, min=0, max = nrow(pdata)+1)
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
    quant_dataset(dataset_name)

    # save to bulk datasets to indicate that has been quantified
    add_bulk_dataset(dataset_name, dataset_dir, data_dir)

    # re-enable inputs
    shinyjs::enable(selector = 'input')
    progress$inc(1)
  })


  observeEvent(dsForm$run_anal(), {
    # visual that running
    disable(selector = 'input')

    progress <- Progress$new(session, min=0, max = 2)
    progress$set(message = "Differential expression", value = 1)
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

    # run
    diff_expr(eset, data_dir = fastq_dir, anal_name = anal_name, prev_anal = list(pdata = pdata))

    # add to analious anals
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

validate_pdata <- function(pdata) {
  group <- pdata$Group
  group <- group[!is.na(group)]

  if (length(unique(group)) != 2) {
    msg <- 'Analysis requires test and control groups'

  } else if (length(group) < 3) {
    msg <- 'At least three samples are required for analysis'

  } else {
    msg <- NULL
  }
  return(msg)
}


save_bulk_anals <- function(dataset_name, dataset_dir, anal_name, data_dir) {
  anals_path <- file.path(data_dir, 'bulk', 'anals.rds')
  anals <- readRDS(anals_path)

  anals[nrow(anals)+1, ] <- c(dataset_name, dataset_dir, anal_name)
  saveRDS(anals, anals_path)
}


dsForm <- function(input, output, session, data_dir, quant_dataset, msg_quant, msg_anal, is_prev_anal) {

  dataset <- callModule(dsSelectedDataset, 'selected_dataset',
                        data_dir = data_dir,
                        quant_dataset = quant_dataset)


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
                     is_prev_anal = is_prev_anal)


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


dsFormAnal <- function(input, output, session, error_msg, dataset_name, data_dir, is_prev_anal) {


  run_anal <- reactiveVal()
  labels <- list(
    test = reactive(input$test),
    ctrl = reactive(input$ctrl),
    reset = reactive(input$reset)
  )

  observe({
    toggle('anal_buttons_panel', condition = is_prev_anal())
  })

  # analyses (can be multiple) from dataset
  dataset_anals <- reactive({
    dataset_name <- dataset_name()
    req(dataset_name)

    anals <- load_bulk_anals(data_dir)
    c('', anals[anals$dataset_name == dataset_name, 'anal_name'])
  })

  observe({
    updateSelectizeInput(session, 'anal_name', choices = dataset_anals())
  })

  anal_name <- reactive(input$anal_name)

  has_anal_name <- reactive(anal_name() != '')

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

  return(list(
    labels = labels,
    run_anal = run_anal,
    anal_name = anal_name
  ))
}


save_bulk_dataset <- function(dataset_name, dataset_dir, data_dir) {
  datasets_path <- file.path(data_dir, 'bulk', 'datasets.rds')
  datasets <- readRDS(datasets_path)

  datasets[nrow(datasets)+1, ] <- c(dataset_name, dataset_dir)
  saveRDS(datasets, datasets_path)
}

load_bulk_datasets <-function(data_dir) {
  datasets_path <- file.path(data_dir, 'bulk', 'datasets.rds')

  if (file.exists(datasets_path)) {
    datasets <- readRDS(datasets_path)

  } else {
    datasets <- data.frame(matrix(ncol = 2, nrow = 0), stringsAsFactors = FALSE)
    colnames(datasets) <- c("dataset_name", "dataset_dir")
    saveRDS(datasets, datasets_path)
  }

  datasets$value <-  datasets$label <- datasets$dataset_name

  return(datasets)
}



dsSelectedDataset <- function(input, output, session, data_dir, quant_dataset) {


  # get directory with fastqs
  roots <- c('bulk' = file.path(data_dir, 'bulk'))
  shinyFiles::shinyDirChoose(input, "dataset_dir", roots = roots)


  datasets <- reactive({
    quant_dataset()
    load_bulk_datasets(data_dir)
  })

  dataset_name <- reactive(input$dataset_name)

  # is the dataset a quant one?
  is.create <- reactive({
    dataset_name <- dataset_name()
    datasets <- datasets()
    req(dataset_name)

    !dataset_name %in% datasets$dataset_name
  })

  observe({
    req(datasets())
    updateSelectizeInput(session, 'dataset_name', choices = datasets(), server = TRUE)
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


dsQuantDatasetTable <- function(input, output, session, fastq_dir, labels, paired) {

  # things user will update and return
  pdata_r <- reactiveVal()
  pairs_r <- reactiveVal()
  reps_r <- reactiveVal()
  valid_msg <- reactiveVal()


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

    pdata_r(pdata)
    pairs_r(pdata$Pair)
    reps_r(pdata$Replicate)
  })

  # redraw table when quant pdata (otherwise update data using proxy)
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


dsAnalDatasetTable <- function(input, output, session, fastq_dir, labels, data_dir, dataset_name, anal_name) {

  background <- '#e9305d url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAAPklEQVQoU43Myw0AIAgEUbdAq7VADCQaPyww55dBKyQiHZkzBIwQLqQzCk9E4Ytc6KEPMnTBCG2YIYMVpHAC84EnVbOkv3wAAAAASUVORK5CYII=) repeat'

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

    pdata$Group <- group_r()

    return(pdata)
  })

  # path to previous analysis
  anal_path <- reactive({
    anal_name <- anal_name()
    anals <- load_bulk_anals(data_dir)
    dataset_dir <- anals$dataset_dir[anals$dataset_name == dataset_name() & anals$anal_name == anal_name()]
    if (!length(dataset_dir) | anal_name == '') return(NULL)

    anal_file <- paste0('diff_expr_symbol_', anal_name(), '.rds')
    file.path(data_dir, 'bulk', dataset_dir, anal_file)
  })

  is_prev_anal <- reactive(is.null(anal_path()))

  # pdata from previous analysis
  analysed_pdata <- reactive({
    anal_path <- anal_path()
    if (is.null(anal_path)) return(NULL)

    anal <- readRDS(anal_path)
    anal$pdata
  })


  # reset when new eset or analysis name
  observe({
    eset <- eset()
    req(eset)

    pdata <- data.frame(Group = NA, Title = colnames(eset), row.names = colnames(eset), stringsAsFactors = FALSE)
    group <- rep(NA, nrow(pdata))

    pdata_r(pdata)
    group_r(group)
  })

  # update groups if have analysed pdata
  observe({
    anal_pdata <- analysed_pdata()
    pdata <- pdata_r()
    pdata[row.names(anal_pdata), 'Group'] <- anal_pdata$group
    group <- pdata$Group
    group_r(group)
  })

  # redraw table when new pdata (otherwise update data using proxy)
  output$pdata <- DT::renderDataTable({
    dummy_pdata <- data.frame(Group = NA, Title = NA, stringsAsFactors = FALSE)

    DT::datatable(
      pdata_r(),
      class = 'cell-border dt-fake-height',
      rownames = FALSE,
      escape = FALSE, # to allow HTML in table
      options = list(
        columnDefs = list(list(className = 'dt-nopad', targets = 0)),
        scrollY = FALSE,
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
    eset = eset,
    is_prev_anal = is_prev_anal
  ))

}


dsLabelAnalRows <- function(input, output, session, fastq_dir) {


}



#' Logic for Drugs page
#' @export
#' @keywords internal
drugsPage <- function(input, output, session, new_anal, data_dir) {

  # the form area inputs/results
  form <- callModule(drugsForm, 'form',
                     new_anal = new_anal,
                     data_dir = data_dir)

  # the output table
  callModule(drugsTable, 'table',
             query_res = form$query_res,
             drug_study = form$drug_study,
             cells = form$cells,
             sort_by = form$sort_by,
             show_clinical = form$show_clinical)




}

# Logic for form on drugs page
#' @export
#' @keywords internal
drugsForm <- function(input, output, session, new_anal, data_dir) {

  querySignature <- callModule(querySignature, 'signature',
                               new_anal = new_anal,
                               data_dir = data_dir)

  drugStudy <- callModule(selectedDrugStudy, 'drug_study',
                          anal = querySignature$anal)

  advancedOptions <- callModule(advancedOptions, 'advanced',
                                cmap_res = querySignature$cmap_res,
                                l1000_res = querySignature$l1000_res,
                                drug_study = drugStudy$drug_study,
                                show_advanced = drugStudy$show_advanced)

  query_res <- reactive({
    drug_study <- drugStudy$drug_study()
    cmap_res <- querySignature$cmap_res()
    l1000_res <- querySignature$l1000_res()

    if (drug_study == 'CMAP02') return(cmap_res)
    if (drug_study == 'L1000') return(l1000_res)

    return(NULL)
  })




  return(list(
    query_res = query_res,
    drug_study = drugStudy$drug_study,
    cells = advancedOptions$cells,
    sort_by = advancedOptions$sort_by,
    show_clinical = drugStudy$show_clinical
  ))


}

#' Logic for selected drug study
#' @export
#' @keywords internal
selectedDrugStudy <- function(input, output, session, anal) {


  drug_study <- reactive(input$study)

  # boolean for advanced options
  show_advanced <- reactive({
    input$advanced %% 2 != 0
  })

  observe({
    req(anal())
    updateSelectizeInput(session, 'study', choices = c('CMAP02', 'L1000'), selected = NULL)
  })

  # toggle for clinical status
  show_clinical <- reactive({
    input$clinical %% 2 != 0
  })
  observe({
    toggleClass('advanced', 'btn-primary', condition = show_advanced())
    toggleClass('clinical', 'btn-primary', condition = show_clinical())
  })



  return(list(
    drug_study = drug_study,
    show_clinical = show_clinical,
    show_advanced = show_advanced
  ))

}

load_bulk_anals <- function(data_dir) {
  anals_path <- file.path(data_dir, 'bulk', 'anals.rds')

  if (file.exists(anals_path)) {
    anals <- readRDS(anals_path)

  } else {
    anals <- data.frame(matrix(ncol = 3, nrow = 0), stringsAsFactors = FALSE)
    colnames(anals) <- c("dataset_name", "dataset_dir", "anal_name")
    saveRDS(anals, anals_path)
  }

  anals$label <- anals$anal_name
  anals$value <- 1:nrow(anals)

  return(anals)

}

querySignature <- function(input, output, session, new_anal, data_dir) {


  cmap_res <- reactiveVal()
  l1000_res <- reactiveVal()

  # reload query choices if quant analysis
  anals <- reactive({
    new_anal()
    load_bulk_anals(data_dir)
  })

  observe({
    anals <- anals()
    req(anals)
    updateSelectizeInput(session, 'query', choices = anals, server = TRUE)
  })

  anal <- reactive({
    row_num <- input$query
    anals <- anals()
    req(row_num, anals)

    anals[row_num, ]
  })


  # paths to analysis and drug query results
  res_paths <- reactive({
    anal <- anal()
    dataset_dir <- file.path(data_dir, 'bulk', anal$dataset_dir)
    anal_name <- anal$anal_name

    list(
      anal = file.path(dataset_dir, paste0('diff_expr_symbol_', anal_name, '.rds')),
      cmap = file.path(dataset_dir, paste0('cmap_res_', anal_name, '.rds')),
      l1000 = file.path(dataset_dir, paste0('l1000_res_', anal_name, '.rds'))
    )
  })


  # get saved cmap/l1000 query results
  observe({
    res_paths <- res_paths()

    # load if available
    if (file.exists(res_paths$cmap)) {
      cmap_res <- readRDS(res_paths$cmap)
      l1000_res <- readRDS(res_paths$l1000)


    } else {
      # otherwise run
      disable('query')

      progress <- Progress$new(session, min = 0, max = 4)
      progress$set(message = "Querying drugs", value = 1)
      on.exit(progress$close())

      cmap_path  <- system.file('extdata', 'cmap_es_ind.rds', package = 'drugseqr', mustWork = TRUE)
      l1000_path <- system.file('extdata', 'l1000_es.rds', package = 'drugseqr', mustWork = TRUE)

      cmap_es  <- readRDS(cmap_path)
      progress$inc(1)
      l1000_es <- readRDS(l1000_path)
      progress$inc(1)


      # get dprime effect size values for analysis
      anal <- readRDS(res_paths$anal)
      dprimes <- get_dprimes(anal)

      # get correlations between query and drug signatures
      cmap_res <- query_drugs(dprimes, cmap_es)
      l1000_res <- query_drugs(dprimes, l1000_es)
      progress$inc(1)

      saveRDS(cmap_res, res_paths$cmap)
      saveRDS(l1000_res, res_paths$l1000)
      enable('query')
    }

    cmap_res(cmap_res)
    l1000_res(l1000_res)
  })



  return(list(
    cmap_res = cmap_res,
    l1000_res = l1000_res,
    anal = anal
  ))

}

#' Logic for advanced options for selectedDrugStudy
#' @export
#' @keywords internal
advancedOptions <- function(input, output, session, cmap_res, l1000_res, drug_study, show_advanced) {

  # available cell lines
  cmap_cells <- unique(cmap_annot$cell)
  l1000_cells <- unique(l1000_annot$cell)

  # update choices for cell lines based on selected study
  cell_choices <- shiny::reactive({
    req(drug_study())
    if (drug_study() == 'L1000') return(l1000_cells)
    else if (drug_study() == 'CMAP02') return(cmap_cells)
  })

  #  toggle  showing advanced options
  shiny::observe({
    toggle('advanced-panel', condition = show_advanced(), anim = TRUE)
  })

  # update choices for cell lines
  shiny::observe({
    shiny::updateSelectizeInput(session, 'cells', choices = cell_choices(), selected = NULL)
  })

  return(list(
    cells = reactive(input$cells),
    sort_by = reactive(input$sort_by)
  ))

}

#' Logic for drug table
#' @export
#' @keywords internal
#' @importFrom magrittr "%>%"
drugsTable <- function(input, output, session, query_res, drug_study, cells, show_clinical, sort_by) {

  # will update with proxy to analent redraw
  dummy_table <- data.frame('Correlation' = NA,
                            'Compound' = NA,
                            'Clinical Phase' = NA,
                            'External Links' = NA,
                            'MOA' = NA,
                            'Target' = NA,
                            'Disease Area' = NA,
                            'Indication' = NA,
                            'Vendor' = NA,
                            'Catalog #' = NA,
                            'Vendor Name' = NA, check.names = FALSE)
  dummy_rendered <- reactiveVal(FALSE)

  # get either cmap or l1000 annotations
  drug_annot <- reactive({
    drug_study <- drug_study()
    req(drug_study)

    if (drug_study == 'CMAP02') return(cmap_annot)
    else if (drug_study == 'L1000') return(l1000_annot)
  })

  # add annotations to query result
  query_table_full <- reactive({
    query_res <- query_res()
    if (is.null(query_res)) return(NULL)

    drug_annot <- drug_annot()
    req(query_res, drug_annot)
    stopifnot(all.equal(drug_annot$title, names(query_res)))

    tibble::add_column(drug_annot, Correlation = query_res, .before=0)
  })

  # subset to selected cells, summarize by compound, and add html
  query_table_summarised <- reactive({
    query_table_full <- query_table_full()
    if (is.null(query_table_full)) return(NULL)

    query_table <- query_table_full %>%
      limit_cells(cells()) %>%
      summarize_compound() %>%
      add_table_html()
  })

  query_table_final <- reactive({
    query_table <- query_table_summarised()
    if (is.null(query_table)) return(NULL)
    sort_by <- sort_by()

    # subset by clinical phase
    if (show_clinical()) query_table <- dplyr::filter(query_table, !is.na(`Clinical Phase`))

    if (sort_by == 'avg_cor') {
      query_table$Correlation <- gsub('simplot', 'simplot show-meanline', query_table$Correlation)
    }

    # sort as desired
    dplyr::arrange(query_table, !!sym(sort_by)) %>%
      select(-min_cor, -avg_cor)
  })


  # show query data
  output$query_table <- DT::renderDataTable({
    # ellipses for wide columns
    wide_cols <- c('MOA', 'Target', 'Disease Area', 'Indication', 'Vendor', 'Catalog #', 'Vendor Name')
    # -1 needed with rownames = FALSE
    elipsis_targets <- which(colnames(dummy_table) %in% wide_cols) - 1
    dummy_rendered(TRUE)

    DT::datatable(
      dummy_table,
      class = 'cell-border',
      rownames = FALSE,
      selection = 'none',
      escape = FALSE, # to allow HTML in table
      options = list(
        columnDefs = list(list(className = 'dt-nopad sim-cell', height=38, width=120, targets = 0),
                          list(targets = c(4,5,6,7), render = DT::JS(
                            "function(data, type, row, meta) {",
                            "return type === 'display' && data !== null && data.length > 17 ?",
                            "'<span title=\"' + data + '\">' + data.substr(0, 17) + '...</span>' : data;",
                            "}"))),
        ordering=FALSE,
        scrollX = TRUE,
        pageLength = 20,
        paging = TRUE,
        bInfo = 0,
        dom = 'ftp'
      )
    )
  }, server = TRUE)

  # proxy used to replace data
  # low priority to make sure data has been rendered
  proxy <- DT::dataTableProxy("query_table")
  observe({
    req(dummy_rendered())
    query_table <- query_table_final()
    DT::replaceData(proxy, query_table, rownames = FALSE)
  })
}


server <- function(input, output, session) {

  # get arguments from calling function
  # defaults for server
  # base directory contains data_dir folder
  data_dir <- getShinyOption('data_dir', '/srv/shiny-server/drugseqr/data_dir')

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

  dir.create(sc_dir, showWarnings = FALSE)
  dir.create(bulk_dir, showWarnings = FALSE)


  # single cell analysis and options
  scPage <- callModule(scPage, 'sc',
                       sc_dir = sc_dir)

  dsPage <- callModule(dsPage, 'datasets',
                       data_dir = data_dir)


  drugsPage <- callModule(drugsPage, 'drug',
                          new_anal = dsPage$new_anal,
                          data_dir = data_dir)



}
