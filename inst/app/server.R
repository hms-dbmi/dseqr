dsPage <- function(input, output, session, data_dir) {


  new_anal <- reactiveVal()
  new_dataset <- reactiveVal()

  dsForm <- callModule(dsForm, 'form', data_dir,
                       new_dataset = new_dataset)


  observe({
    toggle('new_table_container', condition = dsForm$show_new())
    toggle('prev_table_container', condition = dsForm$show_prev())
  })

  new_pdata <- callModule(dsNewDatasetTable, 'new',
                          fastq_dir = dsForm$fastq_dir,
                          labels = dsForm$new_labels,
                          paired = dsForm$paired)

  dsPrevTable <- callModule(dsPrevDatasetTable, 'prev',
                            fastq_dir = dsForm$fastq_dir,
                            labels = dsForm$prev_labels)

  # run quantification for new dataset
  observeEvent(dsForm$run_quant(), {
    indices_dir <- '/home/alex/Documents/Batcave/zaklab/drugseqr/data-raw/indices/kallisto'

    # setup
    pdata <- new_pdata()
    paired <- dsForm$paired()
    fastq_dir <- dsForm$fastq_dir()
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
    run_kallisto_bulk(indices_dir = indices_dir,
                      data_dir = fastq_dir,
                      pdata = pdata,
                      paired = paired,
                      updateProgress = updateProgress)

    # generate eset and save
    progress$set(message = 'Annotating dataset')
    eset <- load_seq(fastq_dir)
    new_dataset(dataset_name)

    # save in file to indicate that has been quantified
    quant_path <- file.path(data_dir, 'quantified.rds')
    names(fastq_dir) <- dataset_name

    quant <- readRDS(quant_path)
    saveRDS(c(quant, fastq_dir), quant_path)

    # re-enable inputs
    shinyjs::enable(selector = 'input')
    progress$inc(1)
  })

  observeEvent(dsForm$run_diff(), {

    eset <- dsPrevTable$eset()
    pdata <- dsPrevTable$pdata()
    fastq_dir <- dsForm$fastq_dir()
    anal_name <- dsForm$anal_name()
    dataset_name <- dsForm$dataset_name()
    req(eset, pdata, fastq_dir, anal_name)

    # setup for non-interactive differential expression
    colnames(pdata) <- tolower(colnames(pdata))
    pdata <- data.frame(pdata, row.names = pdata$title)

    diff_expr(eset, data_dir = fastq_dir, anal_name = anal_name, prev_anal = list(pdata = pdata))

    # add to previously analysis
    anal_path <- file.path(fastq_dir, paste0('diff_expr_symbol_', anal_name, '.rds'))
    names(anal_path) <- anal_name

    anals_path <- file.path(data_dir, 'analysed.rds')
    anals <- readRDS(anals_path)
    anals[[dataset_name]] <- c(anals[[dataset_name]], anal_path)
    saveRDS(anals, anals_path)

    new_anal(anal_path)
  })


  return(list(
    new_anal = new_anal,
    data_dir = dsForm$fastq_dir
  ))


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

  new <- callModule(dsFormNew, 'new_dataset',
                    fastq_dir = dataset$fastq_dir)

  prev <- callModule(dsFormPrev, 'prev_dataset')


  return(list(
    fastq_dir = dataset$fastq_dir,
    paired = new$paired,
    new_labels = new$labels,
    prev_labels = prev$labels,
    run_quant = new$run_quant,
    run_diff = prev$run_diff,
    dataset_name = dataset$dataset_name,
    anal_name = prev$anal_name,
    show_new = show_new,
    show_prev = show_prev
  ))

}

dsFormNew <- function(input, output, session, fastq_dir) {

  paired <- callModule(dsEndType, 'end_type',
                       fastq_dir = fastq_dir)

  labels <- callModule(dsLabelNewRows, 'label_rows',
                       paired = paired)

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
    labels = labels,
    run_quant = run_quant
  ))
}


dsFormPrev <- function(input, output, session) {


  run <- reactiveVal()
  labels <- callModule(dsLabelPrevRows, 'label_rows')

  anal_name <- reactive(input$anal_name)

  observeEvent(input$run_diff, {

    if (anal_name() == '') {
      return(NULL)
    }

    run(input$run_diff)
  })

  return(list(
    labels = labels,
    run_diff = run,
    anal_name = anal_name
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

    return(dir)
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


dsNewDatasetTable <- function(input, output, session, fastq_dir, labels, paired) {

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


dsPrevDatasetTable <- function(input, output, session, fastq_dir, labels) {

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

  # reset when new eset
  observe({
    eset <- eset()
    req(eset)

    pdata <- tibble::tibble(Group = NA, Title = colnames(eset))
    group <- rep(NA, nrow(pdata))

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
      options = list(
        columnDefs = list(list(className = 'dt-nopad', targets = 0)),
        scrollY = FALSE,
        paging = FALSE,
        bInfo = 0
      )
    )
  })


  html_pata <- reactive({

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
    req(html_pata())
    DT::replaceData(proxy, html_pata(), rownames = FALSE)
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


dsLabelPrevRows <- function(input, output, session, fastq_dir) {

  return(list(
    test = reactive(input$test),
    ctrl = reactive(input$ctrl),
    reset = reactive(input$reset),
    rep = reactive(input$rep),
    pair = reactive(input$pair)
  ))
}



#' Logic for Drugs page
#' @export
#' @keywords internal
drugsPage <- function(input, output, session, new_anal, bulk_dir) {

  # the form area inputs/results
  form <- callModule(drugsForm, 'form',
                     new_anal = new_anal,
                     bulk_dir = bulk_dir)

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
drugsForm <- function(input, output, session, new_anal, bulk_dir) {

  querySignature <- callModule(querySignature, 'signature',
                               new_anal = new_anal,
                               bulk_dir = bulk_dir)

  drugStudy <- callModule(selectedDrugStudy, 'drug_study',
                          study_choices = querySignature$study_choices)


  advancedOptions <- callModule(advancedOptions, 'advanced',
                                querySignature$cmap_res,
                                querySignature$l1000_res,
                                drugStudy$drug_study,
                                drugStudy$show_advanced)

  query_res <- reactive({
    drug_study <- drugStudy$drug_study()
    cmap_res = querySignature$cmap_res()
    l1000_res = querySignature$l1000_res()
    req(drug_study, cmap_res, l1000_res)

    if (drug_study == 'CMAP02') return(cmap_res)
    if (drug_study == 'L1000') return(l1000_res)

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
selectedDrugStudy <- function(input, output, session, study_choices) {


  drug_study <- reactive(input$study)

  # boolean for advanced options
  show_advanced <- reactive({
    input$advanced %% 2 != 0
  })

  observe({
    req(study_choices())
    updateSelectizeInput(session, 'study', choices = study_choices())
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

querySignature <- function(input, output, session, new_anal, bulk_dir) {
  anals_path <- file.path(bulk_dir, 'analysed.rds')
  if (!file.exists(anals_path)) saveRDS(list(), anals_path)

  cmap_res <- reactiveVal()
  l1000_res <- reactiveVal()
  study_choices <- reactiveVal()
  anals <- reactiveVal()

  # reset choices if new analysis
  observe({
    new_anal()
    anals(readRDS(anals_path))
    study_choices(NULL)
  })


  observe({
    req(anals())
    updateSelectizeInput(session, 'query', choices = anals())
  })

  observe({
    query_path <- input$query
    req(query_path)
    query_dir <- dirname(query_path)

    anal_names <- sapply(anals(), names)
    anal_dirs <- unlist(anals())
    anal_name <- anal_names[which(anal_dirs == query_path)]

    # load or run drug queries
    cmap_res_path <- file.path(query_dir, paste0('cmap_res_', anal_name, '.rds'))
    l1000_res_path <- file.path(query_dir, paste0('l1000_res_', anal_name, '.rds'))

    if (file.exists(cmap_res_path)) {
      cmap_res <- readRDS(cmap_res_path)
      l1000_res <- readRDS(l1000_res_path)

    } else {

      # load drug studies
      cmap_path <- system.file('extdata', 'cmap_es_ind.rds', package = 'drugseqr', mustWork = TRUE)
      cmap_es <- readRDS(cmap_path)

      l1000_path <- system.file('extdata', 'l1000_es.rds', package = 'drugseqr', mustWork = TRUE)
      l1000_es <- readRDS(l1000_path)

      # get dprime effect size values for analysis
      anal <- readRDS(query_path)
      dprimes <- get_dprimes(anal)

      # get correlations between query and drug signatures
      cmap_res <- query_drugs(dprimes, cmap_es)
      l1000_res <- query_drugs(dprimes, l1000_es)

      saveRDS(cmap_res, cmap_res_path)
      saveRDS(l1000_res, l1000_res_path)
    }


    cmap_res(cmap_res)
    l1000_res(l1000_res)
    study_choices(c('CMAP02', 'L1000'))

  })


  return(list(
    cmap_res = cmap_res,
    l1000_res = l1000_res,
    study_choices = study_choices
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

  # will update with proxy to prevent redraw
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
    drug_annot <- drug_annot()
    req(query_res, drug_annot)
    stopifnot(all.equal(drug_annot$title, names(query_res)))

    tibble::add_column(drug_annot, Correlation = query_res, .before=0)
  })

  # subset to selected cells, summarize by compound, and add html
  query_table_summarised <- reactive({
    query_table_full <- query_table_full()
    req(query_table_full)

    query_table <- query_table_full %>%
      limit_cells(cells()) %>%
      summarize_compound() %>%
      add_table_html()
  })

  query_table_final <- reactive({
    query_table <- query_table_summarised()
    req(query_table)

    # subset by clinical phase
    if (show_clinical()) query_table <- dplyr::filter(query_table, !is.na(`Clinical Phase`))

    # sort as desired
    dplyr::arrange(query_table, !!sym(sort_by())) %>%
      select(-min_cor, -avg_cor)
  })


  # show query data
  output$query_table <- DT::renderDataTable({
    # redraw if drug study changes
    drug_study()

    # ellipses for wide columns
    wide_cols <- c('MOA', 'Target', 'Disease Area', 'Indication', 'Vendor', 'Catalog #', 'Vendor Name')
    # -1 needed with rownames = FALSE
    elipsis_targets <- which(colnames(dummy_table) %in% wide_cols) - 1

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
  proxy <- DT::dataTableProxy("query_table")
  observe({
    query_table <- query_table_final()
    req(query_table)
    DT::replaceData(proxy, query_table, rownames = FALSE)
  })
}


server <- function(input, output, session) {

  # get arguments from calling function
  # defaults for server
  # base directory contains data_dir folder
  base_dir <- getShinyOption('base_dir', '/srv/shiny-server/drugseqr/data_dir')

  # for testing don't seem to be able to pass arguments as options
  if (isTRUE(getOption('shiny.testmode'))) {

    # reset data for testing
    data_dir <- 'tests/data/test'
    static_dir <- 'tests/data/static'
    unlink(data_dir, recursive = TRUE)
    dir.create(data_dir)
    file.copy(list.files(static_dir, full.names = TRUE), data_dir, recursive = TRUE)
  }

  # start in base directory so that everything is relative to this locally and on server
  setwd(base_dir)

  sc_dir <- 'single-cell'
  bulk_dir <- 'bulk'

  if (!dir.exists(sc_dir)) dir.create(sc_dir)
  if (!dir.exists(bulk_dir)) dir.create(bulk_dir)


  # single cell analysis and options
  scPage <- callModule(scPage, 'sc',
                       data_dir = sc_dir)

  dsPage <- callModule(dsPage, 'datasets',
                       data_dir = bulk_dir)


  drugsPage <- callModule(drugsPage, 'drug',
                          new_anal = dsPage$new_anal,
                          bulk_dir = bulk_dir)



}
