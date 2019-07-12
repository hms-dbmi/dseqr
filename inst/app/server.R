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

#' Logic for Drugs page
drugsPage <- function(input, output, session) {

  # the form area inputs/results
  form <- callModule(drugsForm, 'form')

  # the output table
  callModule(drugsTable, 'table',
             cmap_res = form$cmap_res,
             l1000_res = form$l1000_res,
             drug_study = form$drug_study,
             cells = form$cells,
             show_clinical = form$show_clinical)




}

# Logic for form on drugs page
drugsForm <- function(input, output, session) {

  drugStudy <- callModule(selectedDrugStudy, 'drug_study')

  advancedOptions <- callModule(advancedOptions, 'advanced',
                                drugStudy$cmap_res,
                                drugStudy$l1000_res,
                                drugStudy$drug_study,
                                drugStudy$show_advanced)




  return(list(
    cmap_res = drugStudy$cmap_res,
    l1000_res = drugStudy$l1000_res,
    drug_study = drugStudy$drug_study,
    cells = advancedOptions$cells,
    show_clinical = drugStudy$show_clinical
  ))


}

#' Logic for selected drug study
selectedDrugStudy <- function(input, output, session) {

  data_dir <- '/home/alex/Documents/Batcave/zaklab/drugseqr/data-raw/bulk/example-data/IBD'

  cmap_res <- readRDS(file.path(data_dir, 'cmap_res.rds'))
  l1000_res <- readRDS(file.path(data_dir, 'l1000_res.rds'))

  drug_study <- reactive(input$study)

  # boolean for advanced options
  show_advanced <- reactive({
    input$advanced %% 2 != 0
  })


  study_choices <- c('CMAP02', 'L1000')
  updateSelectizeInput(session, 'study', choices = study_choices, selected = 'CMAP02')


  # toggle for clinical status
  show_clinical <- reactive({
    input$clinical %% 2 != 0
  })
  observe({
    toggleClass('advanced', 'btn-primary', condition = show_advanced())
    toggleClass('clinical', 'btn-primary', condition = show_clinical())
  })



  return(list(
    cmap_res = cmap_res,
    l1000_res = l1000_res,
    drug_study = reactive(input$study),
    show_clinical = show_clinical,
    show_advanced = show_advanced
  ))

}

#' Logic for advanced options for selectedDrugStudy
advancedOptions <- function(input, output, session, cmap_res, l1000_res, drug_study, show_advanced) {


  null_cmap <- is.null(cmap_res)
  null_l1000 <- is.null(l1000_res)


  # available cell lines
  if (!null_cmap)
    cmap_cells  <- unique(gsub('^[^_]+_([^_]+)_.+?$', '\\1', names(cmap_res)))

  if (!null_l1000)
    l1000_cells <- unique(gsub('^[^_]+_([^_]+)_.+?$', '\\1', names(l1000_res)))

  #  toggle button styling and showing advanced options
  shiny::observe({
    toggle('advanced-panel', condition = show_advanced(), anim = TRUE)
  })


  # update choices for cell lines
  shiny::observe({
    req(drug_study())
    if (drug_study() == 'L1000') {
      cell_choices <- l1000_cells

    } else if (drug_study() == 'CMAP02') {
      cell_choices <- cmap_cells
    }
    shiny::updateSelectizeInput(session, 'cells', choices = cell_choices, selected = NULL)
  })

  return(list(
    cells = reactive(input$cells)
  ))

}

#' Logic for drug table
drugsTable <- function(input, output, session, cmap_res, l1000_res, drug_study, cells, show_clinical) {

  # generate table to display
  query_res <- shiny::reactive({
    drug_study <- drug_study()

    req(drug_study)
    if (drug_study == 'L1000') {
      query_res <- study_table(l1000_res, 'L1000', cells())

    } else if (drug_study == 'CMAP02') {
      query_res <- study_table(cmap_res, 'CMAP02', cells())
    }

    # for removing entries without a clinical phase
    if (show_clinical()) {
      query_res <- tibble::as_tibble(query_res)
      query_res <- dplyr::filter(query_res, !is.na(.data$`Clinical Phase`))
    }

    return(query_res)
  })

  # show query data
  output$query_res <- DT::renderDataTable({

    query_res <- query_res()
    wide_cols <- c('MOA', 'Target', 'Disease Area', 'Indication', 'Vendor', 'Catalog #', 'Vendor Name')
    # -1 needed with rownames = FALSE
    elipsis_targets <- which(colnames(query_res) %in% wide_cols) - 1

    DT::datatable(
      query_res,
      class = 'cell-border',
      rownames = FALSE,
      selection = 'none',
      escape = FALSE, # to allow HTML in table
      options = list(
        columnDefs = list(list(className = 'dt-nopad sim-cell', height=38, width=120, targets = 0),
                          list(targets = elipsis_targets, render = DT::JS(
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
