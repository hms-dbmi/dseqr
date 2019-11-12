#' Logic Datasets page
#' @export
#' @keywords internal
dsPage <- function(input, output, session, data_dir, sc_dir, indices_dir) {

  new_anal <- reactiveVal()
  new_dataset <- reactiveVal()
  msg_quant <- reactiveVal()
  msg_anal <- reactiveVal()


  dsForm <- callModule(dsForm, 'form',
                       data_dir = data_dir,
                       sc_dir = sc_dir,
                       new_dataset = new_dataset,
                       msg_quant = msg_quant,
                       msg_anal = msg_anal,
                       new_anal = new_anal)



  callModule(dsMDSplotly, 'mds_plotly',
             data_dir = data_dir,
             dataset_dir = dsForm$dataset_dir,
             anal_name = dsForm$anal_name,
             new_anal = new_anal)


  # toggle tables
  observe({
    toggle('quant_table_container', condition = dsForm$show_quant())
    toggle('anal_table_container', condition = dsForm$show_anal() & !dsForm$is.explore())
    toggle('explore_table_container', condition = dsForm$is.explore())
  })

  # toggle plots
  observe({
    toggle('mds_plotly_container', condition = dsForm$show_anal() & !dsForm$is.explore())
    toggle('gene_plotly_container', condition = dsForm$is.explore())
  })

  dsQuantTable <- callModule(dsQuantTable, 'quant',
                             fastq_dir = dsForm$fastq_dir,
                             labels = dsForm$quant_labels,
                             paired = dsForm$paired)

  dsAnalTable <- callModule(dsAnalTable, 'anal',
                            eset = dsForm$eset,
                            labels = dsForm$anal_labels,
                            data_dir = data_dir,
                            dataset_dir = dsForm$dataset_dir,
                            anal_name = dsForm$anal_name)

  dsExploreTable <- callModule(dsExploreTable, 'explore',
                               eset = dsForm$eset,
                               labels = dsForm$explore_labels,
                               data_dir = data_dir,
                               dataset_dir = dsForm$dataset_dir)

  callModule(dsGenePlotly, 'gene_plotly',
             eset = dsForm$eset,
             explore_genes = dsForm$explore_genes,
             pdata = dsExploreTable$pdata,
             dataset_name = dsForm$dataset_name)

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
    # disable inputs
    shinyjs::disable(selector = 'input')

    # setup
    is.sc <- dsForm$is.sc()
    is.cellranger <- dsForm$is.cellranger()

    pdata <- dsQuantTable$pdata()
    dataset_dir <- dsForm$dataset_dir()
    dataset_name <- dsForm$dataset_name()
    fastq_dir <- file.path(data_dir, dataset_dir)


    # Create a Progress object
    progress <- Progress$new(session, min=0, max = ifelse(is.sc, 8, nrow(pdata)+1))
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())

    if (is.sc) {
      progress$set(message = "Quantifying files", value = 0)
      progress$set(value = 1)
      if (!is.cellranger) run_kallisto_scseq(indices_dir, fastq_dir)

      progress$set(message = "Loading and QC", value = 2)
      type <- ifelse(is.cellranger, 'cellranger', 'kallisto')
      scseq <- load_scseq(fastq_dir, project = dataset_name, type = type)
      scseq <- scseq[, scseq$whitelist]
      gc()

      progress$set(message = "Preprocessing", value = 3)
      scseq <- preprocess_scseq(scseq)
      gc()

      progress$set(message = "Clustering", value = 4)
      scseq <- add_scseq_clusters(scseq)
      gc()

      progress$set(message = "Reducing dimensions", value = 5)
      scseq <- run_umap(scseq)
      gc()

      progress$set(message = "Getting markers", value = 6)
      markers <- get_scseq_markers(scseq)
      gc()

      progress$set(message = "Saving", value = 7)
      anal <- list(scseq = scseq, markers = markers, annot = names(markers))
      save_scseq_data(anal, dataset_name, sc_dir)

      # get and save cluster stats for selectizeInputs
      get_cluster_stats(sc_dir, dataset_name, scseq)

    } else {
      # Create a callback function to update progress.
      progress$set(message = "Quantifying files", value = 0)
      updateProgress <- function(amount = NULL, detail = NULL) {
        progress$inc(amount = amount, detail = detail)
      }

      # setup bulk
      pdata <- dsQuantTable$pdata()
      paired <- dsForm$paired()


      # quantification
      run_kallisto_bulk(indices_dir = indices_dir,
                        data_dir = fastq_dir,
                        pdata = pdata,
                        paired = paired,
                        updateProgress = updateProgress)

      # generate eset and save
      progress$set(message = 'Annotating dataset')
      eset <- load_seq(fastq_dir)

      # save to bulk datasets to indicate that has been quantified
      save_bulk_dataset(dataset_name, dataset_dir, data_dir)

    }
    # trigger to update rest of app
    new_dataset(dataset_name)

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
    new_dataset = new_dataset,
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
    file.path(data_dir, dataset_dir, group_file)
  })

  anal <- reactive({
    anal_path <- anal_path()
    req(file.exists(anal_path))
    readRDS(anal_path)
  })

  # MDS plot
  output$plotly <- plotly::renderPlotly({
    mds <- anal()$mds
    plotlyMDS(scaling = mds$scaling, scaling_sva = mds$scaling_sva)
  })


}

#' Logic for Dataset Gene plotly
#' @export
#' @keywords internal
dsGenePlotly <- function(input, output, session, eset, explore_genes, pdata, dataset_name) {

  # MDS plot
  output$plotly <- plotly::renderPlotly({
    # need at least two groups
    pdata <- pdata()
    pdata <- pdata[!is.na(pdata$Group), ]
    req(length(unique(pdata$Group)) > 1)

    # need eset and at least one gene
    eset <- eset()
    explore_genes <- explore_genes()
    req(eset, explore_genes)

    plotlyGene(eset, explore_genes, pdata, dataset_name())
  })


}


#' Logic for Datasets form
#' @export
#' @keywords internal
dsForm <- function(input, output, session, data_dir, sc_dir, new_dataset, msg_quant, msg_anal, new_anal) {

  dataset <- callModule(dsDataset, 'selected_dataset',
                        data_dir = data_dir,
                        new_dataset = new_dataset)


  eset <- reactive({
    fastq_dir <- dataset$fastq_dir()
    eset_path <- file.path(fastq_dir, 'eset.rds')

    req(file.exists(eset_path))
    readRDS(eset_path)
  })


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
                      error_msg = msg_quant,
                      is.sc = dataset$is.sc)

  anal <- callModule(dsFormAnal, 'anal_form',
                     data_dir = data_dir,
                     sc_dir = sc_dir,
                     error_msg = msg_anal,
                     dataset_name = dataset$dataset_name,
                     dataset_dir = dataset$dataset_dir,
                     new_anal = new_anal,
                     new_dataset = new_dataset,
                     eset = eset)



  return(list(
    fastq_dir = dataset$fastq_dir,
    paired = quant$paired,
    quant_labels = quant$labels,
    anal_labels = anal$labels,
    explore_labels = anal$explore_labels,
    explore_genes = anal$explore_genes,
    run_quant = quant$run_quant,
    run_anal = anal$run_anal,
    dataset_name = dataset$dataset_name,
    dataset_dir = dataset$dataset_dir,
    anal_name = anal$anal_name,
    show_quant = show_quant,
    show_anal = show_anal,
    is.sc = dataset$is.sc,
    is.cellranger = dataset$is.cellranger,
    is.explore = anal$is.explore,
    eset = eset
  ))

}

#' Logic for selected dataset part of dsFrom
#' @export
#' @keywords internal
dsDataset <- function(input, output, session, data_dir, new_dataset) {


  # get directory with fastqs
  roots <- c('data_dir' = data_dir)
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
    updateSelectizeInput(session, 'dataset_name', choices = rbind(rep(NA, 5), datasets()), server = TRUE)
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
      req(!'integer' %in% class(input$dataset_dir))
      type <- input$dataset_dir$path[[2]]
      dir <- shinyFiles::parseDirPath(roots, input$dataset_dir)
      dir <- as.character(dir)
      dir <- gsub(paste0(data_dir, '/'), '', dir)
    }
    else {
      datasets <- datasets()
      req(datasets)
      dir <- datasets[datasets$dataset_name == dataset_name, 'dataset_dir']
    }
    return(dir)
  })


  fastq_dir <- reactive(file.path(data_dir, dataset_dir()))

  is.sc <- reactive(grepl('^single-cell/', dataset_dir()))
  is.cellranger <- reactive(check_is_cellranger(fastq_dir()))

  observe(if (is.cellranger()) standardize_cellranger(fastq_dir()))

  return(list(
    fastq_dir = fastq_dir,
    dataset_name = dataset_name,
    dataset_dir = dataset_dir,
    is.sc = is.sc,
    is.cellranger = is.cellranger,
    is.create = is.create
  ))

}


#' Logic for dataset quantification part of dsForm
#' @export
#' @keywords internal
dsFormQuant <- function(input, output, session, fastq_dir, error_msg, is.sc) {


  paired <- callModule(dsEndType, 'end_type',
                       fastq_dir = fastq_dir,
                       is.sc = is.sc)

  observe(shinyjs::toggleClass("pair", 'disabled', condition = !paired()))

  observe({
    toggle('bulk_controls', condition = !is.sc())
  })


  reset <- reactive(input$reset)
  rep <- reactive(input$rep)
  pair <- reactive(input$pair)


  observe({
    error_msg <- error_msg()
    toggleClass('quant_labels', 'has-error', condition = !is.null(error_msg))
    html('error_msg', html = error_msg)

  })

  quantModal <- function(is.sc, paired) {
    end_type <- ifelse(paired, 'pair ended', 'single ended')

    if (!is.sc) {
      UI <- withTags({
        dl(
          dt('End type:'),
          dd(paste('Experiment is selected as', end_type)),
          hr(),
          dt('Replicates:'),
          dd('If any - e.g. same sample sequenced in replicate. These will be treated as a single library.')
        )
      })

    } else {
      UI <- withTags({
        dl(
          dt('Are you sure?'),
          dd('This will take a while.')
        )
      })
    }

    modalDialog(
      UI,
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
    showModal(quantModal(is.sc(), paired()))
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
dsEndType <- function(input, output, session, fastq_dir, is.sc) {


  # get fastq files in directory
  fastq_files <- reactive({
    fastq_dir <- fastq_dir()
    req(fastq_dir)

    list.files(fastq_dir, '.fastq.gz$')
  })

  # auto detected if paired
  detected_paired <- reactive({
    if (is.sc()) return(TRUE)
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
dsFormAnal <- function(input, output, session, data_dir, sc_dir, error_msg, dataset_name, dataset_dir, new_anal, new_dataset, eset) {


  run_anal <- reactiveVal()
  labels <- list(
    test = reactive(input$test),
    ctrl = reactive(input$ctrl),
    reset = reactive(input$reset)
  )

  # toggle analysis type
  is.explore <- reactive({
    input$anal_type == 'exploratory'
  })

  observe({
    toggle('diff_panel', condition = !is.explore())
    toggle('explore_panel', condition = is.explore())
  })

  # toggle cell-type deconvolution
  show_decon <- reactive(input$show_decon %% 2 != 0)

  observe({
    toggleClass(id = "show_decon", 'btn-primary', condition = show_decon())
  })

  deconForm <- callModule(deconvolutionForm, 'decon',
                          show_decon = show_decon,
                          new_dataset = new_dataset,
                          sc_dir = sc_dir)


  # logic for group name buttons
  explore_group_name <- reactive(input$explore_group_name)

  observeEvent(input$grouped, {
    updateTextInput(session, 'explore_group_name', value = '')
  })

  explore_labels <- list(
    grouped = reactive(input$grouped),
    reset = reactive(input$reset_explore),
    explore_group_name = explore_group_name
  )

  anal_name <- reactive(input$anal_name)
  has_anal_name <- reactive(anal_name() != '')

  # gene choices
  observe({
    updateSelectizeInput(session, 'explore_genes', choices = c(NA, row.names(eset())), server = TRUE)
  })


  # -------------------------------
  # Differential Expression Section
  # -------------------------------
  # TODO: refactor into seperate module

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

  # show anal buttons if new analysis name
  observe({
    toggle('anal_buttons_panel', condition = !is_prev_anal() & isTruthy(anal_name()))
  })

  # enable download results if previous analysis
  observe({
    toggleState('download', condition = is_prev_anal())
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
    dataset_dir <- dataset_dir()

    req(anal_name, dataset_dir)

    # path to analysis result
    anal_file <- paste0('diff_expr_symbol_', anal_name, '.rds')
    anal_path <- file.path(data_dir, dataset_dir, anal_file)
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
    explore_labels = explore_labels,
    explore_genes = reactive(input$explore_genes),
    run_anal = run_anal,
    anal_name = anal_name,
    is.explore = is.explore
  ))
}

#' Logic for deconvolution form
#' @export
#' @keywords internal
deconvolutionForm <- function(input, output, session, show_decon, new_dataset, sc_dir) {
  exclude_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))

  # show deconvolution form toggle
  observe({
    toggle(id = "decon_form", anim = TRUE, condition = show_decon())
  })

  # available single cell datasets for deconvolution
  ref_anals <- reactive({

    # reactive to new sc datasets
    new_dataset()

    # make sure integrated rds exists
    int_path <- file.path(sc_dir, 'integrated.rds')
    if (!file.exists(int_path)) saveRDS(NULL, int_path)

    # use saved anals as options
    integrated <- readRDS(file.path(sc_dir, 'integrated.rds'))
    individual <- setdiff(list.files(sc_dir), c(integrated, 'integrated.rds'))
    return(individual)
  })



  # update reference dataset choices
  observe({
    ref_anals <- ref_anals()
    req(ref_anals)
    updateSelectizeInput(session, 'decon_anal', choices = c('', ref_anals))
  })

  annot <- reactive({
    anal_name <- input$decon_anal
    req(anal_name)
    annot_path <- scseq_part_path(sc_dir, anal_name, 'annot')
    readRDS(annot_path)
  })

  # update exclude cluster choices
  exclude_choices <- reactive({
    clusters <- annot()
    anal_name <- input$decon_anal
    get_cluster_choices(clusters, anal_name, sc_dir)
  })

  observe({
    choices <- exclude_choices()
    updateSelectizeInput(session, 'exclude_clusters', choices = choices, options = exclude_options, server = TRUE)
  })

  observeEvent(input$submit_decon, {
    browser()

  })
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
      if (!length(fastqs)) fastqs <- NA
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
dsAnalTable <- function(input, output, session, eset, labels, data_dir, dataset_dir, anal_name) {

  background <- '#337ab7 url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAAPklEQVQoU43Myw0AIAgEUbdAq7VADCQaPyww55dBKyQiHZkzBIwQLqQzCk9E4Ytc6KEPMnTBCG2YIYMVpHAC84EnVbOkv3wAAAAASUVORK5CYII=) repeat'

  # things user will update and return
  pdata_r <- reactiveVal()
  group_r <- reactiveVal()


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
    file.path(data_dir, dataset_dir, group_file)
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
      options = list(
        columnDefs = list(list(className = 'dt-nopad', targets = 0)),
        scrollX = TRUE,
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


#' Logic for differential expression analysis table
#' @export
#' @keywords internal
dsExploreTable <- function(input, output, session, eset, labels, data_dir, dataset_dir) {

  # colors
  group_colors <- RColorBrewer::brewer.pal(8, 'Set2')
  ncolors <- length(group_colors)

  # things user will update and return
  pdata_r <- reactiveVal()
  group_r <- reactiveVal()
  name_r <- reactiveVal()
  table_rendered <- reactiveVal()


  # path to saved pdata
  pdata_path <- reactive(file.path(data_dir, dataset_dir(), 'pdata_explore.rds'))


  # pdata that gets returned with group column
  returned_pdata <- reactive({
    pdata_path <- pdata_path()
    pdata <- pdata_r()
    req(pdata, pdata_path)

    pdata$Group <- group_r()
    pdata$`Group name` <- name_r()

    saveRDS(pdata, pdata_path)
    return(pdata)
  })

  # reset when new eset or analysis name
  observe({
    eset <- eset()
    pdata_path <- pdata_path()
    req(eset, pdata_path)

    pdata <- Biobase::pData(eset) %>%
      dplyr::select(-lib.size, -norm.factors) %>%
      tibble::add_column(Group = NA, 'Group name' = NA, Title = colnames(eset), .before = 1)

    group <- name <- rep(NA, nrow(pdata))
    # load pdata from previous if available
    if (file.exists(pdata_path)) {
      saved <- readRDS(pdata_path)
      group <- saved$Group
      name  <- saved$`Group name`
    }

    table_rendered(FALSE)

    pdata_r(pdata)
    group_r(group)
    name_r(name)
  })


  # redraw table when new pdata (otherwise update data using proxy)
  output$pdata <- DT::renderDataTable({

    table_rendered(TRUE)

    DT::datatable(
      pdata_r(),
      class = 'cell-border dt-fake-height',
      rownames = FALSE,
      escape = FALSE, # to allow HTML in table
      options = list(
        columnDefs = list(list(className = 'dt-nopad', targets = 0)),
        scrollX = TRUE,
        paging = FALSE,
        bInfo = 0
      )
    )
  })


  html_pdata <- reactive({

    # things that trigger update
    pdata <- returned_pdata()
    group <- pdata$Group
    name  <- pdata$`Group name`
    req(pdata)


    # update pdata Group column
    not.na <- !is.na(group)
    group_nums  <- unique(group[not.na])
    group_names <- unique(name[not.na])
    ind <- order(group_nums)

    group_nums  <- group_nums[ind]
    group_names <- group_names[ind]

    for (i in seq_along(group_nums)) {
      group_num <- group_nums[i]
      group_name <- group_names[i]

      # plotly color bug when two groups
      color <- group_colors[group_num]

      rows <- which(group == group_num)
      pdata[rows, 'Group'] <- paste('<div style="background-color:', color, ';"></div>')
      pdata[rows, 'Group name'] <- group_name
    }

    return(pdata)
  })


  # proxy used to replace data
  proxy <- DT::dataTableProxy("pdata")
  shiny::observe({
    req(html_pdata(), table_rendered())
    DT::replaceData(proxy, html_pdata(), rownames = FALSE)
  })



  # click 'grouped'
  shiny::observeEvent(labels$grouped(), {
    group_name <- labels$explore_group_name()
    req(group_name)
    group <- group_r()
    name <- name_r()

    rows  <- input$pdata_rows_selected
    group_num <- length(unique(setdiff(group, NA))) + 1


    group[rows] <- group_num
    name[rows] <- group_name
    group_r(group)
    name_r(name)
  })



  # click 'Reset'
  shiny::observeEvent(labels$reset(), {
    req(labels$reset())
    group <- group_r()
    clear <- rep(NA, length(group))
    group_r(clear)
    name_r(clear)

  })


  return(list(
    pdata = returned_pdata
  ))

}
