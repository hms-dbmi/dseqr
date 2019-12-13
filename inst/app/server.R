#' Logic Datasets page
#' @export
#' @keywords internal
dsPage <- function(input, output, session, data_dir, sc_dir, bulk_dir, indices_dir) {

  new_anal <- reactiveVal()
  new_dataset <- reactiveVal()
  msg_quant <- reactiveVal()

  # paths to objects used for normalization
  eset_paths <- reactive({
    fastq_dir <- dsForm$fastq_dir()
    types <- c('eset', 'pdata_explore', 'vsd', 'svobj', 'numsv', 'adj')
    eset_paths <- file.path(fastq_dir, paste0(types, '.rds'))
    names(eset_paths) <- types
    return(eset_paths)
  })


  # explore_eset used for all plots
  explore_eset <- reactive({
    # update when change dataset or manually trigger
    dirs <- eset_paths()
    input$run_norm_and_sva

    # minimum need saved eset and pdata_explore
    need <- dirs[c('eset', 'pdata_explore')]
    req(all(file.exists(need)))

    # also need that pdata_explore has more than two groups
    pdata <- readRDS(dirs['pdata_explore'])
    keep <- row.names(pdata)[!is.na(pdata$Group)]
    pdata <- pdata[keep, ]
    req(length(unique(pdata$Group)) > 1)

    # subset eset and add explore_pdata
    eset <- readRDS(dirs['eset'])
    eset <- eset[, keep]
    pdata$group <- pdata$`Group name`
    Biobase::pData(eset) <- pdata

    # rlog normalize
    eset <- add_vsd(eset, vsd_path = dirs['vsd'])

    # adjust for pairs/surrogate variables
    svobj_dir <- dirs['svobj']
    numsv_dir <- dirs['numsv']
    svobj <- ifelse(file.exists(svobj_dir), readRDS(svobj_dir), list(sv = NULL))
    numsv <- ifelse(file.exists(numsv_dir), readRDS(numsv_dir), 0)

    # adjusted path is dynamic
    adj_file <- ifelse(numsv > 0, paste0('adjusted_', numsv, 'svs.rds'), 'adjusted.rds')
    adj_path <- gsub('adj.rds$', adj_file, dirs['adj'])

    eset <- add_adjusted(eset, svobj, num_svs, adj_path = adj_path)

    # use SYMBOL as annotation
    # keep unique symbol based on row IQRs
    eset <- iqr_replicates(eset)

    return(eset)

  })


  dsForm <- callModule(dsForm, 'form',
                       data_dir = data_dir,
                       sc_dir = sc_dir,
                       bulk_dir = bulk_dir,
                       new_dataset = new_dataset,
                       msg_quant = msg_quant,
                       new_anal = new_anal,
                       explore_eset = explore_eset)

  callModule(dsMDSplotly, 'mds_plotly', explore_eset = explore_eset)


  # toggle tables
  observe({
    toggle('quant_table_container', condition = dsForm$show_quant())
    toggle('anal_table_container', condition = dsForm$show_anal())
  })

  # toggle plots
  sel_genes <- reactive(length(dsForm$explore_genes() > 0))

  observe({
    toggle('mds_plotly_container', condition = !sel_genes() & !dsForm$show_dtangle())
    toggle('gene_plotly_container', condition = !dsForm$show_dtangle() & sel_genes())
    toggle('cells_plotly_container', condition = dsForm$show_dtangle())
  })

  dsQuantTable <- callModule(dsQuantTable, 'quant',
                             fastq_dir = dsForm$fastq_dir,
                             labels = dsForm$quant_labels,
                             paired = dsForm$paired)


  dsExploreTable <- callModule(dsExploreTable, 'explore',
                               eset_paths = eset_paths,
                               labels = dsForm$explore_labels,
                               data_dir = data_dir,
                               dataset_dir = dsForm$dataset_dir)


  callModule(dsGenePlotly, 'gene_plotly',
             eset = explore_eset,
             explore_genes = dsForm$explore_genes,
             dataset_name = dsForm$dataset_name)

  callModule(dsCellsPlotly, 'cells_plotly',
             dtangle_est = dsForm$dtangle_est,
             pdata = dsExploreTable$pdata,
             dataset_name = dsForm$dataset_name)

  observe({
    msg_quant(dsQuantTable$valid_msg())
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
    eset <- explore_eset()
    fastq_dir <- dsForm$fastq_dir()

    dataset_name <- dsForm$dataset_name()
    dataset_dir <- dsForm$dataset_dir()
    anal_name <- dsForm$anal_name()
    contrast <- gsub('_vs_', '-', anal_name)

    svobj <- dsForm$svobj()
    num_svs <- dsForm$num_svs()

    req(eset, fastq_dir, anal_name, dataset_dir)

    # run differential expression
    progress$set(message = "Differential expression", value = 1)
    browser()
    anal <- diff_expr(eset,
                      data_dir = fastq_dir,
                      anal_name = anal_name,
                      contrast = contrast,
                      svobj = svobj,
                      num_svs = num_svs,
                      prev_anal = list(pdata = Biobase::pData(eset)))

    # run pathway analysis
    # progress$set(message = "Pathway analysis", value = 2)
    # rna_seq <- 'lib.size' %in% colnames(Biobase::pData(eset))
    # path_anal <- diff_path(eset, prev_anal = anal, data_dir = fastq_dir, anal_name = anal_name, rna_seq = rna_seq, NI = 24)

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
dsMDSplotly <- function(input, output, session, explore_eset) {

  # MDS plot
  plotly_fun <- reactive({
    eset <- explore_eset()
    pdata <- Biobase::pData(eset)

    # setup group factor and colors
    group <- pdata$`Group name`
    group_order <- order(unique(pdata$Group))
    group_levels <- unique(group)[group_order]
    group <- factor(group, levels = group_levels)
    group_colors <- RColorBrewer::brewer.pal(8, 'Set2')[seq_along(group_levels)]

    vsd <- Biobase::assayDataElement(eset, 'vsd')
    adj <- Biobase::assayDataElement(eset, 'adjusted')
    mds <- get_mds(vsd, adj, group)
    plotlyMDS(mds$scaling, mds$scaling_adj, group_colors = group_colors)
  })


  # not currently downloadable so use default fname_fun and data_fun
  callModule(downloadablePlotly, 'plotly', plotly_fun = plotly_fun)
}


#' Logic for Dataset Gene plotly
#' @export
#' @keywords internal
dsGenePlotly <- function(input, output, session, eset, explore_genes, dataset_name) {

  boxplotly_args <- reactive({
    # need eset and at least one gene
    explore_genes <- explore_genes()
    req(explore_genes)
    eset <- eset()

    # prevent warning when switching between dataset (eset updated before genes)
    req(all(explore_genes %in% row.names(eset)))

    get_boxplotly_gene_args(eset, explore_genes, dataset_name())
  })

  plotly_fun <- reactive({

    args <- boxplotly_args()

    boxPlotly(df = args$df,
              boxgap = args$boxgap,
              boxgroupgap = args$boxgroupgap,
              plot_fname = args$plot_fname,
              ytitle = 'Normalized Expression',
              xtitle = 'Gene')
  })

  fname_fun <- function() {
    paste(dataset_name(), '_', paste(explore_genes(), collapse = '_'), '_', Sys.Date(), ".csv", sep = "")
  }


  data_fun <- function(file) {
    args <- boxplotly_args()
    df <- args$df
    df$color <- NULL
    colnames(df) <- c('Sample', 'Gene', 'Normalized Expression', 'Group')
    write.csv(df, file, row.names = FALSE)
  }

  callModule(downloadablePlotly, 'plotly', plotly_fun = plotly_fun, fname_fun = fname_fun, data_fun = data_fun)



}


#' Logic for Dataset Cell Type deconvolution plotly
#' @export
#' @keywords internal
dsCellsPlotly <- function(input, output, session, dtangle_est, pdata, dataset_name) {

  boxplotly_args <- reactive({
    # need at least two groups
    pdata <- pdata()
    pdata <- pdata[!is.na(pdata$Group), ]
    req(length(unique(pdata$Group)) > 1)

    dtangle_est <- dtangle_est()
    req(dtangle_est)

    get_boxplotly_cell_args(pdata, dtangle_est, dataset_name())
  })

  plotly_fun <- reactive({

    args <- boxplotly_args()

    boxPlotly(df = args$df,
              boxgap = args$boxgap,
              boxgroupgap = args$boxgroupgap,
              plot_fname = args$plot_fname,
              ytitle = 'Estimated Proportion',
              xtitle = 'Cluster')
  })

  fname_fun <- function() {
    clusters <- colnames(dtangle_est())
    clusters <- gsub(' ', '', clusters)
    paste(dataset_name(), '_', paste(clusters, collapse = '_'), '_', Sys.Date(), ".csv", sep = "")
  }


  data_fun <- function(file) {
    args <- boxplotly_args()
    df <- args$df
    df$color <- NULL
    colnames(df) <- c('Sample', 'Group', 'Proportion', 'Cluster')
    write.csv(df, file, row.names = FALSE)
  }

  callModule(downloadablePlotly, 'plotly', plotly_fun = plotly_fun, fname_fun = fname_fun, data_fun = data_fun)


}


#' Logic for Datasets form
#' @export
#' @keywords internal
dsForm <- function(input, output, session, data_dir, sc_dir, bulk_dir, new_dataset, msg_quant, new_anal, explore_eset) {

  dataset <- callModule(dsDataset, 'selected_dataset',
                        data_dir = data_dir,
                        sc_dir = sc_dir,
                        bulk_dir = bulk_dir,
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
                      error_msg = msg_quant,
                      is.sc = dataset$is.sc)

  anal <- callModule(dsFormAnal, 'anal_form',
                     data_dir = data_dir,
                     sc_dir = sc_dir,
                     bulk_dir = bulk_dir,
                     dataset_name = dataset$dataset_name,
                     dataset_dir = dataset$dataset_dir,
                     new_anal = new_anal,
                     new_dataset = new_dataset,
                     explore_eset = explore_eset)



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
    dtangle_est = dataset$dtangle_est,
    show_dtangle = dataset$show_dtangle,
    selected_nsv = dataset$selected_nsv
  ))

}

#' Logic for selected dataset part of dsFrom
#' @export
#' @keywords internal
dsDataset <- function(input, output, session, sc_dir, bulk_dir, data_dir, new_dataset) {


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

  is.existing <- reactive({
    dataset_name <- dataset_name()
    datasets <- datasets()
    dataset_name %in% datasets$dataset_name
  })

  observe({
    req(datasets())
    updateSelectizeInput(session, 'dataset_name', choices = rbind(rep(NA, 5), datasets()), server = TRUE)
  })

  # open selector if creating
  observe({
    req(is.create())
    if (is.create()) {
      shinyjs::click('dataset_dir')
    }
  })

  # enable dataset button for existing datasets
  observe({
    toggleState('run_norm_and_sva', condition = is.existing())
    toggleState('show_nsv', condition = is.existing())
    toggleState('show_dtangle', condition = is.existing())
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


  # toggle cell-type deconvolution
  show_dtangle <- reactive(input$show_dtangle %% 2 != 0)

  observe({
    toggleClass(id = "show_dtangle", 'btn-primary', condition = show_dtangle())
  })

  dtangleForm <- callModule(dtangleForm, 'dtangle',
                            show_dtangle = show_dtangle,
                            new_dataset = new_dataset,
                            sc_dir = sc_dir,
                            bulk_dir = bulk_dir,
                            dataset_name = dataset_name,
                            explore_eset = explore_eset,
                            dataset_dir = dataset_dir)

  return(list(
    fastq_dir = fastq_dir,
    dataset_name = dataset_name,
    dataset_dir = dataset_dir,
    is.sc = is.sc,
    is.cellranger = is.cellranger,
    is.create = is.create,
    dtangle_est = dtangleForm$dtangle_est,
    show_dtangle = show_dtangle,
    selected_nsv = reactive(input$selected_nsv)
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
dsFormAnal <- function(input, output, session, data_dir, sc_dir, bulk_dir, dataset_name, dataset_dir, new_anal, new_dataset, explore_eset) {
  contrast_options <- list(render = I('{option: bulkContrastOptions, item: bulkContrastItem}'))

  svobj <- reactiveVal()
  numsv <- reactiveVal()
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


  # remove dataset files that depend on groupings
  observe({
    input$grouped
    input$reset_explore

    data_dir <- file.path(data_dir, isolate(dataset_dir()))
    remove_dataset_files(data_dir)
    numsv(0)
    svobj(NULL)
  })


  # gene choices
  observe({
    updateSelectizeInput(session, 'explore_genes', choices = c(NA, row.names(explore_eset())), server = TRUE)
  })

  # path to saved svobj
  svobj_path <- reactive(file.path(data_dir, dataset_dir(), 'svobj.rds'))
  numsv_path <- reactive(file.path(data_dir, dataset_dir(), 'numsv.rds'))

  # obdate svobj and numsv when paths update
  observe({
    svobj_path <- svobj_path()
    req(file.exists(svobj_path))
    svobj(readRDS(svobj_path))
  })

  observe({
    numsv_path <- numsv_path()
    if (file.exists(numsv_path)) {
      numsv(readRDS(numsv_path))
    } else {
      numsv(0)
    }
  })

  # run surrogate variable analysis on click
  observeEvent(input$run_sva, {


    # remove previously adjusted data
    data_dir <- file.path(data_dir, dataset_dir())
    remove_dataset_files(data_dir, c('^adjusted_\\d+svs.rds$', '^numsv.rds$'))
    numsv(0)

    eset <- explore_eset()
    pdata <- Biobase::pData(eset)
    group <- pdata$group
    req(length(unique(group)) > 1)

    mods <- get_mods(eset)
    rna_seq <- 'lib.size' %in% colnames(pdata)
    svobj_new <- run_sva(mods, eset, rna_seq = rna_seq)

    saveRDS(svobj_new, svobj_path())
    svobj(svobj_new)
  })

  # update number of surrogate variables slider
  observe({
    svobj <- svobj()
    updateSliderInput(session, 'num_svs', value = numsv(), min = 0, max = svobj$n.sv)
  })

  observeEvent(input$num_svs, {
    numsv(input$num_svs)
    saveRDS(input$num_svs, numsv_path())
  }, ignoreInit = TRUE)

  # update button icon with selected number of surrogate variables
  observe({
    updateActionButton(session, 'show_nsv', label = htmltools::doRenderTags(tags$span(input$num_svs, class='fa fa-fw')))
  })


  # -------------------------------
  # Differential Expression Section
  # -------------------------------
  # TODO: refactor into seperate module

  anal_name <- reactive({
    contrast_groups <- input$contrast_groups
    req(length(contrast_groups == 2))

    return(paste0(contrast_groups[1], '_vs_', contrast_groups[2]))

  })
  has_anal_name <- reactive(isTruthy(anal_name()))

  # analysis info table (can be multiple) from dataset
  dataset_anals <- reactive({
    # reload if new analysis
    new_anal()
    dataset_name <- dataset_name()
    req(dataset_name)

    anals <- load_bulk_anals(data_dir)
    c('', anals[anals$dataset_name == dataset_name, 'anal_name'])
  })


  # group levels used for selecting test and control groups
  group_levels <- reactive({
    eset <- explore_eset()
    pdata <- Biobase::pData(eset)
    group <- pdata$`Group name`
    group_order <- order(unique(pdata$Group))
    group_levels <- unique(group)[group_order]

    group_colors <- RColorBrewer::brewer.pal(8, 'Set2')

    data.frame(
      name = group_levels,
      value = group_levels,
      color = group_colors[seq_along(group_levels)], stringsAsFactors = FALSE
    )
  })

  observe({
    updateSelectizeInput(session, 'contrast_groups', choices = group_levels(), server = TRUE, options = contrast_options)
  })

  is_prev_anal <- reactive({
    anal_name() %in% setdiff(dataset_anals(), '')
  })


  # enable download results if previous analysis
  observe({
    toggleState('download', condition = is_prev_anal())
  })

  # enable running analysis
  full_contrast <- reactive(length(input$contrast_groups) == 2)
  observe({
    toggleState('run_sva', condition = full_contrast())
  })




  observeEvent(input$run_anal, {
    # check for analysis name
    if (!has_anal_name()) {
      addClass('run_anal_container', 'has-error')
      html('run_anal_help', 'Select two groups.')
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
    is.explore = is.explore,
    num_svs = reactive(input$num_svs),
    svobj = svobj
  ))
}

#' Logic for deconvolution form
#' @export
#' @keywords internal
dtangleForm <- function(input, output, session, show_dtangle, new_dataset, sc_dir, bulk_dir, explore_eset, dataset_dir, dataset_name) {
  include_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))
  input_ids <- c('include_clusters', 'dtangle_anal', 'submit_dtangle')

  dtangle_est <- reactiveVal()

  # show deconvolution form toggle
  observe({
    toggle(id = "dtangle_form", anim = TRUE, condition = show_dtangle())
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
    updateSelectizeInput(session, 'dtangle_anal', choices = c('', ref_anals))
  })

  annot <- reactive({
    anal_name <- input$dtangle_anal
    req(anal_name)
    annot_path <- scseq_part_path(sc_dir, anal_name, 'annot')
    readRDS(annot_path)
  })

  # update exclude cluster choices
  include_choices <- reactive({
    clusters <- annot()
    anal_name <- input$dtangle_anal
    get_cluster_choices(clusters, anal_name, sc_dir)
  })

  observe({
    choices <- include_choices()
    updateSelectizeInput(session, 'include_clusters', choices = choices, options = include_options, server = TRUE)
  })

  # scseq for deconvolution
  scseq <- reactive({
    anal_name <- input$dtangle_anal
    scseq_path <- scseq_part_path(sc_dir, anal_name, 'scseq')
    readRDS(scseq_path)
  })

  observeEvent(input$submit_dtangle, {

    # disable inputs
    toggleAll(input_ids)

    # Create a Progress object
    progress <- Progress$new(session, min=0, max = 4)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())

    progress$set(message = "Deconvoluting", value = 0)
    progress$set(value = 1)

    anal_name <- input$dtangle_anal
    dataset_name <- dataset_name()

    # get names of clusters
    include_clusters <- input$include_clusters
    include_choices <- include_choices()
    row.names(include_choices) <- include_choices$value
    include_names <- include_choices[include_clusters, ]$name

    # require at least two labeled groups
    eset <- explore_eset()
    pdata <- Biobase::pData(eset)

    # if select none deconvolute using all clusters
    if (!length(include_clusters))
      include_clusters <- as.character(include_choices$value)

    # subset to selected clusters
    scseq <- scseq()
    scseq <- scseq[, scseq$seurat_clusters %in% include_clusters]

    # get DESeq2::vst normalized values from eset
    vsd <- Biobase::assayDataElement(eset, 'vsd')

    # common genes only
    commongenes <- intersect (rownames(vsd), rownames(scseq))
    vsd <- vsd[commongenes, ]
    scseq <- scseq[commongenes, ]

    # quantile normalize scseq and rnaseq dataset
    progress$set(value = 2)
    y <- cbind(as.matrix(scseq[['SCT']]@data), vsd)

    y <- limma::normalizeBetweenArrays(y)
    y <- t(y)


    progress$set(value = 3)
    # indicies for cells in each included cluster
    pure_samples <- list()
    for (i in seq_along(include_clusters))
      pure_samples[[include_names[i]]] <-
      which(scseq$seurat_clusters == include_clusters[i])

    # markers for each included cluster
    marker_list = dtangle::find_markers(y,
                                        pure_samples = pure_samples,
                                        data_type = "rna-seq",
                                        marker_method='ratio')

    # use markers in top 10th quantile with a minimum of 3
    q = 0.1
    quantiles = lapply(marker_list$V,function(x) quantile(x,1-q))
    K = length(pure_samples)
    n_markers = sapply(seq_len(K),function(i){
      max(3, which(marker_list$V[[i]] > quantiles[[i]]))
    })

    # run deconvolution and get get proportion estimates
    marks <- marker_list$L
    dc <- dtangle::dtangle(y,
                           pure_samples = pure_samples,
                           n_markers = n_markers,
                           data_type = 'rna-seq',
                           markers = marks)

    dc <- dc$estimates[colnames(eset), ]
    dtangle_est(dc)
    toggleAll(input_ids)
    progress$set(value = 4)
  })


  return(list(
    dtangle_est = dtangle_est
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
dsExploreTable <- function(input, output, session, eset_paths, labels, data_dir, dataset_dir) {

  # colors
  group_colors <- RColorBrewer::brewer.pal(8, 'Set2')
  ncolors <- length(group_colors)

  # things user will update and return
  pdata_r <- reactiveVal()
  group_r <- reactiveVal()
  name_r <- reactiveVal()
  table_rendered <- reactiveVal()


  # pdata that gets returned with group column
  returned_pdata <- reactive({

    dirs <- eset_paths()
    pdata_path <- dirs['pdata_explore']
    pdata <- pdata_r()
    req(pdata)

    pdata$Group <- group_r()
    pdata$`Group name` <- name_r()

    saveRDS(pdata, pdata_path)
    return(pdata)
  })

  # initialize saved pdata_explore
  observe({
    dirs <- eset_paths()
    pdata_path <- dirs['pdata_explore']

    if (!file.exists(pdata_path)) {
      eset <- readRDS(dirs['eset'])

      pdata <- Biobase::pData(eset) %>%
        tibble::add_column(Group = NA, 'Group name' = NA, Title = colnames(eset), .before = 1)

      saveRDS(pdata, pdata_path)
    }
  })


  # reset when new dataset selected
  observe({
    dirs <- eset_paths()
    pdata_path <- dirs['pdata_explore']
    req(file.exists(pdata_path))

    pdata <- readRDS(pdata_path)
    table_rendered(FALSE)

    pdata_r(pdata)
    group_r(pdata$Group)
    name_r(pdata$`Group name`)
  })


  # redraw table when new pdata (otherwise update data using proxy)
  output$pdata <- DT::renderDataTable({

    pdata <- pdata_r()
    hide_target <- which(colnames(pdata) %in% c('lib.size', 'norm.factors', 'pair')) - 1

    table_rendered(TRUE)

    DT::datatable(
      pdata,
      class = 'cell-border dt-fake-height',
      rownames = FALSE,
      escape = FALSE, # to allow HTML in table
      options = list(
        columnDefs = list(
          list(className = 'dt-nopad', targets = 0),
          list(targets = hide_target, visible=FALSE)),
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


server <- function(input, output, session) {

  # get arguments from calling function
  # defaults for testing

  # base directory contains data_dir folder
  data_dir <- getShinyOption('data_dir', 'tests/data/test/example')

  # path where pert queries will be stored
  pert_query_dir <- getShinyOption('pert_query_dir', '/srv/drugseqr/pert_query_dir')

  # path where pert signatures will be stored
  pert_signature_dir <- getShinyOption('pert_signature_dir', '/srv/drugseqr/pert_signature_dir')

  # path where kallisto index is downloaded and stored
  indices_dir <- getShinyOption('indices_dir', '/srv/drugseqr/indices')

  if (!dir.exists(pert_query_dir)) dir.create(pert_query_dir)
  if (!dir.exists(pert_signature_dir)) dir.create(pert_signature_dir)

  # for testing don't seem to be able to pass arguments as options
  if (isTRUE(getOption('shiny.testmode'))) {
    # reset data for testing
    data_dir <- 'tests/data/test/example'
    static_dir <- 'tests/data/static/example'
    unlink(data_dir, recursive = TRUE)
    dir.create(data_dir)
    file.copy(list.files(static_dir, full.names = TRUE), data_dir, recursive = TRUE)
  }

  sc_dir <- file.path(data_dir, 'single-cell')
  bulk_dir <- file.path(data_dir, 'bulk')

  dir.create(sc_dir, showWarnings = FALSE)
  dir.create(bulk_dir, showWarnings = FALSE)



  dsPage <- callModule(dsPage, 'bulk',
                       data_dir = data_dir,
                       sc_dir = sc_dir,
                       bulk_dir = bulk_dir,
                       indices_dir = indices_dir)

  scPage <- callModule(scPage, 'sc',
                       sc_dir = sc_dir,
                       new_dataset = dsPage$new_dataset)

  drugsPage <- callModule(drugsPage, 'drug',
                          new_anal = dsPage$new_anal,
                          data_dir = data_dir,
                          pert_query_dir = pert_query_dir,
                          pert_signature_dir = pert_signature_dir)

  pathPage <- callModule(pathPage, 'pathways',
                         new_anal = dsPage$new_anal,
                         data_dir = data_dir)



}
