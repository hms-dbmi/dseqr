#' Logic Bulk Data page
#' @export
#' @keywords internal
bulkPage <- function(input, output, session, data_dir, sc_dir, bulk_dir, indices_dir) {

  new_dataset <- reactiveVal()
  msg_quant <- reactiveVal()

  eset <- reactive(readRDS(file.path(bulkForm$dataset_dir(), 'eset.rds')))


  explore_eset <- exploreEset(eset = eset,
                              dataset_dir = bulkForm$dataset_dir,
                              explore_pdata = dsExploreTable$pdata,
                              numsv = bulkForm$numsv_r,
                              svobj = bulkForm$svobj_r)


  bulkForm <- callModule(bulkForm, 'form',
                         data_dir = data_dir,
                         sc_dir = sc_dir,
                         bulk_dir = bulk_dir,
                         new_dataset = new_dataset,
                         msg_quant = msg_quant,
                         new_anal = bulk_anal,
                         explore_eset = explore_eset,
                         enable_sva = dsExploreTable$enable_sva)

  callModule(bulkMDSplotly, 'mds_plotly',
             explore_eset = explore_eset,
             dataset_name = bulkForm$dataset_name,
             numsv = bulkForm$numsv_r)


  # toggle tables
  observe({
    toggle('quant_table_container', condition = bulkForm$show_quant())
    toggle('anal_table_container', condition = bulkForm$show_anal())
  })

  # toggle plots
  sel_genes <- reactive(length(bulkForm$explore_genes() > 0))

  observe({
    toggle('mds_plotly_container', condition = !sel_genes() & !bulkForm$show_dtangle())
    toggle('gene_plotly_container', condition = sel_genes() & !bulkForm$show_dtangle())
    toggle('cells_plotly_container', condition = bulkForm$show_dtangle())
  })

  dsQuantTable <- callModule(bulkQuantTable, 'quant',
                             fastq_dir = bulkForm$fastq_dir,
                             labels = bulkForm$quant_labels,
                             paired = bulkForm$paired)


  dsExploreTable <- callModule(bulkExploreTable, 'explore',
                               eset = eset,
                               labels = bulkForm$explore_labels,
                               data_dir = data_dir,
                               dataset_dir = bulkForm$dataset_dir,
                               dataset_name = bulkForm$dataset_name,
                               svobj_r = bulkForm$svobj_r,
                               numsv_r = bulkForm$numsv_r)


  callModule(bulkGenePlotly, 'gene_plotly',
             eset = explore_eset,
             explore_genes = bulkForm$explore_genes,
             dataset_name = bulkForm$dataset_name)

  callModule(bulkCellsPlotly, 'cells_plotly',
             dtangle_est = bulkForm$dtangle_est,
             pdata = dsExploreTable$pdata,
             dataset_name = bulkForm$dataset_name)

  observe({
    msg_quant(dsQuantTable$valid_msg())
  })



  # run quantification for quant dataset
  observeEvent(bulkForm$run_quant(), {
    # disable inputs
    shinyjs::disable(selector = 'input')

    # setup
    pdata <- dsQuantTable$pdata()
    dataset_name <- bulkForm$dataset_name()
    fastq_dir <- bulkForm$fastq_dir()

    # Create a Progress object
    progress <- Progress$new(session, min=0, max = nrow(pdata)+1)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())

    # Create a callback function to update progress.
    progress$set(message = "Quantifying files", value = 0)
    updateProgress <- function(amount = NULL, detail = NULL) {
      progress$inc(amount = amount, detail = detail)
    }

    # setup bulk
    pdata <- dsQuantTable$pdata()
    paired <- bulkForm$paired()


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

    # trigger to update rest of app
    new_dataset(dataset_name)

    # re-enable inputs
    shinyjs::enable(selector = 'input')
    progress$inc(1)
  })


  return(list(
    new_dataset = new_dataset
  ))
}


#' Logic for Bulk Data MDS plotly
#' @export
#' @keywords internal
bulkMDSplotly <- function(input, output, session, dataset_name, explore_eset, numsv) {

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

  base_fname <- reactive(paste0(dataset_name(), '_', numsv(), 'SV'))

  data_fun <- function(file) {
    #go to a temp dir to avoid permission issues
    owd <- setwd(tempdir())
    on.exit(setwd(owd))

    base_fname <- base_fname()
    eset <- explore_eset()
    vsd <- Biobase::assayDataElement(eset, 'vsd')
    pdata <- Biobase::pData(eset)

    vsd_fname <- paste0(base_fname, '.csv')
    pdata_fname <- paste0(base_fname, '_pdata.csv')
    write.csv(vsd, vsd_fname)
    write.csv(pdata, pdata_fname)

    #create the zip file
    zip(file, c(vsd_fname, pdata_fname))
  }

  fname_fun <- function() {paste0(base_fname(), '_', Sys.Date(), '.zip')}


  # not currently downloadable so use default fname_fun and data_fun
  callModule(downloadablePlotly, 'plotly',
             plotly_fun = plotly_fun,
             data_fun = data_fun,
             fname_fun = fname_fun,
             title = 'Download adjusted data')
}


#' Logic for Bulk Data gene plotly
#' @export
#' @keywords internal
bulkGenePlotly <- function(input, output, session, eset, explore_genes, dataset_name) {

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


#' Logic for Bulk Data cell type deconvolution plotly
#' @export
#' @keywords internal
bulkCellsPlotly <- function(input, output, session, dtangle_est, pdata, dataset_name) {

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


#' Logic for Bulk Data form
#' @export
#' @keywords internal
bulkForm <- function(input, output, session, data_dir, sc_dir, bulk_dir, new_dataset, msg_quant, new_anal, explore_eset, enable_sva) {

  dataset <- callModule(bulkDataset, 'selected_dataset',
                        data_dir = data_dir,
                        sc_dir = sc_dir,
                        bulk_dir = bulk_dir,
                        new_dataset = new_dataset,
                        explore_eset = explore_eset)


  # show quant, anals or neither
  show_quant <- reactive({
    dataset$dataset_name() != '' && dataset$is.create() && isTruthy(dataset$fastq_dir())
  })

  show_anal <- reactive({
    dataset$dataset_name() != '' && !dataset$is.create()
  })

  observe({
    toggle('quant_dataset_panel', condition = show_quant())
    toggle('anal_dataset_panel', condition = show_anal())
  })

  quant <- callModule(bulkFormQuant, 'quant_form',
                      fastq_dir = dataset$fastq_dir,
                      error_msg = msg_quant)

  anal <- callModule(bulkFormAnal, 'anal_form',
                     data_dir = data_dir,
                     dataset_dir = dataset$dataset_dir,
                     dataset_name = dataset$dataset_name,
                     explore_eset = explore_eset,
                     numsv_r = dataset$numsv_r,
                     svobj_r = dataset$svobj_r,
                     enable_sva = enable_sva)



  return(list(
    fastq_dir = dataset$fastq_dir,
    paired = quant$paired,
    quant_labels = quant$labels,
    anal_labels = anal$labels,
    explore_labels = anal$explore_labels,
    explore_genes = anal$explore_genes,
    run_quant = quant$run_quant,
    dataset_name = dataset$dataset_name,
    dataset_dir = dataset$dataset_dir,
    show_quant = show_quant,
    show_anal = show_anal,
    is.cellranger = dataset$is.cellranger,
    is.explore = anal$is.explore,
    dtangle_est = dataset$dtangle_est,
    show_dtangle = dataset$show_dtangle,
    numsv_r = dataset$numsv_r,
    svobj_r = dataset$svobj_r
  ))

}


#' Logic for selected dataset part of bulkFrom
#' @export
#' @keywords internal
bulkDataset <- function(input, output, session, sc_dir, bulk_dir, data_dir, new_dataset, explore_eset) {

  # get directory with fastqs
  roots <- c('bulk' = bulk_dir)
  shinyFiles::shinyDirChoose(input, "new_dataset_dir", roots = roots)

  # only show nsv/dtangle toggle if existing dataset
  observe({
    toggleSelectizeButtons('dataset_name',
                           button_ids = c('show_nsv', 'show_dtangle'),
                           condition = dataset_exists())
  })


  datasets <- reactive({
    new_dataset()
    load_bulk_datasets(data_dir)
  })

  # is the dataset a quantified one?
  is.create <- reactive({
    dataset_name <- input$dataset_name
    datasets <- datasets()
    req(dataset_name)

    !dataset_name %in% datasets$dataset_name
  })

  dataset_exists <- reactive(isTruthy(input$dataset_name) & !is.create())

  observe({
    req(datasets())
    updateSelectizeInput(session, 'dataset_name', choices = rbind(rep(NA, 5), datasets()), server = TRUE)
  })

  # open selector if creating
  observe({
    req(is.create())
    shinyjs::click('new_dataset_dir')
  })

  # directory with fastq files for quantificant
  fastq_dir <- reactive({
    req(!'integer' %in% class(input$new_dataset_dir))
    dir <- shinyFiles::parseDirPath(roots, input$new_dataset_dir)
    as.character(dir)
  })

  # directory to existing dataset for exploration
  dataset_dir <- reactive({
    req(dataset_exists())
    dataset_name <- input$dataset_name
    datasets <- datasets()
    req(dataset_name, datasets)
    dir <- datasets[datasets$dataset_name == dataset_name, 'dataset_dir']
    file.path(data_dir, dir)
  })

  numsv_path <- reactive(file.path(dataset_dir(), 'numsv.rds'))
  svobj_path <- reactive(file.path(dataset_dir(), 'svobj.rds'))

  # initialize svobj
  svobj_r <- reactiveVal()

  observe({
    svobj <- NULL
    svobj_path <- svobj_path()

    if (file.exists(svobj_path)) {
      svobj <- readRDS(svobj_path)
    }

    svobj_r(svobj)
  })

  # initialize numsv
  numsv_r <- reactiveVal()
  maxsv_r <- reactiveVal()

  # initialize selected and max number of svs
  observe({
    numsv <- maxsv <- 0
    svobj <- svobj_r()
    if (!is.null(svobj$n.sv)) maxsv <- svobj$n.sv

    numsv_path <- numsv_path()
    if (file.exists(numsv_path)) numsv <- readRDS(numsv_path)
    else saveRDS(0, numsv_path)

    maxsv_r(maxsv)
    numsv_r(numsv)
  })

  # update number of surrogate variables slider
  observe({
    updateSliderInput(session, 'selected_nsv', value = numsv_r(), min = 0, max = maxsv_r())
  })

  observeEvent(input$selected_nsv, {
    numsv_r(input$selected_nsv)
    saveRDS(input$selected_nsv, numsv_path())
  }, ignoreInit = TRUE)


  # update button icon with selected number of surrogate variables
  observe({
    updateActionButton(session, 'show_nsv', label = htmltools::doRenderTags(tags$span(input$selected_nsv, class='fa fa-fw')))
  })

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
    dataset_name = reactive(input$dataset_name),
    dataset_dir = dataset_dir,
    is.create = is.create,
    dtangle_est = dtangleForm$dtangle_est,
    show_dtangle = show_dtangle,
    selected_nsv = reactive(input$selected_nsv),
    numsv_r = numsv_r,
    svobj_r = svobj_r
  ))

}


#' Logic for dataset quantification part of bulkForm
#' @export
#' @keywords internal
bulkFormQuant <- function(input, output, session, fastq_dir, error_msg) {


  paired <- callModule(bulkEndType, 'end_type',
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

  quantModal <- function(paired) {
    end_type <- ifelse(paired, 'pair ended', 'single ended')

    UI <- withTags({
      dl(
        dt('Are you sure?'),
        dd('This will take a while.')
      )
    })

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
    showModal(quantModal(paired()))
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

#' Logic for end type selection is bulkFormQuant
#' @export
#' @keywords internal
bulkEndType <- function(input, output, session, fastq_dir) {


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

#' Logic for differential expression analysis part of bulkForm
#' @export
#' @keywords internal
bulkFormAnal <- function(input, output, session, data_dir, dataset_name, dataset_dir, explore_eset, numsv_r, svobj_r, enable_sva) {

  observe({
    toggleState('run_sva', condition = enable_sva())
  })

  # Group labels section
  # ---

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



  # run surrogate variable analysis
  observeEvent(input$run_sva, {


    eset <- explore_eset()
    pdata <- Biobase::pData(eset)
    group <- pdata$group
    req(length(unique(group)) > 1)

    # remove previously adjusted data
    remove_dataset_files(dataset_dir(), exclude = '_0svs.rds$')

    mods <- get_mods(eset)
    rna_seq <- 'lib.size' %in% colnames(pdata)
    svobj <- run_sva(mods, eset, rna_seq = rna_seq)

    # add row names so that can check sync during dataset switch
    row.names(svobj$sv) <- colnames(eset)

    # save current pdata_explore  so that can tell if changed
    file.copy(file.path(dataset_dir(), 'pdata_explore.rds'),
              file.path(dataset_dir(), 'pdata_explore_prev.rds'), overwrite = TRUE)


    # update saved svobj
    saveRDS(svobj, file.path(dataset_dir(), 'svobj.rds'))

    svobj_r(svobj)
    numsv_r(0)

  })

  # Gene choices
  # ---
  observeEvent(dataset_name(), {
    eset <- explore_eset()
    choices <- c(NA, row.names(eset))
    updateSelectizeInput(session, 'explore_genes', choices = choices, server = TRUE)
  })

  # Bulk anal
  # ---
  pdata <- reactive(Biobase::pData(explore_eset()))


  bulkAnal <- callModule(bulkAnal, 'ds',
                         pdata = pdata,
                         dataset_name = dataset_name,
                         eset = explore_eset,
                         svobj = svobj_r,
                         numsv = numsv_r,
                         dataset_dir = dataset_dir)


  return(list(
    labels = labels,
    explore_labels = explore_labels,
    explore_genes = reactive(input$explore_genes),
    run_sva = reactive(input$run_sva)
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

    # get normalized/adjusted values
    adj <- Biobase::assayDataElement(eset, 'adjusted')

    # common genes only
    commongenes <- intersect (rownames(adj), rownames(scseq))
    adj <- adj[commongenes, ]
    scseq <- scseq[commongenes, ]

    # quantile normalize scseq and rnaseq dataset
    progress$set(value = 2)
    y <- cbind(as.matrix(scseq[['SCT']]@data), adj)

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
bulkQuantTable <- function(input, output, session, fastq_dir, labels, paired) {

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
bulkExploreTable <- function(input, output, session, eset, labels, data_dir, dataset_dir, dataset_name, svobj_r, numsv_r) {

  # colors
  group_colors <- RColorBrewer::brewer.pal(8, 'Set2')
  ncolors <- length(group_colors)

  # things user will update and return
  pdata_r <- reactiveVal()
  group_r <- reactiveVal()
  name_r <- reactiveVal()
  table_rendered <- reactiveVal()
  enable_sva <- reactiveVal()

  pdata_path <- reactive(file.path(dataset_dir(), 'pdata_explore.rds'))

  # initial sva btn status
  observe({
    enable <- TRUE
    prev_path <- file.path(dataset_dir(), 'pdata_explore_prev.rds')
    if (file.exists(prev_path)) {
      prev <- readRDS(prev_path)
      pdata <- readRDS(pdata_path())
      enable <- !identical(prev, pdata)
    }

    enable_sva(enable)
  })


  # pdata that gets returned with group column
  returned_pdata <- reactive({
    pdata_path <- pdata_path()
    pdata <- pdata_r()
    req(pdata, pdata_path)

    pdata$Group <- group_r()
    pdata$`Group name` <- name_r()

    # prevent overwriting saved pdata when switch between analyses
    if (file.exists(pdata_path)) {
      saved_pdata <- readRDS(pdata_path)
      req(all(row.names(pdata) == row.names(saved_pdata)))
    }

    saveRDS(pdata, pdata_path)
    return(pdata)
  })

  eset_pdata <- reactive({
    eset <- eset()
    req(eset)

    pdata <- Biobase::pData(eset) %>%
      tibble::add_column(Group = NA, 'Group name' = NA, Title = colnames(eset), .before = 1)

    return(pdata)
  })


  # reset when new eset or analysis name
  observe({
    pdata_path <- pdata_path()
    eset_pdata <- pdata <- eset_pdata()
    req(eset_pdata, pdata_path)

    group <- name <- rep(NA, nrow(eset_pdata))
    # load pdata from previous if available
    if (file.exists(pdata_path)) {
      saved_pdata <- pdata <- readRDS(pdata_path)

      group <- pdata$Group
      name  <- pdata$`Group name`
    }

    table_rendered(FALSE)

    pdata_r(pdata)
    group_r(group)
    name_r(name)
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
    rows  <- input$pdata_rows_selected
    req(group_name, length(rows))

    group <- group_r()
    name <- name_r()

    group_num <- length(unique(setdiff(group, NA))) + 1

    group[rows] <- group_num
    name[rows] <- group_name
    group_r(group)
    name_r(name)

    # remove dataset files since groups changed
    svobj_r(NULL)
    numsv_r(0)
    remove_dataset_files(dataset_dir())
    remove_bulk_anals(dataset_name(), data_dir)

    if (length(group_num) > 1) enable_sva(TRUE)
  })



  # click 'Reset'
  shiny::observeEvent(labels$reset(), {
    req(labels$reset())
    group <- group_r()
    clear <- rep(NA, length(group))
    group_r(clear)
    name_r(clear)

    # remove dataset files since groups changed
    svobj_r(NULL)
    numsv_r(0)
    remove_dataset_files(dataset_dir())
    remove_bulk_anals(dataset_name(), data_dir)
    enable_sva(FALSE)
  })



  return(list(
    pdata = returned_pdata,
    enable_sva = enable_sva
  ))

}

#' Get group levels for bulk data plots
#'
#' @param pdata Data.frame of phenotype data
#' @export
#'
#' @keywords internal
get_group_levels <- function(pdata) {
  group <- pdata$`Group name`
  group_order <- order(unique(pdata$Group))
  unique(group)[group_order]
}

#' Get group colors for bulk data plots
#'
#' @param group_levels result of \link{get_group_levels}
#' @export
#'
#' @keywords internal
get_group_colors <- function(group_levels) {
  RColorBrewer::brewer.pal(8, 'Set2')[seq_along(group_levels)]
}


#' Logic for bulk group analyses for Bulk, Drugs, and Pathways tabs
#' @export
#' @keywords internal
bulkAnal <- function(input, output, session, pdata, dataset_name, eset, numsv, svobj, dataset_dir, is_bulk = function()TRUE) {
  contrast_options <- list(render = I('{option: bulkContrastOptions, item: bulkContrastItem}'))


  # group levels used for selecting test and control groups
  group_levels <- reactive({
    req(is_bulk())
    get_group_levels(pdata())
  })

  group_colors <- reactive(get_group_colors(group_levels()))

  group_choices <- reactive({

    data.frame(
      name = group_levels(),
      value = group_levels(),
      color = group_colors(), stringsAsFactors = FALSE
    )
  })

  observe({
    updateSelectizeInput(session, 'contrast_groups', choices = group_choices(), server = TRUE, options = contrast_options)
  })

  full_contrast <- reactive(length(input$contrast_groups) == 2)

  anal_name <- reactive({
    req(full_contrast())
    groups <- input$contrast_groups
    return(paste0(groups[1], '_vs_', groups[2]))
  })

  # path to lmfit and drug query results
  numsv_str <- reactive(paste0(numsv(), 'svs'))

  lmfit_path <- reactive({
    req(is_bulk())
    lmfit_file <- paste0('lm_fit_', numsv_str(), '.rds')
    file.path(dataset_dir(), lmfit_file)
  })

  drug_paths <- reactive({
    suffix <- paste(anal_name(), numsv_str(), sep = '_')
    get_drug_paths(dataset_dir(), suffix)
  })

  goana_path <- reactive({
    fname <- paste0('goana_', anal_name(), '_', numsv_str(), '.rds')
    file.path(dataset_dir(), fname)
  })

  # do we have lm_fit and drug query results?
  saved_lmfit <- reactive(file.exists(lmfit_path()))
  saved_drugs <- reactive(file.exists(drug_paths()$cmap))

  # load lm_fit if saved or run limma
  disable_inputs <- ''
  lm_fit <- reactive({

    if (saved_lmfit()) {
      lm_fit <- readRDS(lmfit_path())

    } else {

      # visual that running
      toggleAll(disable_inputs)

      progress <- Progress$new(session, min=0, max = 3)
      on.exit(progress$close())

      # check for previous lm_fit
      dataset_dir <- dataset_dir()
      numsv <- numsv()

      # get what need
      eset <- eset()
      svobj <- svobj()
      req(eset, dataset_dir)

      prev_anal <- list(pdata = Biobase::pData(eset))

      # run differential expression
      progress$set(message = "Fitting limma model", value = 1)
      lm_fit <- run_limma(eset,
                          dataset_dir = dataset_dir,
                          svobj = svobj,
                          numsv = numsv,
                          prev_anal = prev_anal)

      # visual that done
      progress$inc(1)
      toggleAll(disable_inputs)
    }

    return(lm_fit)
  })

  drug_queries <- reactive({
    if (!full_contrast()) {
      res <- NULL
    } else if (saved_drugs()) {
      paths <- drug_paths()
      res <- lapply(paths, readRDS)

    } else {
      top_table <- top_table()
      res <- run_drug_queries(top_table, drug_paths(), session)
    }
    return(res)
  })

  # differential expression top table
  top_table <- reactive({
    req(full_contrast())
    lm_fit <- lm_fit()
    groups <- input$contrast_groups

    # loses sync when groups selected and change dataset
    req(all(groups %in% colnames(lm_fit$mod)))

    tt <- get_top_table(lm_fit, groups)
    tt[order(tt$P.Value), ]
  })

  # goana pathway result
  path_res <- reactive({
    req(full_contrast())
    goana_path <- goana_path()

    if (file.exists(goana_path)) {
      goana_res <- readRDS(goana_path)

    } else {
      lm_fit <- lm_fit()
      groups <- input$contrast_groups

      # loses sync when groups selected and change dataset
      req(all(groups %in% colnames(lm_fit$mod)))

      contrast <- paste0(groups[1], '-', groups[2])
      ebfit <- fit_ebayes(lm_fit, contrast)
      goana_res <- limma::goana(ebfit, species = 'Hs', geneid = 'ENTREZID')
      saveRDS(goana_res, goana_path)
    }

    return(goana_res)
  })


  # enable download
  observe({
    toggleState('download', condition = full_contrast())
  })

  dl_fname <- reactive({
    date <- paste0(Sys.Date(), '.zip')

    numsv_str <- paste0(numsv(), 'SV')
    paste('bulk', dataset_name(), anal_name(), numsv_str, date , sep='_')
  })

  data_fun <- function(file) {
    #go to a temp dir to avoid permission issues
    owd <- setwd(tempdir())
    on.exit(setwd(owd))

    tt_fname <- 'top_table.csv'
    go_fname <- 'goana.csv'

    write.csv(top_table(), tt_fname)
    write.csv(path_res(), go_fname)

    #create the zip file
    zip(file, c(tt_fname, go_fname))
  }

  output$download <- downloadHandler(
    filename = function() {
      dl_fname()
    },
    content = data_fun
  )

  return(list(
    name = anal_name,
    contrast_groups = reactive(input$contrast_groups),
    lm_fit = lm_fit,
    top_table = top_table,
    drug_queries = drug_queries,
    path_res = path_res
  ))
}



#' Logic to setup explore_eset for Bulk Data plots
#' @export
#' @keywords internal
exploreEset <- function(eset, dataset_dir, explore_pdata, numsv, svobj) {

  vsd_path <- reactive(file.path(dataset_dir(), 'vsd.rds'))
  adj_path <- reactive(file.path(dataset_dir(), paste0('adjusted_', numsv(), 'svs.rds')))
  keep_path <- reactive(file.path(dataset_dir(), paste0('iqr_keep_', numsv(), 'svs.rds')))


  norm_eset <- reactive({
    # pdata and eset lose sync when switch datasets
    eset <- eset()
    pdata <- explore_pdata()

    # need that pdata_explore has more than two groups
    keep <- row.names(pdata)[!is.na(pdata$Group)]
    pdata <- pdata[keep, ]
    req(length(unique(pdata$Group)) > 1)

    # subset eset and add explore_pdata
    eset <- eset[, keep]
    pdata$group <- pdata$`Group name`
    Biobase::pData(eset) <- pdata

    # rlog normalize
    eset <- add_vsd(eset, vsd_path = vsd_path())
    return(eset)
  })


  # explore_eset used for all plots
  explore_eset <- reactive({
    eset <- norm_eset()
    numsv <- numsv()

    # adjust for pairs/surrogate variables
    svobj <- svobj()

    # can lose sync when switching datasets
    if (!is.null(svobj)) {
      req(numsv <= svobj$n.sv)
      req(row.names(svobj$sv) == colnames(eset))
    }
    eset <- add_adjusted(eset, svobj, numsv, adj_path = adj_path())

    # use SYMBOL as annotation
    # keep unique symbol based on row IQRs
    eset <- iqr_replicates(eset, keep_path = keep_path())

    return(eset)
  })
  #
  return(explore_eset)
}
