#' Logic for Bulk Data Tab
#'
#' to be called with \link[shiny]{callModule}
#'
#' @param input,output,session standard shiny module boilerplate
#' @param data_dir path to folder with application name
#' @param sc_dir sub folder of \code{data_dir} where single-cell data is stored
#' @param bulk_dir sub folder of \code{data_dir} where bulk data is stored
#' @inheritParams run_drugseqr
#'
#' @return list with reactive \code{new_dataset} that is triggered with a new
#'   bulk dataset is added.
#'
#' @export
bulkPage <- function(input, output, session, data_dir, sc_dir, bulk_dir, indices_dir) {

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
                         msg_quant = msg_quant,
                         explore_eset = explore_eset,
                         pdata = dsQuantTable$pdata,
                         indices_dir = indices_dir)

  # mds plotly with different orientations for mobile/desktop

  bulkMDS <- callModule(bulkMDS, 'bulk_mds',
                        explore_eset = explore_eset)

  callModule(bulkMDSplotly, 'mds_plotly_unadjusted',
             explore_eset = explore_eset,
             dataset_name = bulkForm$dataset_name,
             numsv = bulkForm$numsv_r,
             mds = bulkMDS$mds,
             group_colors = bulkMDS$group_colors,
             adjusted = FALSE)

  callModule(bulkMDSplotly, 'mds_plotly_adjusted',
             explore_eset = explore_eset,
             dataset_name = bulkForm$dataset_name,
             numsv = bulkForm$numsv_r,
             mds = bulkMDS$mds,
             group_colors = bulkMDS$group_colors,
             adjusted = TRUE)


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

  up_annot <- callModule(bulkAnnot, 'anal',
                         dataset_name = bulkForm$dataset_name,
                         annot = dsExploreTable$pdata)


  dsExploreTable <- callModule(bulkExploreTable, 'explore',
                               eset = eset,
                               up_annot = up_annot,
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



  return(list(
    new_dataset = bulkForm$new_dataset
  ))
}


#' Logic for Bulk Data MDS data
#'
#' @keywords internal
#' @noRd
bulkMDS <- function(input, output, session, explore_eset) {

  group <- reactive({
    eset <- explore_eset()
    pdata <- Biobase::pData(eset)

    # setup group factor and colors
    group <- pdata$`Group name`
    group_order <- order(unique(pdata$Group))
    group_levels <- unique(group)[group_order]
    factor(group, levels = group_levels)
  })


  group_colors <- reactive({
    get_palette(levels(group()))
  })


  mds <- reactive({
    eset <- explore_eset()
    vsd <- Biobase::assayDataElement(eset, 'vsd')
    adj <- Biobase::assayDataElement(eset, 'adjusted')
    get_mds(vsd, adj, group())
  })

  return(list(
    mds = mds,
    group_colors = group_colors
  ))
}


#' Logic for Bulk Data MDS plotly
#'
#' @keywords internal
#' @noRd
bulkMDSplotly <- function(input, output, session, explore_eset, dataset_name, numsv, mds, group_colors, adjusted) {

  # MDS plot
  plot <- reactive({
    mds <- mds()
    plotlyMDS(mds$scaling, mds$scaling_adj, group_colors = group_colors(), adjusted = adjusted)
  })

  base_fname <- reactive(paste0(dataset_name(), '_', numsv(), 'SV'))

  content <- function(file) {
    #go to a temp dir to avoid permission issues
    owd <- setwd(tempdir())
    on.exit(setwd(owd))

    base_fname <- base_fname()
    eset <- explore_eset()
    vsd <- Biobase::assayDataElement(eset, 'vsd')
    pdata <- Biobase::pData(eset)

    vsd_fname <- paste0(base_fname, '.csv')
    pdata_fname <- paste0(base_fname, '_pdata.csv')
    utils::write.csv(vsd, vsd_fname)
    utils::write.csv(pdata, pdata_fname)

    #create the zip file
    utils::zip(file, c(vsd_fname, pdata_fname))
  }

  filename <- function() {paste0(base_fname(), '_', Sys.Date(), '.zip')}


  # not currently downloadable so use default fname_fun and data_fun
  callModule(shinydlplot::downloadablePlotly, 'plotly',
             plot = plot,
             content = content,
             filename = filename,
             title = 'Download adjusted data')
}


#' Logic for Bulk Data gene plotly
#'
#' @keywords internal
#' @noRd
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

  plot <- reactive({

    args <- boxplotly_args()

    boxPlotly(df = args$df,
              boxgap = args$boxgap,
              boxgroupgap = args$boxgroupgap,
              plot_fname = args$plot_fname,
              ytitle = 'Normalized Expression',
              xtitle = 'Gene')
  })

  filename <- function() {
    paste(dataset_name(), '_', paste(explore_genes(), collapse = '_'), '_', Sys.Date(), ".csv", sep = "")
  }


  content <- function(file) {
    args <- boxplotly_args()
    df <- args$df
    df$color <- NULL
    colnames(df) <- c('Sample', 'Gene', 'Normalized Expression', 'Group')
    utils::write.csv(df, file, row.names = FALSE)
  }

  callModule(shinydlplot::downloadablePlotly,
             'plotly',
             plot = plot,
             filename = filename,
             content = content)



}


#' Logic for Bulk Data cell type deconvolution plotly
#'
#' @keywords internal
#' @noRd
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

  plot <- reactive({

    args <- boxplotly_args()

    boxPlotlyCells(df = args$df,
                   boxgap = args$boxgap,
                   boxgroupgap = args$boxgroupgap,
                   plot_fname = args$plot_fname,
                   ytitle = 'Estimated Proportion',
                   xtitle = 'Cluster')
  })

  filename <- function() {
    clusters <- colnames(dtangle_est())
    clusters <- gsub(' ', '', clusters)
    paste(dataset_name(), '_', paste(clusters, collapse = '_'), '_', Sys.Date(), ".csv", sep = "")
  }


  content <- function(file) {
    args <- boxplotly_args()
    df <- args$df
    df$color <- NULL
    colnames(df) <- c('Sample', 'Group', 'Proportion', 'Cluster')
    utils::write.csv(df, file, row.names = FALSE)
  }

  callModule(shinydlplot::downloadablePlotly, 'plotly', plot = plot, filename = filename, content = content)


}


#' Logic for Bulk Data form
#'
#' @keywords internal
#' @noRd
bulkForm <- function(input, output, session, data_dir, sc_dir, bulk_dir, msg_quant, explore_eset, pdata, indices_dir) {

  dataset <- callModule(bulkDataset, 'selected_dataset',
                        data_dir = data_dir,
                        sc_dir = sc_dir,
                        bulk_dir = bulk_dir,
                        new_dataset = quant$new_dataset,
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
                      data_dir = data_dir,
                      error_msg = msg_quant,
                      dataset_name = dataset$dataset_name,
                      pdata = pdata,
                      fastq_dir = dataset$fastq_dir,
                      indices_dir = indices_dir)


  anal <- callModule(bulkFormAnal, 'anal_form',
                     data_dir = data_dir,
                     dataset_dir = dataset$dataset_dir,
                     dataset_name = dataset$dataset_name,
                     explore_eset = explore_eset,
                     numsv_r = dataset$numsv_r,
                     svobj_r = dataset$svobj_r)



  return(list(
    fastq_dir = dataset$fastq_dir,
    paired = quant$paired,
    quant_labels = quant$labels,
    new_dataset = quant$new_dataset,
    anal_labels = anal$labels,
    explore_genes = anal$explore_genes,
    dataset_name = dataset$dataset_name,
    dataset_dir = dataset$dataset_dir,
    show_quant = show_quant,
    show_anal = show_anal,
    is.explore = anal$is.explore,
    dtangle_est = dataset$dtangle_est,
    show_dtangle = dataset$show_dtangle,
    numsv_r = dataset$numsv_r,
    svobj_r = dataset$svobj_r
  ))

}


#' Logic for selected dataset part of bulkFrom
#'
#' @keywords internal
#' @noRd
bulkDataset <- function(input, output, session, sc_dir, bulk_dir, data_dir, new_dataset, explore_eset) {

  # get directory with fastqs
  roots <- c('bulk' = bulk_dir)
  shinyFiles::shinyDirChoose(input, "new_dataset_dir", roots = roots)

  # only show nsv/dtangle toggle if existing dataset
  observe({
    shinypanel::toggleSelectizeButtons('dataset_name',
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
    req(!methods::is(input$new_dataset_dir, 'integer'))
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
    new_dataset()
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
#'
#' @keywords internal
#' @noRd
bulkFormQuant <- function(input, output, session, error_msg, dataset_name, pdata, fastq_dir, data_dir, indices_dir) {
  quant_inputs <- c('end_type', 'pair', 'rep', 'reset', 'run_quant')

  paired <- bulkEndType(input, output, session, fastq_dir)

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
      tags$dl(
        tags$dt('Are you sure?'),
        tags$dd('This will take a while.')
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

  # Show confirm modal when click Run Qunatification
  observeEvent(input$run_quant, {
    showModal(quantModal(paired()))
  })

  new_dataset <- reactiveVal()

  # run quantification upon confirmation
  observeEvent(input$confirm, {
    removeModal()

    # check for index
    index_path <- rkal::get_kallisto_index(indices_dir)
    if(!length(index_path)) {
      error_msg('No kallisto index. See github README.')
      return(NULL)
    }

    # disable inputs
    disableAll(quant_inputs)

    # setup
    pdata <- pdata()
    paired <- paired()
    dataset_name <- dataset_name()
    fastq_dir <- fastq_dir()

    # Create a Progress object
    progress <- Progress$new(session, min=0, max = nrow(pdata)+1)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())

    # Create a callback function to update progress.
    progress$set(message = "Quantifying files", value = 0)
    updateProgress <- function(amount = NULL, detail = NULL) {
      progress$inc(amount = amount, detail = detail)
    }

    # quantification
    rkal::run_kallisto_bulk(
      indices_dir = indices_dir,
      data_dir = fastq_dir,
      quant_meta = pdata,
      paired = paired,
      updateProgress = updateProgress)

    # generate eset and save
    progress$set(message = 'Annotating dataset')
    eset <- rkal::load_seq(fastq_dir)

    # trigger to update rest of app
    new_dataset(dataset_name)

    # re-enable inputs
    enableAll(quant_inputs)
    progress$inc(1)
  })

  return(list(
    new_dataset = new_dataset,
    paired = paired,
    labels = list(
      reset = reset,
      pair = pair,
      rep = rep
    )
  ))
}

#' Logic for end type selection is bulkFormQuant
#'
#' @keywords internal
#' @noRd
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
    fastq_id1s <- rkal::get_fastq_id1s(file.path(fastq_dir, fastqs))
    rkal::detect_paired(fastq_id1s)
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
#'
#' @keywords internal
#' @noRd
bulkFormAnal <- function(input, output, session, data_dir, dataset_name, dataset_dir, explore_eset, numsv_r, svobj_r) {


  # run surrogate variable analysis if required
  observeEvent(explore_eset(), {

    eset <- explore_eset()
    pdata <- Biobase::pData(eset)
    req(eset)

    prev_path <- file.path(dataset_dir(), 'pdata_explore_prev.rds')
    pdata_path <- file.path(dataset_dir(), 'pdata_explore.rds')

    if (file.exists(prev_path)) {
      prev <- readRDS(prev_path)
      saved <- readRDS(pdata_path)
      changed <- check_bulk_changed(prev, saved)
      req(changed)
    }

    group <- pdata$group
    req(length(unique(group)) > 1)

    # Create a Progress object
    progress <- Progress$new(session, min=0, max = 2)
    on.exit(progress$close())
    progress$set(message = "Running SVA", value = 1)

    mods <- crossmeta::get_sva_mods(eset)
    svobj <- crossmeta::run_sva(mods, eset)

    # add row names so that can check sync during dataset switch
    row.names(svobj$sv) <- colnames(eset)

    # save current pdata_explore  so that can tell if changed
    file.copy(pdata_path, prev_path, overwrite = TRUE)


    # update saved svobj
    saveRDS(svobj, file.path(dataset_dir(), 'svobj.rds'))

    svobj_r(svobj)
    numsv_r(NULL)

    progress$set(value = 2)

  })

  # Gene choices
  # ---
  observe({
    dataset_name()
    eset <- explore_eset()
    choices <- c(NA, row.names(eset))
    updateSelectizeInput(session, 'explore_genes', choices = choices, server = TRUE)
  })

  # Bulk anal
  # ---
  pdata <- reactive({
    req(explore_eset())
    Biobase::pData(explore_eset())
  })


  bulkAnal <- callModule(bulkAnal, 'ds',
                         pdata = pdata,
                         dataset_name = dataset_name,
                         eset = explore_eset,
                         svobj = svobj_r,
                         numsv = numsv_r,
                         dataset_dir = dataset_dir)


  return(list(
    explore_genes = reactive(input$explore_genes)
  ))
}

#' Logic for deconvolution form
#'
#' @keywords internal
#' @noRd
dtangleForm <- function(input, output, session, show_dtangle, new_dataset, sc_dir, bulk_dir, explore_eset, dataset_dir, dataset_name) {
  include_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))
  input_ids <- c('include_clusters', 'dtangle_dataset', 'submit_dtangle')

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

    # exclude individual without scseq (e.g. folder with fastq.gz files only)
    has.scseq <- sapply(individual, function(ind) any(list.files(file.path(sc_dir, ind)) == 'scseq.rds'))
    individual <- individual[unlist(has.scseq)]
    return(individual)
  })



  # update reference dataset choices
  observe({
    ref_anals <- ref_anals()
    req(ref_anals)
    updateSelectizeInput(session, 'dtangle_dataset', choices = c('', ref_anals))
  })

  annot <- reactive({
    anal_name <- input$dtangle_dataset
    req(anal_name)
    annot_path <- scseq_part_path(sc_dir, anal_name, 'annot')
    readRDS(annot_path)
  })

  # update exclude cluster choices
  include_choices <- reactive({
    clusters <- annot()
    dataset_dir <- file.path(sc_dir, input$dtangle_dataset)
    get_cluster_choices(clusters, dataset_dir = dataset_dir)
  })

  observe({
    choices <- include_choices()
    updateSelectizeInput(session, 'include_clusters', choices = choices, options = include_options, server = TRUE)
  })

  # scseq for deconvolution
  scseq <- reactive({
    dataset_name <- input$dtangle_dataset
    load_scseq(file.path(sc_dir, dataset_name))
  })

  observeEvent(input$submit_dtangle, {

    # disable inputs
    disableAll(input_ids)

    # Create a Progress object
    progress <- Progress$new(session, min=0, max = 4)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())

    progress$set(message = "Deconvoluting", value = 0)
    progress$set(value = 1)

    dtangle_dataset <- input$dtangle_dataset

    # get names of clusters
    include_clusters <- input$include_clusters
    include_choices <- include_choices()

    # if select none deconvolute using all clusters
    if (!length(include_clusters))
      include_clusters <- as.character(include_choices$value)

    include_names <- include_choices$name[as.numeric(include_clusters)]

    eset <- explore_eset()
    pdata <- Biobase::pData(eset)

    # subset to selected clusters
    scseq <- scseq()
    scseq <- scseq[, scseq$cluster %in% include_clusters]

    # get normalized/adjusted values
    adj <- Biobase::assayDataElement(eset, 'adjusted')

    # common genes only
    commongenes <- intersect(rownames(adj), rownames(scseq))
    adj <- adj[commongenes, ]
    scseq <- scseq[commongenes, ]

    # quantile normalize scseq and rnaseq dataset
    progress$set(value = 2)
    y <- cbind(as.matrix(SingleCellExperiment::logcounts(scseq)), adj)

    y <- limma::normalizeBetweenArrays(y)
    y <- t(y)

    progress$set(value = 3)
    # indicies for cells in each included cluster
    pure_samples <- list()
    for (i in seq_along(include_clusters))
      pure_samples[[include_names[i]]] <- which(scseq$cluster == include_clusters[i])

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
    enableAll(input_ids)
    progress$set(value = 4)
  })


  return(list(
    dtangle_est = dtangle_est
  ))
}


#' Logic for dataset quantification table
#'
#' @keywords internal
#' @noRd
bulkQuantTable <- function(input, output, session, fastq_dir, labels, paired) {

  # things user will update and return
  pdata_r <- reactiveVal()
  pairs_r <- reactiveVal()
  reps_r <- reactiveVal()
  valid_msg <- reactiveVal()
  is_rendered <- reactiveVal(FALSE)


  # colors
  background <- 'url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAAPklEQVQoU43Myw0AIAgEUbdAq7VADCQaPyww55dBKyQiHZkzBIwQLqQzCk9E4Ytc6KEPMnTBCG2YIYMVpHAC84EnVbOkv3wAAAAASUVORK5CYII=) repeat'

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

    pdata$Pair <- pdata$Replicate <- NA_character_
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
    rep_nums <- as.numeric(rep_nums)
    rep_colors <- get_palette(rep_nums)
    for (rep_num in rep_nums) {
      color <- rep_colors[rep_num]
      rows <- which(reps == rep_num)
      pdata[rows, 'Replicate'] <- paste('<div style="background-color:', color, ';"></div>')
    }

    # update pdata Pair column
    if (paired()) {
      pair_nums <- sort(unique(setdiff(pairs, NA)))
      pair_nums <- as.numeric(pair_nums)
      pair_colors <- get_palette(pair_nums)
      for (pair_num in pair_nums) {
        color <- pair_colors[pair_num]
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
    msg <- rkal::validate_pairs(pairs, rows, reps)
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
    msg <- rkal::validate_reps(pairs, rows, reps)
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
#'
#' @keywords internal
#' @noRd
bulkExploreTable <- function(input, output, session, eset, up_annot, data_dir, dataset_dir, dataset_name, svobj_r, numsv_r) {

  # things user will update and return
  pdata_r <- reactiveVal()
  pdata_path <- reactive(file.path(dataset_dir(), 'pdata_explore.rds'))

  eset_pdata <- reactive({
    eset <- eset()
    req(eset)

    pdata <- Biobase::pData(eset) %>%
      tibble::add_column(Group = NA, 'Group name' = NA, Title = colnames(eset),  Pair = NA, .before = 1)

    return(pdata)
  })

  # update when upload new annotation
  observeEvent(up_annot(), {
    up <- up_annot()
    pdata_path <- pdata_path()
    req(up)


    if (file.exists(pdata_path)) {
      prev <- readRDS(pdata_path)
      changed <- check_bulk_changed(prev, up)

      if (changed) {
        remove_dataset_files(dataset_dir())
        svobj_r(NULL)
        numsv_r(NULL)
      }
    }

    pdata_r(up)
    saveRDS(up, pdata_path)
  })


  # reset when new eset
  observe({
    pdata_path <- pdata_path()
    pdata <- eset_pdata()
    req(pdata, pdata_path)

    # load pdata from previous if available
    if (file.exists(pdata_path)) {
      pdata <- readRDS(pdata_path)

      # TODO remove
      # need for legacy purposes
      pair <- pdata$Pair
      if (is.null(pair)) pair <- pdata$pair
      if (is.null(pair)) pair <- NA
      pdata$Pair <- NULL
      pdata <- tibble::add_column(pdata, Pair = pair, .after = 'Title')
    }

    pdata_r(pdata)
  })



  # format pdata for table
  html_pdata <- reactive({

    # things that trigger update
    pdata <- pdata_r()
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

    group_colors <- get_palette(group_nums)
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

  # redraw table when new pdata
  output$pdata <- DT::renderDataTable({

    pdata <- html_pdata()
    hide_target <- which(colnames(pdata) %in% c('lib.size', 'norm.factors', 'pair')) - 1

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


  return(list(
    pdata = pdata_r
  ))

}

#' Logic for downloading and uploading bulk annotation
#'
#' @keywords internal
#' @noRd
bulkAnnot <- function(input, output, session, dataset_name, annot) {


  observeEvent(input$click_up, {
    shinyjs::click('up_annot')
    error_msg(NULL)
  })

  observeEvent(input$click_dl, {
    shinyjs::click('dl_annot')
    error_msg(NULL)
  })

  fname <- reactive(paste0(dataset_name(), '_annot.csv'))

  output$dl_annot <- downloadHandler(
    filename = fname,
    content = function(con) {

      utils::write.csv(format_dl_annot(annot()), con, row.names = FALSE)
    }
  )

  # uploaded annotation
  up_annot <- reactiveVal()
  error_msg <- reactiveVal()

  # reset when change dataset
  observe({
    dataset_name()
    up_annot(NULL)
  })

  observe({
    msg <- error_msg()
    html('error_msg', html = msg)
    toggleClass('validate-up', 'has-error', condition = isTruthy(msg))
  })

  observeEvent(input$up_annot, {
    ref <- annot()
    req(ref)

    infile <- input$up_annot
    if (!isTruthy(infile)){
      res <- msg <- NULL

    } else {
      res <- utils::read.csv(infile$datapath, check.names = FALSE, stringsAsFactors = FALSE)
      msg <- validate_up_annot(res, annot())

      res <- if (is.null(msg)) format_up_annot(res, ref) else NULL
    }

    error_msg(msg)
    up_annot(res)
  })

  return(up_annot)

}


#' Logic for bulk group analyses for Bulk, Drugs, and Pathways tabs
#'
#' @keywords internal
#' @noRd
bulkAnal <- function(input, output, session, pdata, dataset_name, eset, numsv, svobj, dataset_dir, is_bulk = function()TRUE) {
  contrast_options <- list(render = I('{option: bulkContrastOptions, item: bulkContrastItem}'))
  input_ids <- c('download', 'contrast_groups')


  # group levels used for selecting test and control groups
  group_levels <- reactive({
    req(is_bulk())
    get_group_levels(pdata())
  })

  group_colors <- reactive(get_palette(group_levels()))

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

  go_path <- reactive({
    fname <- paste0('go_', anal_name(), '_', numsv_str(), '.rds')
    file.path(dataset_dir(), fname)
  })

  goana_path <- reactive({
    fname <- paste0('goana_', anal_name(), '_', numsv_str(), '.rds')
    file.path(dataset_dir(), fname)
  })

  kegg_path <- reactive({
    fname <- paste0('kegg_', anal_name(), '_', numsv_str(), '.rds')
    file.path(dataset_dir(), fname)
  })

  kegga_path <- reactive({
    fname <- paste0('kegga_', anal_name(), '_', numsv_str(), '.rds')
    file.path(dataset_dir(), fname)
  })

  # do we have drug query results?
  saved_drugs <- reactive(file.exists(drug_paths()$cmap))

  # load lm_fit if saved or run limma
  lm_fit <- reactive({

    if (file.exists(lmfit_path())) {
      lm_fit <- readRDS(lmfit_path())

    } else {

      # visual that running
      disableAll(input_ids)

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
      eset <- crossmeta::run_limma_setup(eset, prev_anal)
      lm_fit <- crossmeta::run_limma(eset,
                                     svobj = svobj,
                                     numsv = numsv,
                                     filter = FALSE)

      save_lmfit(lm_fit, dataset_dir, numsv = numsv)

      # visual that done
      progress$inc(1)
      enableAll(input_ids)
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
      tt <- top_table()
      if (!isTruthy(tt)) return(NULL)

      disableAll(input_ids)
      progress <- Progress$new(session, min = 0, max = 3)
      progress$set(message = "Querying drugs", value = 1)
      on.exit(progress$close())

      es <- load_drug_es()
      progress$inc(1)
      res <- run_drug_queries(tt, drug_paths(), es)
      progress$inc(1)
      enableAll(input_ids)
    }
    return(res)
  })

  # differential expression top table
  top_table <- reactive({
    if (!full_contrast()) return(NULL)
    lm_fit <- lm_fit()
    groups <- input$contrast_groups

    # loses sync when groups selected and change dataset
    if (!all(groups %in% colnames(lm_fit$mod))) return(NULL)

    crossmeta::get_top_table(lm_fit, groups)
  })

  # go/kegg pathway result
  path_res <- reactive({
    req(full_contrast())
    go_path <- go_path()
    kegg_path <- kegg_path()
    goana_path <- goana_path()
    kegga_path <- kegga_path()

    if (file.exists(kegga_path)) {
      res <- list(
        go = readRDS(go_path),
        kg = readRDS(kegg_path),
        goana = readRDS(goana_path),
        kegga = readRDS(kegga_path))

    } else {
      lm_fit <- lm_fit()

      # visual that running
      disableAll(input_ids)
      progress <- Progress$new(session, min=0, max = 2)
      on.exit(progress$close())
      progress$set(message = "Running pathway analysis", value = 1)
      groups <- input$contrast_groups

      # loses sync when groups selected and change dataset
      req(all(groups %in% colnames(lm_fit$mod)))

      groups <- make.names(groups)

      contrast <- paste0(groups[1], '-', groups[2])
      ebfit <- crossmeta::fit_ebayes(lm_fit, contrast)
      res <- get_path_res(ebfit,
                          go_path = go_path,
                          kegg_path = kegg_path,
                          goana_path = goana_path,
                          kegga_path = kegga_path)
      progress$inc(1)
      enableAll(input_ids)
    }

    return(res)
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
    go_fname <- 'cameraPR_go.csv'
    kg_fname <- 'cameraPR_kegg.csv'
    goana_fname <- 'goana.csv'
    kegga_fname <- 'kegga.csv'

    path_res <- path_res()
    utils::write.csv(top_table(), tt_fname)
    utils::write.csv(path_res$go, go_fname)
    utils::write.csv(path_res$kg, kg_fname)
    utils::write.csv(path_res$kegga, kegga_fname)
    utils::write.csv(path_res$goana, goana_fname)

    #create the zip file
    utils::zip(file, c(tt_fname, go_fname, kg_fname, goana_fname, kegga_fname))
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
#'
#' @keywords internal
#' @noRd
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
    req(length(keep) > 2)

    # determine if this is rna seq data
    rna_seq <- 'norm.factors' %in% colnames(Biobase::pData(eset))

    # subset eset and add group/pair
    eset <- eset[, keep]
    pdata$group <- pdata$`Group name`

    pair <- pdata$Pair
    if (any(!is.na(pair))) pdata$pair <- pair

    Biobase::pData(eset) <- pdata

    # filter rows by expression
    if (rna_seq) eset <- rkal::filter_genes(eset)

    # rlog normalize
    eset <- crossmeta::add_vsd(eset, rna_seq = rna_seq, vsd_path = vsd_path())
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
      req(identical(row.names(svobj$sv), colnames(eset)))
    }

    eset <- crossmeta::add_adjusted(eset, svobj, numsv, adj_path = adj_path())

    # use SYMBOL as annotation
    # keep unique symbol based on row IQRs
    eset <- crossmeta::iqr_replicates(eset, keep_path = keep_path())

    return(eset)
  })
  return(explore_eset)
}

