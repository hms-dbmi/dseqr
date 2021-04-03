#' Logic for Single Cell Tab
#'
#' @inheritParams bulkPage
#' @export
#'
#' @return Called with \link[shiny]{callModule} to generate logic for
#'   single-cell tab.
#'
scPage <- function(input, output, session, sc_dir, indices_dir, is_mobile) {

  # the analysis and options
  scForm <- callModule(scForm, 'form',
                       sc_dir = sc_dir,
                       indices_dir = indices_dir,
                       is_mobile = is_mobile)

  # cluster plot in top right
  callModule(scClusterPlot, 'cluster_plot',
             scseq = scForm$scseq,
             annot = scForm$annot,
             selected_cluster = scForm$selected_cluster,
             dataset_name = scForm$dataset_name,
             cluster_plot = scForm$cluster_plot,
             is_mobile = is_mobile)

  # cluster comparison plots ---

  callModule(scMarkerPlot, 'marker_plot_cluster',
             scseq = scForm$scseq,
             custom_metrics = scForm$custom_metrics,
             selected_feature = scForm$clusters_gene,
             dataset_name = scForm$dataset_name,
             plots_dir = scForm$plots_dir,
             feature_plot_clusters = scForm$feature_plot_clusters)

  callModule(scBioGpsPlot, 'biogps_plot',
             selected_gene = scForm$clusters_gene,
             species = scForm$species)


  callModule(scRidgePlot, 'ridge_plot',
             selected_gene = scForm$clusters_gene,
             selected_cluster = scForm$clusters_cluster,
             scseq = scForm$scseq,
             annot = scForm$annot,
             plots_dir =scForm$plots_dir)


  # sample comparison plots ---

  callModule(scSampleMarkerPlot, 'left',
             selected_gene = scForm$samples_gene,
             plot_fun = scForm$samples_pfun_left)


  callModule(scSampleMarkerPlot, 'right',
             selected_gene = scForm$samples_gene,
             plot_fun = scForm$samples_pfun_right)

  callModule(scSampleMarkerPlot, 'right_bottom',
             selected_gene = scForm$samples_gene,
             plot_fun = scForm$samples_pfun_right_bottom)




  observe({
    toggle(id = "sample_comparison_row",  condition = scForm$comparison_type() == 'samples')
    toggle(id = "cluster_comparison_row", condition = scForm$comparison_type() == 'clusters')
  })


  observe({
    toggle(id = 'biogps_container', condition = !scForm$show_ridge())
    toggle(id = 'ridge_container', condition = scForm$show_ridge())
  })

  return(NULL)
}


#' Logic for form on Single Cell Exploration page
#'
#' @keywords internal
#' @noRd
scForm <- function(input, output, session, sc_dir, indices_dir, is_mobile) {

  # updates if new integrated or subset dataset
  new_dataset <- reactiveVal()
  observe(new_dataset(scIntegration()))
  observe(new_dataset(scSubset()))

  # directory with cluster resolution independent stuff
  dataset_dir <- reactive({
    dataset <- scDataset$dataset_name()
    if (is.null(dataset)) return(NULL)
    file.path(sc_dir, dataset)
  })

  # directory with cluster resolution dependent stuff
  resoln_dir <- reactive({
    dataset_dir <- dataset_dir()
    resoln <- resoln()
    req(dataset_dir, resoln)
    resoln_dir <- file.path(dataset_dir, paste0('snn', resoln))
    if (!dir.exists(resoln_dir)) dir.create(resoln_dir)
    return(resoln_dir)
  })

  # directory for caching plots in
  plots_dir <- reactive({
    req(resoln_dir())
    dir <- file.path(resoln_dir(), 'plots')
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    return(dir)
  })

  annot <- reactiveVal()
  annot_path <- reactive(file.path(resoln_dir(), 'annot.qs'))
  observe(annot(qs::qread(annot_path())))

  observe(toggle('form_container', condition = scDataset$dataset_exists()))

  # update scseq with cluster changes (from resolution)
  scseq_clusts <- reactive({
    scseq <- scDataset$scseq()
    if (is.null(scseq)) return(NULL)

    clusters <- scResolution$clusters()
    resoln <-  resoln()
    plots_dir <- plots_dir()

    if (!is.null(clusters)) scseq$cluster <- clusters

    return(scseq)
  })

  # update scseq with annotation changes and custom metrics
  resoln <- reactive(scResolution$resoln())

  scseq <- reactive({
    scseq <- scseq_clusts()
    annot <- annot()
    metrics <- scClusterGene$saved_metrics()

    if (!isTruthy(annot) | !isTruthy(scseq)) return(NULL)
    if (!is.null(metrics)) try(scseq@colData <- cbind(scseq@colData, metrics), silent = TRUE)
    try(levels(scseq$cluster) <- annot, silent = TRUE)

    return(scseq)
  })

  qc_metrics <- reactive({
    scseq <- scseq()
    if(is.null(scseq)) return(NULL)
    metrics <- scseq@colData
    qc <- colnames(metrics)

    names(qc) <- sapply(metrics, class)

    qc <- qc[names(qc) %in% c('numeric', 'logical')]
    qc <- qc[!grepl('^sum$|^total$|^subsets|^percent|^sizeFactor', qc)]
  })

  # show the toggle if dataset is integrated
  observe({
    toggle(id = "comparison_toggle_container",  condition = scDataset$is_integrated())
  })



  # show appropriate inputs based on comparison type
  observe({
    toggle(id = "cluster_comparison_inputs",  condition = comparisonType() == 'clusters')
    toggle(id = "sample_comparison_inputs",  condition = comparisonType() == 'samples')
    toggle(id = "labels_comparison_inputs",  condition = comparisonType() == 'labels')
  })

  selected_cluster <- reactiveVal('')
  observe({
    type <- comparisonType()
    req(type)

    old <- isolate(selected_cluster())
    new <- switch(type,
                  'clusters' = scClusterComparison$selected_cluster(),
                  'samples' = scSampleComparison$selected_cluster(),
                  'labels' = scLabelsComparison$selected_cluster())

    if (is.null(new)) new <- ''
    if (new != old) selected_cluster(new)
  })


  # the dataset and options
  scDataset <- callModule(scSelectedDataset, 'dataset',
                          sc_dir = sc_dir,
                          new_dataset = new_dataset,
                          indices_dir = indices_dir)

  scClusterPlots <- clusterPlots(plots_dir, scseq = scseq_clusts)


  # label transfer between datasets
  # show/hide label transfer forms
  observe({
    toggle(id = "label-resolution-form", anim = TRUE, condition = scDataset$show_label_resoln())
  })
  scLabelTransfer <- callModule(labelTransferForm, 'transfer',
                                sc_dir = sc_dir,
                                dataset_dir = dataset_dir,
                                resoln_dir = resoln_dir,
                                resoln_name = scResolution$resoln_name,
                                annot_path = annot_path,
                                datasets = scDataset$datasets,
                                dataset_name = scDataset$dataset_name,
                                scseq = scDataset$scseq,
                                species = scDataset$species,
                                clusters = scResolution$clusters,
                                show_label_resoln = scDataset$show_label_resoln)

  # adjust resolution of dataset
  scResolution <- callModule(resolutionForm, 'resolution',
                             sc_dir = sc_dir,
                             resoln_dir = resoln_dir,
                             dataset_dir = dataset_dir,
                             dataset_name = scDataset$dataset_name,
                             scseq = scDataset$scseq,
                             snn_graph = scDataset$snn_graph,
                             annot_path = annot_path,
                             show_label_resoln = scDataset$show_label_resoln)

  # dataset integration
  scIntegration <- callModule(integrationForm, 'integration',
                              sc_dir = sc_dir,
                              datasets = scDataset$datasets,
                              selected_dataset = scDataset$dataset_name,
                              show_integration = scDataset$show_integration)

  # dataset subset
  scSubset <- callModule(subsetForm, 'subset',
                         sc_dir = sc_dir,
                         scseq = scseq,
                         datasets = scDataset$datasets,
                         selected_dataset = scDataset$dataset_name,
                         show_subset = scDataset$show_subset,
                         is_integrated = scDataset$is_integrated)


  # comparison type
  comparisonType <- callModule(comparisonType, 'comparison',
                               scseq = scseq,
                               is_integrated = scDataset$is_integrated)



  # the selected cluster/gene for cluster comparison
  scClusterComparison <- callModule(clusterComparison, 'cluster',
                                    sc_dir = sc_dir,
                                    dataset_dir = dataset_dir,
                                    dataset_name = scDataset$dataset_name,
                                    resoln_dir = resoln_dir,
                                    resoln = resoln,
                                    scseq = scseq,
                                    annot_path = annot_path,
                                    annot = annot,
                                    ref_preds = scLabelTransfer$pred_annot,
                                    clusters = scResolution$clusters)

  scClusterGene <- callModule(selectedGene, 'gene_clusters',
                              scseq = scseq_clusts,
                              dataset_name = scDataset$dataset_name,
                              resoln_name = scResolution$resoln_name,
                              resoln_dir = resoln_dir,
                              is_integrated = scDataset$is_integrated,
                              selected_markers = scClusterComparison$selected_markers,
                              selected_cluster = scClusterComparison$selected_cluster,
                              qc_metrics = qc_metrics,
                              type = 'clusters')



  # the selected clusters/gene for sample comparison
  scSampleComparison <- callModule(scSampleComparison, 'sample',
                                   input_scseq = scseq,
                                   input_annot = annot,
                                   dataset_dir = dataset_dir,
                                   resoln_dir = resoln_dir,
                                   plots_dir = plots_dir,
                                   feature_plot = scClusterPlots$feature_plot_samples,
                                   dataset_name = scDataset$dataset_name,
                                   sc_dir = sc_dir,
                                   is_integrated = scDataset$is_integrated,
                                   show_dprimes = scSampleGene$show_dprimes,
                                   comparison_type = comparisonType,
                                   exclude_ambient = scSampleGene$exclude_ambient,
                                   applied = scResolution$applied,
                                   is_mobile = is_mobile)

  scSampleGene <- callModule(selectedGene, 'gene_samples',
                             scseq = scDataset$scseq,
                             dataset_name = scDataset$dataset_name,
                             resoln_name = scResolution$resoln_name,
                             resoln_dir = resoln_dir,
                             is_integrated = scDataset$is_integrated,
                             selected_markers = scSampleComparison$top_table,
                             selected_cluster = scSampleComparison$selected_cluster,
                             type = 'samples',
                             ambient = scSampleComparison$ambient)



  scLabelsComparison <- callModule(scLabelsComparison, 'labels',
                                   cluster_choices = scSampleComparison$cluster_choices)



  return(list(
    scseq = scseq,
    samples_gene = scSampleGene$selected_gene,
    clusters_gene = scClusterGene$selected_gene,
    custom_metrics = scClusterGene$custom_metrics,
    show_ridge = scClusterGene$show_ridge,
    has_replicates = scSampleComparison$has_replicates,
    samples_pfun_left = scSampleComparison$pfun_left,
    samples_pfun_right = scSampleComparison$pfun_right,
    samples_pfun_right_bottom = scSampleComparison$pfun_right_bottom,
    clusters_cluster = scClusterComparison$selected_cluster,
    samples_cluster = scSampleComparison$selected_cluster,
    labels_cluster = scLabelsComparison$selected_cluster,
    selected_cluster = selected_cluster,
    comparison_type = comparisonType,
    dataset_name = scDataset$dataset_name,
    species = scDataset$species,
    plots_dir = plots_dir,
    feature_plot_clusters = scClusterPlots$feature_plot_clusters,
    cluster_plot = scClusterPlots$cluster_plot,
    annot = annot
  ))
}

clusterPlots <- function(plots_dir, scseq) {


  # create feature and cluster plot for future updating
  feature_plot_path <- reactive(file.path(plots_dir(), 'feature_plot.qs'))

  feature_plot_clusters <- reactive({
    plot_path <- feature_plot_path()
    plot <- qread.safe(plot_path)

    if (is.null(plot)) {
      scseq <- scseq()
      if (is.null(scseq)) return(NULL)
      plot <- plot_feature(scseq, row.names(scseq)[1])
      qs::qsave(plot, plot_path)
    }
    return(plot)
  })

  # this is ugly but avoids error from clusters/samples updating same plot
  feature_plot_samples <- reactive({
    plot_path <- feature_plot_path()
    plot <- qread.safe(plot_path)

    if (is.null(plot)) {
      scseq <- scseq()
      if (is.null(scseq)) return(NULL)
      plot <- plot_feature(scseq, row.names(scseq)[1])
      qs::qsave(plot, plot_path)
    }
    return(plot)
  })

  cluster_data_path <- reactive(file.path(plots_dir(), 'cluster_data.qs'))
  cluster_plot <- reactive({
    data_path <- cluster_data_path()
    plot_data <- qread.safe(data_path)

    if (!is.null(plot_data)) {
      plot <- plot_cluster(plot_data = plot_data)

    } else {
      scseq <- scseq()
      if (is.null(scseq)) return(NULL)
      plot <- plot_cluster(scseq)
      qs::qsave(plot$data, data_path)
    }

    return(plot)
  })

  return(list(
    cluster_plot = cluster_plot,
    feature_plot_samples = feature_plot_samples,
    feature_plot_clusters = feature_plot_clusters
  ))
}



#' Logic for selected dataset part of scForm
#'
#' @keywords internal
#' @noRd
scSelectedDataset <- function(input, output, session, sc_dir, new_dataset, indices_dir) {
  dataset_inputs <- c('selected_dataset', 'show_integration', 'show_label_resoln')
  options <- list(create = TRUE,
                  placeholder = 'Type name to add new single-cell dataset',
                  render = I('{option: scDatasetOptions, item: scDatasetItem}'),
                  searchField = c('optgroup', 'label'))

  # get directory with fastqs/h5 files
  roots <- c('single-cell' = sc_dir)
  shinyFiles::shinyDirChoose(input, "new_dataset_dir", roots = roots, restrictions = get_exclude_dirs(sc_dir))


  dataset_exists <- reactive(isTruthy(input$selected_dataset) & !is.create())

  dataset_name <- reactive({
    if (!dataset_exists()) return(NULL)
    ds <- datasets()
    ds <- ds$name[ds$value == input$selected_dataset]
  })

  dataset_dir <- reactive(file.path(sc_dir, dataset_name()))
  snn_path <- reactive(file.path(dataset_dir(), 'snn_graph.qs'))

  # load scseq
  scseq <- reactive({
    if (!isTruthy(dataset_name())) return(NULL)
    disableAll(dataset_inputs)
    scseq <- load_scseq(dataset_dir())
    enableAll(dataset_inputs)
    return(scseq)
  })

  # load snn graph
  snn_graph <- reactive({
    snn_path <- snn_path()

    if (file.exists(snn_path)) {
      snn_graph <- qs::qread(snn_path)

    } else {
      disableAll(dataset_inputs)
      scseq <- scseq()
      if (!isTruthy(scseq)) return(NULL)

      snn_graph <- get_snn_graph(scseq)
      qs::qsave(snn_graph, snn_path)
      enableAll(dataset_inputs)
    }

    return(snn_graph)
  })

  is_integrated <- reactive({
    dataset_name <- dataset_name()
    req(dataset_name)
    integrated <- qread.safe(file.path(sc_dir, 'integrated.qs'))
    return(dataset_name %in% integrated)
  })

  species <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)
    scseq@metadata$species
  })

  prev_datasets <- reactiveVal()
  datasets <- reactive({
    # reactive to new single cell datasets
    new_dataset()
    datasets <- get_sc_dataset_choices(sc_dir)
    prev <- isolate(prev_datasets())
    curr <- isolate(input$selected_dataset)

    datasets <- keep_curr_selected(datasets, prev, curr)
    prev_datasets(datasets)
    return(datasets)
  })



  # update previously selected dataset on-file if changes
  prev_path <- file.path(sc_dir, 'prev_dataset.qs')

  observe({
    sel <- dataset_name()
    req(sel)
    qs::qsave(sel, prev_path)
  })


  # are we creating a new dataset?
  is.create <- reactive({
    dataset_name <- input$selected_dataset
    datasets <- datasets()
    if (!isTruthy(dataset_name)) return(FALSE)

    !dataset_name %in% datasets$value
  })

  uploadModal <- function() {
    label <- "Click upload or drag files:"
    label_title <- "Accepts 10X *.fastq.gz or Cell Ranger files (*.h5 or matrix.mtx, barcodes.tsv, and genes.tsv)"
    label <- tags$span(label,
                       title = label_title,
                       span(class = "hover-info",
                            icon("info", "fa-fw")))

    modalDialog(
      fileInput(session$ns('up_raw'), label=label, buttonLabel = 'upload', accept = c('.h5', '.tsv', '.fastq.gz'), multiple = TRUE),
      title = 'Upload or Select Existing?',
      size = 's',
      footer = tagList(
        modalButton("Cancel"),
        actionButton(session$ns("click_existing"), "Select Existing", class = 'pull-left btn-warning')
      ),
      easyClose = FALSE,
    )
  }


  # open shinyFiles selector if creating
  observe({
    req(is.create())
    showModal(uploadModal())
  })


  # move uploaded to destination
  observeEvent(input$up_raw, {
    df <- input$up_raw
    sel <- input$selected_dataset
    req(df, sel)

    dataset_dir <- file.path(sc_dir, input$selected_dataset)
    dir.create(dataset_dir, showWarnings = FALSE)

    for (i in 1:nrow(df)) {
      dpath <- df$datapath[i]
      fpath <- file.path(dataset_dir, df$name[i])
      file.move(dpath, fpath)
    }

    removeModal()
    Sys.sleep(1)
    new_dataset_dir(dataset_dir)
  })

  observeEvent(input$click_existing, {
    removeModal()
    Sys.sleep(1)
    shinyjs::click('new_dataset_dir')
  })

  # get path to dir with new dataset files
  new_dataset_dir <- reactiveVal()
  observe({
    new_dataset_dir <- input$new_dataset_dir

    # need selected subfolder
    # will be integer on create
    req(!methods::is(new_dataset_dir, 'integer'))

    dir <- shinyFiles::parseDirPath(roots, new_dataset_dir)
    new_dataset_dir(as.character(dir))
  })

  observeEvent(input$selected_dataset, {
    new_dataset_dir(NULL)
  })

  # ask for confirmation after folder selection
  observeEvent(new_dataset_dir(), showModal(confirmModal(session)))


  metric_choices <- c('low_lib_size',
                      'low_n_features',
                      'high_subsets_mito_percent',
                      'low_subsets_ribo_percent',
                      'high_doublet_score')

  # run single-cell quantification
  quants <- reactiveValues()
  pquants <- reactiveValues()
  deselect_dataset <- reactiveVal(0)
  observeEvent(input$confirm_quant, {

    metrics <- input$qc_metrics
    # none, all, all and none: can't combine
    if (length(metrics) > 1 && !all(metrics %in% metric_choices)) return(NULL)

    if (!isTruthy(metrics)) metrics <- 'none'
    if (metrics[1] == 'all') metrics <- metric_choices

    removeModal()

    fastq_dir <- new_dataset_dir()
    dataset_name <- input$selected_dataset
    azimuth_ref <- input$azimuth_ref

    if (metrics[1] == 'all and none') {
      opts <- list(
        list(dataset_name = paste0(dataset_name, '_QC0'),
             metrics = NULL,
             founder = dataset_name),
        list(dataset_name = paste0(dataset_name, '_QC1'),
             metrics = metric_choices,
             founder = dataset_name))


    } else {
      if (metrics[1] == 'none') metrics <- NULL
      opts <- list(
        list(dataset_name = dataset_name,
             metrics = metrics,
             founder = NULL))
    }

    quants[[dataset_name]] <- callr::r_bg(
      func = run_load_raw_scseq,
      package = 'dseqr',
      args = list(
        opts = opts,
        fastq_dir = fastq_dir,
        sc_dir = sc_dir,
        indices_dir = indices_dir,
        azimuth_ref = azimuth_ref
      )
    )

    progress <- Progress$new(max=10*length(opts))
    msg <- paste(stringr::str_trunc(dataset_name, 33), "import:")
    progress$set(message = msg, value = 0)
    pquants[[dataset_name]] <- progress

    deselect_dataset(deselect_dataset()+1)
  })


  observe({
    invalidateLater(5000, session)
    handle_sc_progress(quants, pquants, new_dataset)

  })


  observe({
    datasets <- datasets()
    datasets <- datasets_to_list(datasets)
    updateSelectizeInput(session, 'selected_dataset', selected = isolate(input$selected_dataset), choices = datasets, options = options)
  })

  observeEvent(deselect_dataset(), {
    req(deselect_dataset())
    datasets <- datasets()
    datasets <- datasets_to_list(datasets)
    updateSelectizeInput(session, 'selected_dataset', choices = datasets, options = options)
  })

  # show/hide integration/label-transfer forms
  show_integration <- reactive(input$show_integration %% 3 == 2)
  show_subset <- reactive(input$show_integration %% 3 == 1)
  show_label_resoln <- reactive(input$show_label_resoln %% 2 == 1)

  observe({
    icon <- icon('object-ungroup', 'far fa-fw')
    if (show_integration()) icon <- icon('object-group', 'far fa-fw')

    updateActionButton(session, 'show_integration', icon = icon)
  })

  # hide integration/label-transfer buttons no dataset
  observe({
    toggle('show_label_resoln-parent', condition = dataset_exists())
  })


  return(list(
    dataset_name = dataset_name,
    scseq = scseq,
    snn_graph = snn_graph,
    datasets = datasets,
    show_integration = show_integration,
    show_subset = show_subset,
    show_label_resoln = show_label_resoln,
    is_integrated = is_integrated,
    dataset_exists = dataset_exists,
    species = species
  ))
}

# modal to confirm adding single-cell dataset
confirmModal <- function(session, type = c('quant', 'subset')) {
  qc <- NULL

  if (type[1] == 'quant') {
    label <- 'Quantify'
    id <- 'confirm_quant'
    qc <- selectizeInput(
      session$ns('qc_metrics'),
      HTML('Select <a href="https://docs.dseqr.com/docs/single-cell/quality-control/" target="_blank">QC</a> metrics:'),
      choices = c('all and none', 'all', 'none', metric_choices),
      selected = 'all and none',
      multiple = TRUE)

  } else {
    label <- paste(toupper(substr(type, 1, 1)), substr(type, 2, nchar(type)), sep="")
    id <- paste0('confirm_', type)
  }
  azi <- selectizeInput(
    session$ns('azimuth_ref'),
    HTML('Select <a href="https://azimuth.hubmapconsortium.org/" target="_blank">Azimuth</a> reference:'),
    choices = c('', 'human_pbmc'),
    options = list(placeholder = 'optional'))

  UI <- div(qc, azi)

  modalDialog(
    UI,
    title = 'Create new single-cell dataset?',
    size = 's',
    footer = tagList(
      modalButton("Cancel"),
      actionButton(session$ns(id), label, class = 'pull-left btn-warning')
    )
  )
}




#' Logic for selecting cluster to plot label origin for integrated dataset
#'
#' @keywords internal
#' @noRd
scLabelsComparison <- function(input, output, session, cluster_choices) {
  contrast_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))

  observe({
    updateSelectizeInput(session, 'selected_cluster',
                         choices = rbind(NA, cluster_choices()),
                         options = contrast_options, server = TRUE)
  })

  return(list(
    selected_cluster = reactive(input$selected_cluster)
  ))
}

#' Take currently selected row and keep it in the same position
#'
#' Used to stop dataset change when adding new datasets
#'
#' @param datasets current data.frame of single-cell datasets
#' @param prev previous data.frame of single-cell datasets
#' @param curr name of currently selected row from \code{prev}
#'
#' @return \code{datasets} with \code{curr} at same row as it is in \code{prev}
#'
keep_curr_selected <- function(datasets, prev, curr) {

  if (!isTruthy(prev) || !isTruthy(curr)) return(datasets)

  # get currently selected row
  curr <- as.numeric(curr)
  curr_text <- do.call(paste0, prev[curr, ])

  # position in new datasets
  new_posn <- which(curr_text == do.call(paste0, datasets))

  # move so that row at new_posn is at curr
  idx <- seq_len(nrow(datasets))
  idx_new <- replace(idx, c(curr, new_posn), c(new_posn, curr))
  datasets <- datasets[idx_new, ]
  datasets$value <- idx
  return(datasets)

}




#' Logic for single-cell sample comparison plots
#'
#' setup to allow for ggplot/plotly
#'
#'
#' @keywords internal
#' @noRd
scSampleMarkerPlot <- function(input, output, session, selected_gene, plot_fun, is_mobile) {

  res <- reactive({
    gene <- selected_gene()
    req(gene)
    suppressMessages(plot_fun()(gene))
  })

  height <- reactive(ifelse(is_plotly(), 1, res()$height))

  is_plotly <- reactive(methods::is(res(), 'plotly'))


  filename <- function() {
    fname <- plot()$labels$title
    fname <- gsub(':', '', fname)
    paste0(fname, '.csv')
  }

  plot <- reactive({
    if (!is_plotly()) return(res()$plot)
    return(NULL)
  })


  content <- function(file) {
    d <- plot()$data

    # clean up data for ridgeplots
    if (!'TSNE1' %in% colnames(d)) {
      d <- plot()$data[, c('x', 'y')]
      colnames(d) <- c(selected_gene(), 'sample')
    }

    utils::write.csv(d, file)
  }


  callModule(shinydlplot::downloadablePlot,
             "plot",
             plot = plot,
             filename = filename,
             content = content,
             height = height)


  output$plotly <- plotly::renderPlotly(if (is_plotly()) res() else NULL)
}

#' Logic for label transfer between datasets
#'
#' @keywords internal
#' @noRd
labelTransferForm <- function(input, output, session, sc_dir, dataset_dir, resoln_dir, resoln_name, annot_path, datasets, dataset_name, scseq, species, clusters, show_label_resoln) {
  label_transfer_inputs <- c('transfer_study', 'submit_transfer', 'overwrite_annot', 'ref_name', 'resoln', 'apply_update', 'reset_resoln')
  options <- list(render = I('{option: transferLabelOption, item: scDatasetItemDF}'))

  ref_preds <- reactiveVal()
  new_preds <- reactiveVal()
  new_annot <- reactiveVal()

  preds_path <- reactive(file.path(resoln_dir(), 'preds.qs'))


  # saved label transfer predictions
  preds <- reactive({
    new_preds()

    # load previously saved reference preds
    preds_path <- preds_path()
    preds <- if (file.exists(preds_path)) qs::qread(preds_path) else list()
    preds <- validate_preds(preds, sc_dir)
    return(preds)
  })

  observeEvent(resoln_name(), new_preds(NULL))

  # update annotation transfer choices
  observe({
    preds <- preds()

    datasets <- datasets()
    dataset_name <- dataset_name()
    species <- species()
    req(preds, datasets, species)

    choices <- get_label_transfer_choices(datasets, dataset_name, preds, species)
    updateSelectizeInput(session,
                         'ref_name',
                         choices = choices,
                         server = TRUE,
                         selected = isolate(new_preds()),
                         options = options)
  })



  query <- reactive({
    query_path <- scseq_part_path(sc_dir, resoln_name(), 'scseq_sample')
    if (!file.exists(query_path)) return(NULL)

    qs::qread(query_path)
  })

  # submit annotation transfer
  observeEvent(input$ref_name, {

    query_name <- dataset_name()
    ref_name <- input$ref_name
    preds <- preds()

    req(ref_name != 'reset')
    req(query_name, ref_name, preds)
    req(!ref_name %in% names(preds))
    req(show_label_resoln())

    query <- query()
    if (!isTruthy(query)) {
      showModal(warnApplyModal('transfer'))
      updateSelectizeInput(session, 'ref_name', selected = NULL)
      return(NULL)
    }

    disableAll(label_transfer_inputs)

    # Create a Progress object
    progress <- Progress$new()
    progress$set(message = "Transfering labels", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())

    # Create a callback function to update progress.
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 3
      }
      progress$set(value = value, detail = detail)
    }
    n = 3

    # get arguments for SingleR
    tab <- NULL
    senv <- loadNamespace('celldex')
    updateProgress(1/n)

    if (ref_name %in% ls(senv)) {
      ref <- get(ref_name, envir = senv)()
      labels <- ref$label.main
      genes <- 'de'


    } else {
      ref_subname <- get_resoln_name(sc_dir, ref_name)
      ref_path <- scseq_part_path(sc_dir, ref_subname, 'scseq_sample')
      ref_date <- file.info(ref_path)$ctime
      ref <- qs::qread(ref_path)

      # check if ref and query have the same founder
      rfound <- qs::qread(scseq_part_path(sc_dir, ref_name, 'founder'))
      qfound <- qs::qread(scseq_part_path(sc_dir, query_name, 'founder'))

      # use common cells to transfer labels if so
      cells <- intersect(colnames(ref), colnames(query))

      if (identical(qfound, rfound) && length(cells)) {

        ref_cluster <- ref[, cells]$cluster
        query_cluster <- query[, cells]$cluster
        tab <- table(assigned = ref_cluster, cluster = query_cluster)

      } else {
        markers_path <- scseq_part_path(sc_dir, ref_subname, 'top_markers')
        genes <- qs::qread(markers_path)

        # use aggregated reference for speed
        ref_path <- scseq_part_path(sc_dir, ref_subname, 'aggr_ref')
        if (file.exists(ref_path)) {
          ref <- qs::qread(ref_path)

        } else {
          set.seed(100)
          ref <- SingleR::aggregateReference(ref, labels=ref$cluster)
          qs::qsave(ref, ref_path)
        }

        # need until SingleR #77 fixed
        common <- intersect(row.names(ref), row.names(query))
        genes <- lapply(genes, function(x) {lapply(x, function(y) intersect(y, common))})

        labels <- ref$label
      }
    }

    updateProgress(2/n)
    if (is.null(tab)) {
      # take best label for each cluster
      pred <- SingleR::SingleR(test = query, ref = ref, labels = labels, genes = genes)
      tab <- table(assigned = pred$pruned.labels, cluster = query$cluster)
    }

    pred <- row.names(tab)[apply(tab, 2, which.max)]

    # keep track of date that reference was used so that can invalidate if overwritten
    if (exists('ref_date')) names(pred) <- ref_date

    preds_path <- preds_path()
    preds[[ref_name]] <- pred
    qs::qsave(preds, preds_path)

    new_preds(ref_name)
    ref_preds(preds[[ref_name]])

    enableAll(label_transfer_inputs)
  })


  # show transfered labels immediately upon selection if have
  observe({
    query_name <- resoln_name()
    ref_name <- input$ref_name
    req(query_name)

    # load previously saved reference preds
    preds_path <- scseq_part_path(sc_dir, query_name, 'preds')
    preds <- if (file.exists(preds_path)) qs::qread(preds_path) else list()

    # append reset labels
    clusters <- clusters()
    preds$reset <- levels(clusters)

    ref_preds(preds[[ref_name]])
  })



  pred_annot <- reactive({
    # react to new annotation
    new_annot()
    ref_name <- input$ref_name

    ref_preds <- ref_preds()
    query_resoln_name <- resoln_name()
    req(query_resoln_name)

    ref_resoln_name <- get_resoln_name(sc_dir, ref_name)

    # show saved annot if nothing selected or label transfer not open
    if (is.null(ref_preds)) {
      annot <- NULL

    } else if (!isTruthy(ref_name) | !show_label_resoln()) {
      annot_path <- annot_path()
      annot <- qs::qread(annot_path)

    } else {
      annot <- get_pred_annot(ref_preds, ref_resoln_name, query_resoln_name, sc_dir)
    }

    return(annot)
  }) %>% debounce(50)

  # overwrite annotation
  transferModal <- function() {
    modalDialog(
      tags$div('Saved annotation will be overwriten. This action cannot be undone.'),
      title = 'Are you sure?',
      size = 's',
      easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton(session$ns("confirm_overwrite"), "Overwrite", class = 'pull-left btn-warning')
      )
    )
  }

  # Show modal when button is clicked.
  observeEvent(input$overwrite_annot, {
    ref_name <- input$ref_name
    ref_preds <- ref_preds()
    resoln_name <- resoln_name()

    req(resoln_name)

    showModal(transferModal())
  })


  observeEvent(input$confirm_overwrite, {
    removeModal()
    ref_name <- input$ref_name
    ref_resoln_name <- get_resoln_name(sc_dir, ref_name)
    ref_preds <- ref_preds()
    query_resoln_name <- resoln_name()

    req(query_resoln_name)

    pred_annot <- get_pred_annot(ref_preds, ref_resoln_name, query_resoln_name, sc_dir)
    annot_path <- annot_path()
    qs::qsave(pred_annot, annot_path)

    new_annot(pred_annot)
  })




  return(list(
    pred_annot = pred_annot
  ))
}


#' Logic for leiden resolution slider
#'
#' @keywords internal
#' @noRd
resolutionForm <- function(input, output, session, sc_dir, resoln_dir, dataset_dir, dataset_name, scseq, snn_graph, annot_path, show_label_resoln) {
  resolution_inputs <- c('resoln', 'apply_update')

  prev_resoln <- reactiveVal()
  resoln_path <- reactiveVal()
  resoln <- reactiveVal()

  observeEvent(show_label_resoln(), {
    if (!show_label_resoln()) resoln(prev_resoln())
    else (resoln(input[[rname()]]))
  })

  # updateSliderInput removes focus preventing keyboard interaction
  observe({
    clusters <- clusters()
    nclus <- length(levels(clusters))
    type <- ifelse(is_azimuth(), 'nclus_azi', 'nclus')
    shinyjs::html(type, nclus)
  })

  observeEvent(input[[rname()]], {resoln(input[[rname()]])}, ignoreInit = TRUE)

  rname <- reactiveVal('resoln')
  is_azimuth <- reactiveVal(FALSE)

  observe({
    shinyjs::toggle('resoln_container', condition=!is_azimuth())
    shinyjs::toggle('resoln_azi_container', condition=is_azimuth())
  })

  observeEvent(dataset_dir(),  {
    dataset_dir <- dataset_dir()
    req(dataset_dir)

    # restrict resolutions if azimuth
    apath <- file.path(dataset_dir(), 'azimuth_ref.qs')
    axist <- file.exists(apath)
    is_azimuth(axist)

    rname(ifelse(axist, 'resoln_azi', 'resoln'))

    rpath <- file.path(dataset_dir(), 'resoln.qs')
    resoln_path(rpath)
    init <- qread.safe(rpath, 1)
    resoln(init)
    updateSliderInput(session, rname(), value=init)
  }, priority = 1)

  resoln_name <- reactive(file.path(dataset_name(), paste0('snn', resoln())))

  # clusters after change resolution
  clusters_path <- reactive(file.path(resoln_dir(), 'clusters.qs'))

  clusters <- reactive({
    resoln <- resoln()
    clusters_path <- clusters_path()
    clusters <- qread.safe(clusters_path)

    if (!is.null(clusters)) return(clusters)

    g <- snn_graph()
    if (is.null(g)) return(NULL)

    # stop resolution calc when change to dataset with different resolution
    prev_resoln <- prev_resoln()
    if (prev_resoln == resoln) return(NULL)

    clusters <- get_clusters(g, resolution = resoln)
    qs::qsave(clusters, clusters_path)

    # transfer annotation from prev clusters to new
    qs::qsave(levels(clusters), annot_path())
    transfer_prev_annot(resoln, prev_resoln, dataset_name(), sc_dir)

    return(clusters)
  })

  applied <- reactiveVal(FALSE)
  applied_path <- reactive(file.path(resoln_dir(), 'applied.qs'))
  observe(applied(file.exists(applied_path())))

  # update saved and prev resoln if applied
  observe({
    if (applied()) {
      resoln <- resoln()
      prev_resoln(resoln)
      qs::qsave(resoln, isolate(resoln_path()))
    }
  })

  same_resoln <- reactive(!is.null(prev_resoln()) && prev_resoln() == resoln())
  allow_update <- reactive(!applied() | !same_resoln())

  observe({
    toggleClass('apply_update', 'btn-warning', condition = allow_update())
    toggleState('apply_update', condition = allow_update())
    toggleState('reset_resoln', condition = allow_update())
  })

  observeEvent(input$reset_resoln, {
    updateSliderInput(session, rname(), value = prev_resoln())
  })


  observeEvent(input$apply_update, {
    resoln <- resoln()
    disableAll(resolution_inputs)
    progress <- Progress$new(session, min = 0, max = 4)
    progress$set(message = "Applying resolution update:", detail = 'loading', value = 1)
    on.exit(progress$close())

    dataset_name <- dataset_name()

    # need non-loom version
    scseq <- load_scseq_qs(dataset_dir())

    # add new clusters and run post clustering steps
    scseq$cluster <- clusters()
    run_post_cluster(scseq, dataset_name, sc_dir, resoln, progress, 1, reset_annot = FALSE)

    # mark as previously applied and re-enable
    applied(TRUE)
    qs::qsave(TRUE, applied_path())
    qs::qsave(resoln, resoln_path())
    enableAll(resolution_inputs)
  })


  return(list(
    clusters = clusters,
    resoln = resoln,
    resoln_name = resoln_name,
    applied = applied
  ))

}



#' Logic for subsetting a datatset
#'
#' @keywords internal
#' @noRd
subsetForm <- function(input, output, session, sc_dir, scseq, datasets, show_subset, selected_dataset, cluster_choices, is_integrated) {
  type <- name <- NULL
  contrastOptions <- list(render = I('{option: contrastOptions, item: contrastItem}'))

  subset_name <- reactive(input$subset_name)
  new_dataset <- reactiveVal()

  subset_inputs <- c('subset_name',
                     'submit_subset',
                     'subset_clusters',
                     'toggle_exclude',
                     'click_up')

  # show/hide integration forms
  observe({
    toggle(id = "subset-form", anim = TRUE, condition = show_subset())
  })

  is_include <- reactive({ input$toggle_exclude %% 2 != 0 })

  cluster_choices <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)
    annot <- levels(scseq$cluster)
    cluster_choices <- get_cluster_choices(annot, scseq = scseq)
    cluster_choices$value <- paste(selected_dataset(), cluster_choices$value, sep = '_')
    return(cluster_choices)
  })

  metric_choices <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)
    get_metric_choices(scseq)
  })

  exclude_choices <- reactive({
    choices <- cluster_choices()
    if (is.null(choices)) return(NULL)
    choices$type <- 'Cluster'

    metric_choices <- metric_choices()
    if (!is.null(metric_choices)) {
      metric_choices$type <- 'Metrics'
      choices <- rbind(metric_choices, choices)
      choices$pspace <- strrep('&nbsp;&nbsp;', max(0, 2 - nchar(choices$pcells)))
    }

    return(choices)
  })

  # change UI of exclude toggle
  observe({
    toggleClass(id = 'toggle_icon', 'fa-plus text-success', condition = is_include())
    toggleClass(id = 'toggle_icon', 'fa-minus text-warning', condition = !is_include())
  })


  # run integration
  subsets <- reactiveValues()
  psubsets <- reactiveValues()

  observeEvent(input$submit_subset, {
    browser()
    showModal(confirmModal(session, 'subset'))
  })

  observeEvent(input$confirm_subset, {

    subset_clusters <- input$subset_clusters
    subset_name <- input$subset_name
    azimuth_ref <- input$azimuth_ref
    cluster_choices <- cluster_choices()
    metric_choices <- metric_choices()
    is_include <- is_include()
    from_dataset <- selected_dataset()
    is_integrated <- is_integrated()
    hvgs <- hvgs()

    error_msg <- validate_subset(from_dataset, subset_name, subset_clusters, is_include, hvgs)

    if (is.null(error_msg)) {
      removeClass('name-container', 'has-error')

      exclude_clusters <- intersect(cluster_choices$value, subset_clusters)
      subset_metrics <- intersect(metric_choices$value, subset_clusters)

      if (is_include && length(exclude_clusters)) {
        exclude_clusters <- setdiff(cluster_choices$value, exclude_clusters)
      }

      founder <- get_founder(sc_dir, from_dataset)
      dataset_name <- subsets_name <- paste(founder, subset_name, sep = '_')

      # need exclude by cell name if integrated
      # because current clusters don't correspond to original clusters

      exclude_cells <- NULL
      if (is_integrated && length(exclude_clusters))
        exclude_cells <- get_exclude_cells(scseq(), exclude_clusters)

      subsets[[dataset_name]] <- callr::r_bg(
        func = subset_saved_scseq,
        package = 'dseqr',
        args = list(
          sc_dir = sc_dir,
          founder = founder,
          from_dataset = from_dataset,
          dataset_name = dataset_name,
          exclude_clusters = exclude_clusters,
          exclude_cells = exclude_cells,
          subset_metrics = subset_metrics,
          is_include = is_include,
          is_integrated = is_integrated,
          hvgs = hvgs,
          azimuth_ref = azimuth_ref
        )
      )

      progress <- Progress$new(max=ifelse(is_integrated, 9, 8))
      msg <- paste(stringr::str_trunc(dataset_name, 33), "subset:")
      progress$set(message = msg, value = 0)
      psubsets[[dataset_name]] <- progress


      # clear inputs
      updateTextInput(session, 'subset_name', value = '')

    } else {
      # show error message
      html('error_msg', html = error_msg)
      addClass('name-container', class = 'has-error')
    }
  })


  observe({
    invalidateLater(5000, session)
    handle_sc_progress(subsets, psubsets, new_dataset)
  })

  # upload custom HVGs
  hvgs <- reactiveVal()
  observeEvent(input$click_up, {
    hvgs(NULL)
    shinyjs::click('up_hvgs')
  })

  observe({
    infile <- input$up_hvgs
    req(infile)

    # make sure HVGs in scseq
    up_hvgs <- readLines(infile$datapath)
    genes <- isolate(row.names(scseq()))

    if (sum(up_hvgs %in% genes))
      hvgs(readLines(infile$datapath))
  })

  # make upload green when have data
  observe(toggleClass(id = "click_up", 'btn-success', condition = isTruthy(hvgs())))


  # update exclude clusters
  observe({

    updateSelectizeInput(session, 'subset_clusters',
                         choices = exclude_choices(),
                         selected = isolate(input$subset_clusters),
                         options = contrastOptions,
                         server = TRUE)
  })

  return(new_dataset)
}



#' Logic for integration form toggled by showIntegration
#'
#' @keywords internal
#' @noRd
integrationForm <- function(input, output, session, sc_dir, datasets, show_integration, selected_dataset) {
  type <- name <- NULL
  excludeOptions <- list(render = I('{option: excludeOptions, item: excludeOptions}'))

  integration_inputs <- c('ctrl_integration',
                          'integration_name',
                          'submit_integration',
                          'test_integration',
                          'integration_types',
                          'subset_clusters',
                          'click_up',
                          'click_dl',
                          'toggle_exclude')


  integration_name <- reactive(input$integration_name)

  # datasets() with server side selectize causes bug
  integration_choices <- reactive({
    ds <- datasets()
    if (!nrow(ds)) return(NULL)
    int  <- qread.safe(file.path(sc_dir, 'integrated.qs'))
    prev <- qread.safe(file.path(sc_dir, 'prev_dataset.qs'))
    ds <- ds[!ds$name %in% int & !ds$type %in% 'Previous Session', ]

    choices <- ds %>%
      dplyr::group_by(type) %>%
      dplyr::summarise(names = list(name))

    names(choices$names) <- choices$type
    choices$names
  })

  ctrl <- reactiveVal()
  test <- reactiveVal()
  new_dataset <- reactiveVal()
  selected_datasets <- reactive(c(test(), ctrl()))


  is_include <- reactive({ input$toggle_exclude %% 2 != 0 })
  have_contrast <- reactive(length(test()) && length(ctrl()))
  allow_pairs <- reactive(length(selected_datasets()) > 2 && have_contrast())
  allow_integration <- reactive(have_contrast() && isTruthy(input$integration_name) && isTruthy(input$integration_types))

  # show/hide integration forms
  observe({
    toggle(id = "integration-form", anim = TRUE, condition = show_integration())
  })

  observe(ctrl(input$ctrl_integration))
  observe(test(input$test_integration))

  # update test dataset choices
  observe({
    choices <- integration_choices()
    for (i in seq_along(choices)) {
      group <- choices[[i]]
      group <- group[!group %in% ctrl()]
      if (length(group)==1) group <- list(group)
      choices[[i]] <- group
    }

    updateSelectizeInput(session, 'test_integration', choices = choices, selected = isolate(test()))
  })

  # update control dataset choices
  observe({
    choices <- integration_choices()
    for (i in seq_along(choices)) {
      group <- choices[[i]]
      group <- group[!group %in% test()]
      if (length(group)==1) group <- list(group)
      choices[[i]] <- group
    }

    updateSelectizeInput(session, 'ctrl_integration', choices = choices, selected = isolate(ctrl()))
  })



  # hide pairing/subset if not enough datasets
  observe({
    toggleState(id = "click_dl", condition = allow_pairs())
    toggleState(id = "click_up", condition = allow_pairs())
    toggleState(id = 'submit_integration', condition = allow_integration())
  })


  # show cluster type choices if enough datasets
  observe(toggle(id = 'integration_types', condition = have_contrast()))


  # show exclude/exclude and new dataset only if something selected
  observe(toggle(id = 'exclude-container', condition = have_contrast()))
  observe(toggle(id = 'name-container', condition = have_contrast()))


  exclude_choices <- reactive({
    selected <- selected_datasets()
    colors <- get_palette(selected)
    choices <- get_exclude_choices(selected, sc_dir, colors)
    return(choices)
  })

  # update subset clusters
  observe({
    updateSelectizeInput(session, 'subset_clusters', choices = exclude_choices(),selected = isolate(input$subset_clusters), options = excludeOptions, server = TRUE)
  })

  # upload/download pairs
  pairs <- reactiveVal()

  observeEvent(input$click_up, {
    pairs(NULL)
    shinyjs::click('up_pairs')
  })

  observeEvent(input$click_dl, {
    shinyjs::click('dl_samples')
  })

  samples <- reactive({
    dataset_names <- c(test(), ctrl())
    req(dataset_names)
    data.frame(sample = dataset_names, pair = NA)
  })


  output$dl_samples <- downloadHandler(
    filename = 'samples.csv',
    content = function(con) {
      utils::write.csv(samples(), con, row.names = FALSE)
    }
  )

  observe({
    infile <- input$up_pairs
    req(infile)
    pairs(utils::read.csv(infile$datapath, row.names = 'sample'))
  })

  # make upload green when have data
  observe(toggleClass(id = "click_up", 'btn-success', condition = isTruthy(pairs())))

  # change UI of exclude toggle
  observe({
    toggleClass(id = 'toggle_icon', 'fa-plus text-success', condition = is_include())
    toggleClass(id = 'toggle_icon', 'fa-minus text-warning', condition = !is_include())
  })

  # run integration
  pintegs <- reactiveValues()
  integs <- reactiveValues()

  use_azimuth <- reactive('Azimuth' %in% input$integration_types)

  observe({
    toggle('azimuth_ref_container', condition = use_azimuth())
  })

  observeEvent(input$submit_integration, {

    subset <- exclude_clusters <- input$subset_clusters

    if (is_include() && length(subset)) {
      choices <- exclude_choices()
      exclude_clusters <- setdiff(choices$value, subset)
    }

    test <- test()
    ctrl <- ctrl()
    pairs <- pairs()

    integration_types <- input$integration_types
    integration_name <- input$integration_name
    azimuth_ref <- input$azimuth_ref

    if (!use_azimuth()) azimuth_ref <- NULL
    else if (is.null(azimuth_ref)) error_msg <- 'Select Azimuth reference'

    error_msg <- validate_integration(test, ctrl, pairs)

    if (is.null(error_msg)) {
      # clear error and disable button
      removeClass('name-container', 'has-error')

      integs[[integration_name]] <- callr::r_bg(
        func = run_integrate_saved_scseqs,
        package = 'dseqr',
        args = list(
          sc_dir = sc_dir,
          test = test,
          ctrl = ctrl,
          integration_name = integration_name,
          integration_types = integration_types,
          exclude_clusters = exclude_clusters,
          pairs = pairs,
          azimuth_ref = azimuth_ref
        )
      )

      progress <- Progress$new(max=8*length(integration_types))
      msg <- paste(stringr::str_trunc(integration_name, 33), "integration:")
      progress$set(message = msg, value = 0)
      pintegs[[integration_name]] <- progress

    } else {
      # show error message
      pairs(NULL)
      html('error_msg', html = error_msg)
      addClass('name-container', class = 'has-error')
    }

  })

  # progress monitoring of integration
  observe({
    invalidateLater(5000, session)
    handle_sc_progress(integs, pintegs, new_dataset)
  })

  return(new_dataset)
}



#' Logic for comparison type toggle for integrated datasets
#'
#' @keywords internal
#' @noRd
comparisonType <- function(input, output, session, scseq, is_integrated) {

  # always show clusters if not integrated
  observe({
    if(!is_integrated())
      shinyWidgets::updateRadioGroupButtons(session, 'comparison_type', selected = 'clusters')
  })

  return(reactive(input$comparison_type))
}

#' Logic for cluster comparison input
#'
#' @keywords internal
#' @noRd
clusterComparison <- function(input, output, session, sc_dir, dataset_dir, dataset_name, resoln_dir, resoln, scseq, annot_path, annot, ref_preds, clusters) {
  cluster_inputs <- c('selected_cluster', 'rename_cluster', 'show_contrasts', 'show_rename')

  contrast_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))
  selected_cluster <- reactiveVal()
  markers <- reactiveVal(list())

  # things that return for plotting
  selected_markers <- reactiveVal(NULL)

  show_contrasts <- reactive({ input$show_contrasts %% 2 != 0 })
  show_rename <- reactive((input$rename_cluster + input$show_rename) %% 2 != 0)

  test_cluster <- reactive({
    test_cluster <- input$selected_cluster
    req(test_cluster)
    gsub('-vs-.+?$', '', test_cluster)
  })

  # update data.frame for cluster/contrast choices
  choices <- reactive({
    clusters <- annot()
    req(clusters)
    req(scseq())

    if (show_contrasts()) {
      test <- isolate(test_cluster())
      choices <- get_contrast_choices(clusters, test)

    } else {
      choices <- get_cluster_choices(clusters, with_all=TRUE, scseq = scseq())
    }

    return(choices)
  })


  observe({
    ref_preds <- ref_preds()
    annot_path <- annot_path()
    if (!isTruthy(annot_path)) annot(NULL)
    else if (!is.null(ref_preds)) annot(ref_preds)
    else annot(qs::qread(annot_path))
  })

  sel_d <- reactive(input$selected_cluster) %>% debounce(20)

  observe({
    sel <- sel_d()
    prev <- isolate(selected_cluster())

    no.prev <- is.null(prev)
    is.new <- !is.null(sel) && sel != prev
    is.flip <- !show_contrasts() & grepl('-vs-', sel)

    if ((no.prev || is.new) & !is.flip) {
      selected_cluster(sel)
    }
  })

  observeEvent(input$show_contrasts, {
    if (!show_contrasts()) {
      test <- test_cluster()
      prev <- selected_cluster()

      if (prev != test) {

        selected_cluster(test)
      }
    }
  })


  # reset if switch dataset or resolution
  observeEvent(resoln_dir(), {
    markers(list())
    selected_cluster(NULL)
    annot(NULL)
    selected_markers(NULL)
  }, priority = 1)


  # show/hide rename and select panel
  observe({
    toggle(id = "rename_panel", condition = show_rename())
    toggle(id = "select_panel", condition = !show_rename())
  })


  # show/hide contrasts
  observe({
    # update icon on toggle
    icon <- ifelse(show_contrasts(), 'chevron-down', 'chevron-right')

    updateActionButton(session, 'show_contrasts', icon = shiny::icon(icon, 'fa-fw'))
    toggleClass(id = "show_contrasts", 'btn-primary', condition = show_contrasts())
  })


  # modify/save annot if rename a cluster
  observeEvent(input$rename_cluster, {
    req(input$new_cluster_name)

    # update reactive annotation
    choices <- choices()
    sel_clust <- selected_cluster()
    sel_idx <- gsub('-vs-\\d+$', '', sel_clust)
    sel_idx <- as.numeric(sel_idx)

    # use currently save annot as reference
    ref_preds <- ref_preds()
    mod_annot <- qs::qread(annot_path())
    mod_annot[sel_idx] <- ref_preds[sel_idx] <- input$new_cluster_name
    mod_annot <- make.unique(mod_annot, '_')

    # save on disc
    qs::qsave(mod_annot, annot_path())

    # update annot and set selected cluster to new name
    annot(mod_annot)
  })


  # update UI for contrast/cluster choices
  observeEvent(choices(), {
    choices <- choices()
    selected <- NULL

    if (!show_contrasts()) {
      choices <- rbind(NA, choices)
      selected <- isolate(selected_cluster())
    }

    updateSelectizeInput(session, 'selected_cluster',
                         choices = choices,
                         selected = selected,
                         options = contrast_options, server = TRUE)
  })


  # update ui for renaming a cluster
  observe({
    choices <- choices()
    name <- choices[choices$value == input$selected_cluster, 'name']

    if (!show_rename())
      updateTextInput(session, 'new_cluster_name', value = '', placeholder = paste('Type new name for', name, '...'))
  })


  # get/load markers if don't have
  observeEvent(input$selected_cluster, {
    sel <- input$selected_cluster
    resoln_dir <- resoln_dir()
    markers <- markers()
    req(sel, resoln_dir)

    if (sel %in% names(markers)) return(NULL)


    if (show_contrasts()) {
      con_markers <- get_contrast_markers(sel, resoln_dir)
      markers <- c(markers, con_markers)

    } else {
      markers_path <- file.path(resoln_dir, paste0('markers_', sel, '.qs'))
      if (!file.exists(markers_path)) {
        showModal(warnApplyModal('markers'))
        return(NULL)
      }
      markers[[sel]] <- qs::qread(markers_path)
    }

    markers(markers)
  })


  observe({
    sel <- debounce(selected_cluster, 20)
    sel <- sel()

    if (!is.null(sel)) {
      new <- markers()[[sel]]

      prev <- isolate(selected_markers())
      if (is.null(prev) || !identical(row.names(new), row.names(prev)))
        selected_markers(new)
    }
  })


  return(list(
    annot = annot,
    selected_markers = selected_markers,
    selected_cluster = selected_cluster
  ))
}


warnApplyModal <- function(type = c('markers', 'transfer', 'select')) {
  text <- switch(
    type[1],
    markers = 'Either apply or reset cluster resolution to sort marker genes by cluster.',
    transfer = 'Either apply or reset cluster resolution to transfer labels.'
  )

  modalDialog(
    text,
    title = 'Apply or Reset',
    size = 's',
    easyClose = TRUE
  )
}




#' Logic for selected gene to show plots for
#'
#' @keywords internal
#' @noRd
selectedGene <- function(input, output, session, dataset_name, resoln_name, resoln_dir, scseq, is_integrated, selected_markers, selected_cluster, type, qc_metrics = function()NULL, ambient = function()NULL) {
  gene_options <- list(render = I('{option: geneChoice, item: geneChoice}'))

  selected_gene <- reactiveVal(NULL)

  # toggle for excluding ambient
  exclude_ambient <- reactive({
    if (is.null(input$exclude_ambient)) return(TRUE)
    input$exclude_ambient %% 2 != 1
  })

  observe({
    toggleClass('exclude_ambient', class = 'btn-primary', condition = exclude_ambient())
  })


  # toggle for ridgeline
  gene_selected <- reactive({
    sel <- input$selected_gene
    scseq <- scseq()
    if (!isTruthyAll(sel, scseq)) return(FALSE)
    sel %in% row.names(scseq)
  })

  show_ridge <- reactive(input$show_ridge %% 2 != 1 | !gene_selected())
  observe(toggleClass(id = "show_ridge", 'btn-primary', condition = !show_ridge()))

  # toggle for showing custom metric
  show_custom_metric <- reactive(type != 'samples' && (input$show_custom_metric %%2 != 0))

  observe({
    toggle('custom_metric_panel', anim = TRUE, condition = show_custom_metric())
    if (show_custom_metric() & have_metric()) selected_gene(input$custom_metric)
  })

  observe(if (!show_custom_metric()) selected_gene(isolate(input$selected_gene)))
  observe(toggleClass('show_custom_metric', class = 'btn-primary', condition = show_custom_metric()))

  # toggle for showing dprimes plot
  show_dprimes <- reactive(type == 'samples' && (input$show_dprimes %%2 != 0))
  observe(toggleClass(id = "show_dprimes", 'btn-primary', condition = show_dprimes()))


  # hide buttons when not valid
  observe({
    toggle('show_dprimes-parent', condition = gene_selected() & is_integrated())
    toggle('genecards-parent', condition = gene_selected())
    toggle('show_ridge-parent', condition = gene_selected())
  })

  saved_metrics <- reactiveVal()
  custom_metrics <- reactiveVal()
  exist_metric_names <- reactive(c(row.names(scseq()),
                                   colnames(scseq()@colData),
                                   colnames(saved_metrics())
  ))


  have_metric <- reactive(input$custom_metric %in% colnames(custom_metrics()))
  allow_save <- reactive(!input$custom_metric %in% row.names(scseq()))

  observe(toggleState('save_custom_metric', condition = allow_save()))

  # for updating plot of current custom metric
  observe({
    metric <- input$custom_metric
    metric <- gsub('^ | $', '', metric)
    scseq <- scseq()
    req(metric, scseq, show_custom_metric())

    if (metric %in% exist_metric_names()) {
      selected_gene(metric)
      return(NULL)
    }

    res <- suppressWarnings(validate_metric(metric, scseq))
    res.na <- all(is.na(res[[1]]))
    req(!res.na)

    if (methods::is(res, 'data.frame')) {

      prev <- custom_metrics()
      if (!is.null(prev) && nrow(prev) != nrow(res)) {
        custom_metrics(NULL)
        return(NULL)
      }

      if (!is.null(prev)) res <- cbind(prev, res)
      res <- res[, unique(colnames(res)), drop = FALSE]
      custom_metrics(res)
      selected_gene(metric)
    }
  })

  # for saving custom metric for future sessions
  metrics_path <- reactive(file.path(resoln_dir(), 'saved_metrics.qs'))
  observe(saved_metrics(qread.safe(metrics_path())))

  observeEvent(input$save_custom_metric, {
    metric <- input$custom_metric
    prev <- saved_metrics()
    custom_metrics <- custom_metrics()

    # can remove custom metric by selecting as feature and saving empty
    if (!isTruthy(metric)) {
      feature <- selected_gene()
      req(feature %in% colnames(prev))
      res <- prev[, !colnames(prev) %in% feature, drop = FALSE]
      if (ncol(res) == 0) res <- NULL

    } else {
      res <- custom_metrics[, metric, drop = FALSE]
      if (!is.null(prev)) res <- cbind(prev, res)
    }

    qs::qsave(res, metrics_path())
    saved_metrics(res)
    updateTextInput(session, 'custom_metric', value = '')
  })


  filtered_markers <- reactive({

    markers <- selected_markers()
    if (is.null(markers)) return(NULL)

    ambient <- NULL
    if (exclude_ambient()) ambient <- ambient()

    supress.genes(markers, ambient)
  })

  scseq_genes <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)
    data.frame(1:nrow(scseq), row.names = row.names(scseq))
  })
  species <- reactive(scseq()@metadata$species)

  tx2gene <- reactive(dseqr.data::load_tx2gene(species()))
  # update marker genes based on cluster selection
  gene_choices <- reactive({

    qc_first <- FALSE
    markers <- filtered_markers()
    selected_cluster <- selected_cluster()
    qc_metrics <- qc_metrics()

    # will error if labels
    # also prevents intermediate redraws
    if (is.null(markers) & isTruthy(selected_cluster)) return(NULL)

    if (is.null(markers)) {
      qc_first <- TRUE
      markers <- scseq_genes()
    }

    if (is.null(markers) || is.null(scseq())) return(NULL)
    get_gene_choices(markers,
                     qc_metrics = qc_metrics,
                     qc_first = qc_first,
                     species = species(),
                     tx2gene = tx2gene())

  }) %>% debounce(50)

  # click genecards
  observeEvent(input$genecards, {
    gene_link <- paste0('https://www.genecards.org/cgi-bin/carddisp.pl?gene=', input$selected_gene)
    runjs(paste0("window.open('", gene_link, "')"))
  })


  # reset selected gene if dataset or cluster changes
  observe({
    sel <- input$selected_gene
    if (!isTruthy(sel)| !isTruthy(dataset_name())) return(NULL)
    else selected_gene(sel)
  })

  # reset custom metric if dataset changes
  observeEvent(resoln_name(), custom_metrics(NULL))


  # update choices
  observeEvent(gene_choices(), {
    choices <- gene_choices()
    if(is.null(choices)) return(NULL)

    updateSelectizeInput(session, 'selected_gene',
                         choices = rbind(NA, choices, fill = TRUE),
                         server = TRUE,
                         options = gene_options)
  })


  return(list(
    selected_gene = selected_gene,
    exclude_ambient = exclude_ambient,
    show_ridge = show_ridge,
    show_dprimes = show_dprimes,
    custom_metrics = custom_metrics,
    saved_metrics = saved_metrics
  ))


}


#' Logic for cluster plots
#'
#' @keywords internal
#' @noRd
scClusterPlot <- function(input, output, session, scseq, annot, selected_cluster, dataset_name, cluster_plot, is_mobile) {

  filename <- function() {
    paste0(dataset_name(), '_cluster_plot_data_', Sys.Date(), '.csv')
  }


  plot <- reactive({
    plot <- cluster_plot()
    annot <- annot()
    if (is.null(annot)) return(NULL)

    hl <- NULL
    cluster <- selected_cluster()
    cluster <- strsplit(cluster, '-vs-')[[1]]
    nclus <- length(annot)

    if (is_mobile() || length(annot) > 30) {
      annot <- as.character(seq_along(annot))
    }


    if (isTruthy(cluster) && nclus >= as.numeric(cluster))
      hl <- as.numeric(cluster)

    update_cluster_plot(plot, annot, hl)

  })


  content <- function(file) {
    data <- plot()$data
    utils::write.csv(data, file)
  }


  callModule(shinydlplot::downloadablePlot,
             "cluster_plot",
             plot = plot,
             filename = filename,
             content = content)


  return(list(
    plot = plot
  ))
}


#' Logic for marker feature plots
#'
#' @keywords internal
#' @noRd
scMarkerPlot <- function(input, output, session, scseq, selected_feature, dataset_name, plots_dir, feature_plot_clusters, custom_metrics = function()NULL) {


  plot <- reactive({
    feature <- selected_feature()
    scseq <- scseq()
    if (!isTruthy(feature) || !isTruthy(scseq)) return(NULL)

    metrics <- custom_metrics()
    cdata <- scseq@colData

    if (!is.null(metrics) && nrow(cdata) != nrow(metrics)) return(NULL)
    if (!is.null(metrics)) {
      cdata <- cbind(cdata, metrics)
      scseq@colData <- cdata
    }

    is_gene <- feature %in% row.names(scseq)
    is_feature <- feature %in% colnames(cdata)
    if (!is_gene && !is_feature) return(NULL)

    is_log <- is.logical(cdata[[feature]])
    is_num <- is_gene || is.numeric(cdata[[feature]])


    if (is_num) {
      plots_dir <- plots_dir()
      plot <- feature_plot_clusters()
      fdata <- get_feature_data(plots_dir, scseq, feature)
      pl <- update_feature_plot(plot, fdata, feature)

    } else if (is_log) {
      ft <- cdata[[feature]]
      scseq$cluster <- factor(ft, levels = c(FALSE, TRUE))
      ncells <- sum(ft)
      pcells <- round(ncells / length(ft) * 100)

      title <- paste0(feature, ' (', format(ncells, big.mark=","), ' :: ', pcells, '%)')
      pl <- plot_cluster(scseq, label = FALSE, label.index = FALSE, order = TRUE, title = title, cols =  c('lightgray', 'blue'))
    } else {
      pl <- NULL
    }

    return(pl)
  })


  filename <- function() {
    fname <- paste0(dataset_name(),
                    '_', selected_feature(),
                    '_marker_plot_data_', Sys.Date(), '.csv')
    return(fname)
  }

  content <- function(file) {
    plot <- plot()
    utils::write.csv(plot$data, file)
  }

  callModule(shinydlplot::downloadablePlot,
             "marker_plot",
             plot = plot,
             filename = filename,
             content = content)




  return(list(
    plot = plot
  ))
}

#' Logic for BioGPS plot
#'
#' @keywords internal
#' @noRd
scBioGpsPlot <- function(input, output, session, selected_gene, species) {
  SYMBOL <- NULL

  # plot BioGPS data
  output$biogps_plot <- renderPlot({
    species <- species()
    gene <- selected_gene()
    if (!length(gene)) return(NULL)
    if (species == 'Mus musculus') gene <- toupper(gene)
    if (!gene %in% biogps[, SYMBOL]) return(NULL)

    plot_biogps(gene)
  })
}

#' Logic for Ridge plot for clusters
#'
#' @keywords internal
#' @noRd
scRidgePlot <- function(input, output, session, selected_gene, selected_cluster, scseq, annot, plots_dir) {

  height <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) height <- 453
    else height <- length(levels(scseq$cluster))*38
    return(height)
  })


  gene_d <- selected_gene %>% debounce(20)
  clus_d <- selected_cluster %>% debounce(20)

  ridge_data <- reactive({
    gene <- gene_d()
    cluster <- clus_d()
    if (!isTruthy(gene)) return(NULL)
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)
    is.gene <- gene %in% row.names(scseq)
    is.num <- is.gene || is.numeric(scseq@colData[[gene]])
    if (!is.num) return(NULL)

    dat_path <- file.path(plots_dir(), paste(gene, cluster, 'cluster_ridgedat.qs', sep='-'))

    if (file.exists(dat_path)) {
      annot <- annot()
      if (is.null(annot)) return(NULL)
      rdat <- qs::qread(dat_path)
      rdat$clus_levs <- annot()

    } else {
      scseq <- scseq()
      if (is.null(scseq)) return(NULL)
      rdat <- get_ridge_data(gene, scseq, cluster, with_all = TRUE)
      qs::qsave(rdat, dat_path)
    }

    return(rdat)
  }) %>% debounce(20)

  plot <- reactive({
    ridge_data <- ridge_data()
    if (is.null(ridge_data)) return(NULL)
    VlnPlot(ridge_data = ridge_data)
  }) %>% debounce(20)

  content <- function(file){
    d <- ridge_data()$df[, c('x', 'y')]
    colnames(d) <- c(selected_gene(), 'cluster')
    utils::write.csv(d, file, row.names = FALSE)
  }

  filename <- reactive(paste0(selected_gene(), '.csv'))

  callModule(shinydlplot::downloadablePlot,
             "ridge_plot",
             plot = plot,
             filename = filename,
             content = content,
             height = height)
}


#' Logic for single cell cluster analyses for Single Cell, Drugs, and Pathways tabs
#'
#' @keywords internal
#' @noRd
scSampleComparison <- function(input, output, session, dataset_dir, resoln_dir, plots_dir, feature_plot, dataset_name, sc_dir, input_annot = function()NULL, input_scseq = function()NULL, show_dprimes = function()TRUE, is_integrated = function()TRUE, is_sc = function()TRUE, exclude_ambient = function()FALSE, comparison_type = function()'samples', applied = function()TRUE, is_mobile = function()FALSE) {
  contrast_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))
  input_ids <- c('click_dl', 'selected_cluster')

  # need for drugs tab
  scseq <- reactive({
    scseq <- input_scseq()
    if (!is.null(scseq)) return(scseq)

    dataset_dir <- dataset_dir()
    if (!isTruthy(dataset_dir)) return(NULL)

    scseq <- load_scseq(dataset_dir, default_clusters = FALSE)
    attach_clusters(scseq, resoln_dir())
  })

  annot <- reactive({
    annot <- input_annot()
    if (!is.null(annot)) return(annot)

    annot_path <- file.path(resoln_dir(), 'annot.qs')
    if (!file.exists(annot_path)) return(NULL)
    qs::qread(annot_path)
  })


  has_replicates <- reactive(qread.safe(file.path(resoln_dir(), 'has_replicates.qs')))
  species <- reactive(qread.safe(file.path(dataset_dir(), 'species.qs')))


  # update cluster choices in UI
  cluster_choices <- reactive({
    applied()

    integrated <- is_integrated()
    resoln_dir <- resoln_dir()

    if (!isTruthyAll(resoln_dir, integrated)) return(NULL)
    # need to be in sync (don't take from elsewhere)
    annot_path  <- file.path(resoln_dir, 'annot.qs')
    annot <- qread.safe(annot_path)

    if (is.null(annot)) return(NULL)

    # resolution must be applied
    applied <- file.exists(file.path(resoln_dir, 'applied.qs'))
    if (!applied) return(NULL)

    tryCatch({
      get_cluster_choices(clusters = c(annot, 'All Clusters'),
                          sample_comparison = TRUE,
                          resoln_dir = resoln_dir,
                          use_disk = TRUE,
                          top_tables = top_tables(),
                          has_replicates = has_replicates())
    },
    error = function(e) return(NULL))
  }) %>% debounce(20)

  observe({
    updateSelectizeInput(session, 'selected_cluster',
                         choices = rbind(NA, cluster_choices()),
                         options = contrast_options, server = TRUE)
  })


  # path to lmfit, cluster markers, drug query results, and goanna pathway results
  clusters_str <- reactive(collapse_sorted(input$selected_cluster))
  drug_paths <- reactive(get_drug_paths(resoln_dir(), clusters_str()))
  go_path <- reactive(file.path(resoln_dir(), paste0('go_', clusters_str(), '.qs')))
  kegg_path <- reactive(file.path(resoln_dir(), paste0('kegg_', clusters_str(), '.qs')))
  goana_path <- reactive(file.path(resoln_dir(), paste0('goana_', clusters_str(), '.qs')))
  kegga_path <- reactive(file.path(resoln_dir(), paste0('kegga_', clusters_str(), '.qs')))
  top_tables_paths <- reactive(file.path(resoln_dir(), 'top_tables.qs'))


  # require cluster markers and fit result
  lm_fit <- reactive(qread.safe(file.path(resoln_dir(), 'lm_fit_0svs.qs')))
  cluster_markers <- reactive(qread.safe(file.path(resoln_dir(), paste0('markers_', clusters_str(), '.qs'))))

  # plot functions for left
  sel <- reactive(input$selected_cluster)

  scseq_samp <- reactive(downsample_group(scseq()))


  pfun_left <- reactive({
    req(is_integrated())

    function(gene) {
      if(!isTruthy(gene)) return(NULL)

      if (show_dprimes() & is_integrated()) {
        if(is.null(top_tables())) return(NULL)
        annot <- annot()
        sel <- sel()
        amb <- exclude_ambient()
        annot_hash <- digest::digest(c(annot, amb), 'crc32')
        pname <- paste(gene, sel, annot_hash, 'dprimes.qs', sep='-')
        plot_path <- file.path(plots_dir(), pname)

        if (file.exists(plot_path)) {
          pfun <- qs::qread(plot_path)

        } else {
          pfun <- plot_scseq_dprimes(gene, annot, sel, top_tables(), amb)
          qs::qsave(pfun, plot_path)
        }

      } else {
        scseq <- scseq()
        if(is.null(scseq)) return(NULL)
        plots_dir <- plots_dir()
        gene_data <- get_feature_data(plots_dir, scseq, gene)

        # update base feature plot
        plot <- feature_plot()
        plot <- update_feature_plot(plot, gene_data, gene)
        plot <- plot_feature_sample(gene, scseq_samp(), 'test', plot=plot)
        pfun <- list(plot = plot, height = 453)
      }
      return(pfun)
    }
  })



  # plot functions for right
  pfun_right <- reactive({
    req(is_integrated())

    function(gene) {

      scseq <- scseq()
      if (!isTruthyAll(scseq, gene)) return(NULL)

      plots_dir <- plots_dir()
      gene_data <- get_feature_data(plots_dir, scseq, gene)

      # update base feature plot
      plot <- feature_plot()
      plot <- update_feature_plot(plot, gene_data, gene)
      if (!show_dprimes() | !is_integrated()) {
        plot <- plot_feature_sample(gene, scseq_samp(), 'ctrl', plot=plot)
      }
      pfun <- list(plot = plot, height = 453)
      return(pfun)
    }
  })

  pfun_right_bottom <- reactive({

    function(gene) {
      sel <- sel()
      scseq <- scseq()
      is_integrated <- is_integrated()
      default <- list(plot=NULL, height=1)

      if(!isTruthyAll(sel, scseq, gene, is_integrated)) return(default)

      dat_path <- file.path(plots_dir(), paste(gene,  sel, 'sample_ridgedat.qs', sep='-'))
      if (file.exists(dat_path)) {
        ridge_data <- qs::qread(dat_path)
        ridge_data$clus_levs <- annot()


      } else {
        try(ridge_data <- get_ridge_data(gene, scseq, sel, by.sample = TRUE, with_all = TRUE))
        if (is.null(ridge_data)) return(default)
        qs::qsave(ridge_data, dat_path)
      }


      plot <- VlnPlot(ridge_data = ridge_data, with.height = TRUE, is_mobile = is_mobile())
      return(plot)
    }
  }) %>% debounce(20)

  dataset_ambient <- reactive(qread.safe(file.path(dataset_dir(), 'ambient.qs')))

  # differential expression top tables for all clusters
  top_tables <- reactive({

    resoln_dir <- resoln_dir()
    if (is.null(resoln_dir)) return(NULL)
    tts_path <- file.path(resoln_dir, 'top_tables.qs')

    if (file.exists(tts_path)) {
      tts <- qs::qread(tts_path)

    } else {
      fit <- lm_fit()
      dataset_ambient <- dataset_ambient()
      if (is.null(fit) | is.null(dataset_ambient)) return(NULL)

      tts <- list()
      for (i in seq_along(fit)) {

        cluster <- names(fit)[i]
        tt <- crossmeta::get_top_table(fit[[cluster]])
        markers <- get_cluster_markers(cluster, resoln_dir)
        ambient <- decide_ambient(dataset_ambient, tt, markers)
        tt$ambient <- row.names(tt) %in% ambient

        # add ambient-excluded adjusted pvals
        tt$adj.P.Val.Amb[!tt$ambient] <- stats::p.adjust(tt$P.Value[!tt$ambient], method = 'BH')
        if (all(tt$ambient)) tt$adj.P.Val.Amb <- NA

        tts[[cluster]] <- tt
      }

      # add 'All Clusters' result
      annot <-  qs::qread(file.path(resoln_dir, 'annot.qs'))
      all <- as.character(length(annot)+1)
      es <- run_esmeta(tts)
      enids <- extract_enids(tts)
      cols <- colnames(tts[[1]])
      tts[[all]] <- es_to_tt(es, enids, cols)

      qs::qsave(tts, file.path(resoln_dir, 'top_tables.qs'))
    }
    return(tts)
  })


  # top table for selected cluster only
  top_table <- reactive({
    sel <- input$selected_cluster
    if (!isTruthy(sel)) return(NULL)

    top_tables()[[sel]]
  })


  # need ambient for pathway
  ambient <- reactive({tt <- top_table() ; row.names(tt)[tt$ambient]})

  # indicating if 'All Clusters' selected
  is.meta <- reactive(input$selected_cluster == tail(names(top_tables()), 1))

  # drug query results
  drug_queries <- reactive({
    dpaths <- drug_paths()
    saved_drugs <- any(grepl('^cmap_res_', list.files(resoln_dir())))

    if (!isTruthy(input$selected_cluster)) {
      res <- NULL

    } else if (saved_drugs) {
      if (!file.exists(dpaths$cmap)) res <- NULL
      else res <- lapply(dpaths, qs::qread)

    } else {
      # run for all single-cluster comparisons (slowest part is loading es)
      disableAll(input_ids)

      progress <- Progress$new(session, min = 0, max = 3)
      progress$set(message = "Querying drugs", value = 1)
      on.exit(progress$close())

      es <- load_drug_es()
      progress$inc(1)
      tts <- top_tables()
      species <- species()

      for (i in seq_along(tts)) {
        cluster <- names(tts)[i]
        tt <- tts[[cluster]]
        ambient <- row.names(tt)[tt$ambient]
        paths <- get_drug_paths(resoln_dir(), cluster)
        run_drug_queries(tt, paths, es, ambient, species)
      }

      progress$inc(1)
      enableAll(input_ids)
      if (!file.exists(dpaths$cmap)) res <- NULL
      else res <- lapply(dpaths, qs::qread)

    }
    return(res)

  })


  # goana pathway result
  path_res <- reactive({
    go_path <- go_path()
    kegg_path <- kegg_path()
    goana_path <- goana_path()
    kegga_path <- kegga_path()

    if (file.exists(kegga_path)) {
      res <- list(go = qs::qread(go_path),
                  kg = qs::qread(kegg_path),
                  kegga = qs::qread(kegga_path),
                  goana = qs::qread(goana_path))

    } else {
      lm_fit <- lm_fit()
      disableAll(input_ids)
      progress <- Progress$new(session, min = 0, max = 2)
      progress$set(message = "Running pathway analysis", value = 1)
      on.exit(progress$close())

      cluster <- input$selected_cluster
      lm_fit <- lm_fit[[cluster]]

      species <- scseq()@metadata$species
      species <- paste(substring(strsplit(species, ' ')[[1]], 1, 1), collapse = '')

      if (is.meta()) {
        de <- top_table()

      } else {
        # exclude ambient
        ambient <- ambient()
        lm_fit <- within(lm_fit, fit <- fit[!row.names(fit) %in% ambient, ])
        de <- crossmeta::fit_ebayes(lm_fit, 'test-ctrl')
      }

      res <- get_path_res(de,
                          go_path = go_path,
                          kegg_path = kegg_path,
                          goana_path = goana_path,
                          kegga_path = kegga_path,
                          species = species)
      progress$inc(1)
      enableAll(input_ids)
    }

    return(res)
  })

  pairs <- reactive(qread.safe(file.path(dataset_dir(), 'pairs.qs')))

  abundances <- reactive(diff_abundance(scseq(), annot(), pairs()))



  # enable download
  observe({
    toggleState('download', condition = isTruthy(top_table()))
  })

  annot_clusters <- reactive({
    clusts <- input$selected_cluster
    clusts <- as.numeric(clusts)
    req(clusts)

    annot <- gsub(' ', '-', annot())
    clusts <- paste0(c(annot, 'all')[sort(clusts)], collapse = '_')
    return(clusts)
  })

  # name for  downloading
  dl_fname <- reactive({
    date <- paste0(Sys.Date(), '.zip')
    clusts <- annot_clusters()
    snn <- basename(resoln_dir())
    paste('single-cell', dataset_name(), clusts, snn, date , sep='_')
  })

  data_fun <- function(file) {
    #go to a temp dir to avoid permission issues
    owd <- setwd(tempdir())
    on.exit(setwd(owd))

    tt_fname <- 'top_table.csv'
    go_fname <- 'cameraPR_go.csv'
    kg_fname <- 'cameraPR_kegg.csv'
    ab_fname <- 'abundances.csv'
    kegga_fname <- 'kegga.csv'
    goana_fname <- 'goana.csv'

    tt <- top_table()
    if (is.meta()) tt <- tt_to_es(tt)

    pres <- path_res()
    abundances <- abundances()
    tozip <- c()
    tozip <- write.csv.safe(tt, tt_fname, tozip)
    tozip <- write.csv.safe(tt, tt_fname, tozip)
    tozip <- write.csv.safe(pres$go, go_fname, tozip)
    tozip <- write.csv.safe(pres$kg, kg_fname, tozip)
    tozip <- write.csv.safe(pres$kegga, kegga_fname, tozip)
    tozip <- write.csv.safe(pres$goana, goana_fname, tozip)
    tozip <- write.csv.safe(abundances, ab_fname, tozip)

    #create the zip file
    utils::zip(file, tozip)
  }

  output$download <- downloadHandler(
    filename = function() {
      dl_fname()
    },
    content = data_fun
  )

  # download can timeout so get objects before clicking
  observeEvent(input$click_dl, {
    tt <- top_table()
    pres <- path_res()
    shinyjs::click("download")
  })

  observe(toggleState('click_dl', condition = isTruthy(input$selected_cluster)))

  selected_cluster <- reactiveVal()
  observe(selected_cluster(input$selected_cluster))


  return(list(
    ambient = ambient,
    top_table = top_table,
    has_replicates = has_replicates,
    cluster_markers = cluster_markers,
    cluster_choices = cluster_choices,
    drug_queries = drug_queries,
    path_res = path_res,
    selected_cluster = selected_cluster,
    annot_clusters = annot_clusters,
    pfun_left = pfun_left,
    pfun_right = pfun_right,
    pfun_right_bottom = pfun_right_bottom
  ))
}
