#' Logic for Single Cell Exploration page
#' @export
#' @keywords internal
scPage <- function(input, output, session, sc_dir, indices_dir) {

  # the analysis and options
  scForm <- callModule(scForm, 'form',
                       sc_dir = sc_dir,
                       indices_dir = indices_dir)

  # cluster plot in top right
  callModule(scClusterPlot, 'cluster_plot',
             scseq = scForm$scseq,
             selected_cluster = scForm$selected_cluster,
             dataset_name = scForm$dataset_name,
             downloadable = TRUE)

  # cluster comparison plots ---

  callModule(scMarkerPlot, 'marker_plot_cluster',
             scseq = scForm$scseq,
             custom_metrics = scForm$custom_metrics,
             selected_feature = scForm$clusters_gene,
             dataset_name = scForm$dataset_name)

  callModule(scBioGpsPlot, 'biogps_plot',
             selected_gene = scForm$clusters_gene,
             species = scForm$species)


  callModule(scRidgePlot, 'ridge_plot',
             selected_gene = scForm$clusters_gene,
             selected_cluster = scForm$clusters_cluster,
             scseq = scForm$scseq)


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



  # label comparison plot ---
  callModule(scLabelsPlot, 'labels_plot_cluster',
             selected_cluster = scForm$labels_cluster,
             scseq = scForm$scseq,
             sc_dir = sc_dir)

  observe({
    toggle(id = "comparison_row",  condition = isTruthy(scForm$dataset_name()))
  })

  observe({
    toggle(id = "sample_comparison_row",  condition = scForm$comparison_type() == 'samples')
    toggle(id = "cluster_comparison_row", condition = scForm$comparison_type() == 'clusters')
    toggle(id = "labels_comparison_row", condition = scForm$comparison_type() == 'labels')
  })


  observe({
    toggle(id = 'biogps_container', condition = !scForm$show_ridge())
    toggle(id = 'ridge_container', condition = scForm$show_ridge())
  })

  return(NULL)
}


#' Logic for form on Single Cell Exploration page
#' @export
#' @keywords internal
scForm <- function(input, output, session, sc_dir, indices_dir) {

  # updates if new integrated or subset dataset
  new_dataset <- reactiveVal()
  observe(new_dataset(scIntegration()))
  observe(new_dataset(scSubset()))

  # the dataset and options
  scDataset <- callModule(scSelectedDataset, 'dataset',
                          sc_dir = sc_dir,
                          new_dataset = new_dataset,
                          indices_dir = indices_dir)

  observe(toggle('form_container', condition = scDataset$dataset_exists()))


  # update scseq with annotation changes and custom metrics
  scseq <- reactive({
    scseq <- scDataset$scseq()
    annot <- scClusterComparison$annot()
    metrics <- scClusterGene$saved_metrics()

    if (!isTruthy(annot) | !isTruthy(scseq)) return(NULL)
    if (!is.null(metrics)) scseq@colData <- cbind(scseq@colData, metrics)
    levels(scseq$cluster) <- annot

    return(scseq)
  })

  qc_metrics <- reactive({
    scseq <- scseq()
    if(is.null(scseq)) return(NULL)
    metrics <- scseq@colData
    qc <- colnames(metrics)

    names(qc) <- sapply(metrics, class)

    qc <- qc[names(qc) %in% c('numeric', 'logical')]
    qc <- qc[!grepl('^sum$|^total$|^subsets|^percent', qc)]
  })


  #TODO move label transfer and integration into dataset

  # label transfer between datasets
  scLabelTransfer <- callModule(labelTransferForm, 'transfer',
                                sc_dir = sc_dir,
                                datasets = scDataset$datasets,
                                show_label_transfer = scDataset$show_label_transfer,
                                dataset_name = scDataset$dataset_name,
                                scseq = scseq,
                                species = scDataset$species)

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
                         dataset_dir = dataset_dir,
                         show_subset = scDataset$show_subset,
                         is_integrated = scDataset$is_integrated)


  # comparison type
  comparisonType <- callModule(comparisonType, 'comparison',
                               scseq = scseq,
                               is_integrated = scDataset$is_integrated)


  # the selected cluster/gene for cluster comparison
  dataset_dir <- reactive({
    dataset <- scDataset$dataset_name()
    if (is.null(dataset)) return(NULL)
    file.path(sc_dir, dataset)
  })

  scClusterComparison <- callModule(clusterComparison, 'cluster',
                                    dataset_dir = dataset_dir,
                                    scseq = scseq,
                                    annot_path = scDataset$annot_path,
                                    ref_preds = scLabelTransfer)

  scClusterGene <- callModule(selectedGene, 'gene_clusters',
                              scseq = scDataset$scseq,
                              dataset_name = scDataset$dataset_name,
                              dataset_dir = dataset_dir,
                              is_integrated = scDataset$is_integrated,
                              selected_markers = scClusterComparison$selected_markers,
                              selected_cluster = scClusterComparison$selected_cluster,
                              qc_metrics = qc_metrics,
                              type = 'clusters')



  # the selected clusters/gene for sample comparison
  scSampleComparison <- callModule(scSampleComparison, 'sample',
                                   dataset_dir = dataset_dir,
                                   dataset_name = scDataset$dataset_name,
                                   is_integrated = scDataset$is_integrated,
                                   input_scseq = scseq,
                                   comparison_type = comparisonType,
                                   exclude_ambient = scSampleGene$exclude_ambient)

  scSampleGene <- callModule(selectedGene, 'gene_samples',
                             scseq = scDataset$scseq,
                             dataset_name = scDataset$dataset_name,
                             dataset_dir = dataset_dir,
                             is_integrated = scDataset$is_integrated,
                             selected_markers = scSampleComparison$top_table,
                             selected_cluster = scSampleComparison$selected_cluster,
                             type = 'samples',
                             ambient = scSampleComparison$ambient)



  scLabelsComparison <- callModule(scLabelsComparison, 'labels',
                                   cluster_choices = scSampleComparison$cluster_choices)


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
    species = scDataset$species
  ))
}


#' Logic for selected dataset part of scForm
#' @export
#' @keywords internal
scSelectedDataset <- function(input, output, session, sc_dir, new_dataset, indices_dir) {
  dataset_inputs <- c('selected_dataset', 'show_integration', 'show_label_transfer')
  options <- list(render = I('{option: scDatasetOptions, item: scDatasetItem}'))

  # get directory with fastqs
  roots <- c('single-cell' = sc_dir)
  shinyFiles::shinyDirChoose(input, "new_dataset_dir", roots = roots)

  dataset_exists <- reactive(isTruthy(input$selected_dataset) & !is.create())

  dataset_name <- reactive({
    if (!dataset_exists()) return(NULL)
    ds <- datasets()
    ds$name[ds$value == input$selected_dataset]
  })


  # get's used for saving annotation to disc
  annot_path <- reactive({
    scseq_part_path(sc_dir, dataset_name(), 'annot')
  })

  # load annotation for clusters
  annot <- reactive(readRDS(annot_path()))

  # load scseq
  scseq <- reactive({
    if (!isTruthy(dataset_name())) return(NULL)
    disableAll(dataset_inputs)

    dataset_dir <- file.path(sc_dir, dataset_name())
    scseq <- load_scseq(dataset_dir)

    enableAll(dataset_inputs)
    return(scseq)
  })

  is_integrated <- reactive({
    dataset_name <- dataset_name()
    req(dataset_name)
    integrated <- readRDS(file.path(sc_dir, 'integrated.rds'))
    return(dataset_name %in% integrated)
  })

  species <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)
    scseq@metadata$species
  })

  # available single-cell datasets
  datasets <- reactive({
    # reactive to new single cell datasets
    new_dataset()
    get_sc_dataset_choices(sc_dir)
  })


  # update previously selected dataset on-file if changes
  prev_path <- file.path(sc_dir, 'prev_dataset.rds')

  observe({
    sel <- dataset_name()
    req(sel)
    saveRDS(sel, prev_path)
  })


  # are we creating a new dataset?
  is.create <- reactive({
    dataset_name <- input$selected_dataset
    datasets <- datasets()
    if (!isTruthy(dataset_name)) return(FALSE)

    !dataset_name %in% datasets$value
  })

  # open shinyFiles selector if creating
  observe({
    req(is.create())
    shinyjs::click('new_dataset_dir')
  })

  # get path to dir with new dataset files
  new_dataset_dir <- reactive({

    new_dataset_dir <- input$new_dataset_dir

    # need selected subfolder
    # will be integer on create
    req(!'integer' %in% class(new_dataset_dir))

    dir <- shinyFiles::parseDirPath(roots, new_dataset_dir)
    as.character(dir)
  })

  # ask for confirmation after folder selection
  observeEvent(new_dataset_dir(), showModal(quantModal()))


  metric_choices <- c('low_lib_size',
                      'low_n_features',
                      'high_subsets_mito_percent',
                      'low_subsets_ribo_percent',
                      'high_doublet_score')

  # run single-cell quantification
  observeEvent(input$confirm_quant, {

    metrics <- input$qc_metrics
    # none, all, all and none: can't combine
    if (length(metrics) > 1 && !all(metrics %in% metric_choices)) return(NULL)

    if (!isTruthy(metrics)) metrics <- 'none'
    if (metrics == 'none') metrics <- NULL
    if (metrics == 'all') metrics <- metric_choices

    removeModal()
    disableAll(dataset_inputs)

    # Create a Progress object
    progress <- Progress$new(session, min=0, max = 9)
    on.exit(progress$close())

    fastq_dir <- new_dataset_dir()
    dataset_name <- input$selected_dataset

    if (metrics == 'all and none') {
      load_raw_scseq(paste0(dataset_name, '_QC0'), fastq_dir, sc_dir, indices_dir, progress, metrics = NULL)
      load_raw_scseq(paste0(dataset_name, '_QC1'), fastq_dir, sc_dir, indices_dir, progress, metrics = metric_choices, founder = paste0(dataset_name, '_QC0'))

      prev <- paste0(dataset_name, '_QC1')

    } else {
      load_raw_scseq(dataset_name, fastq_dir, sc_dir, indices_dir, progress, metrics = metrics)
      prev <- dataset_name
    }

    saveRDS(prev, prev_path)

    enableAll(dataset_inputs)
    new_dataset(dataset_name)
  })

  # modal to confirm adding single-cell dataset
  quantModal <- function() {
    UI <- selectizeInput(session$ns('qc_metrics'),
                         'Select QC metrics:',
                         choices = c('all and none', 'all', 'none', metric_choices),
                         selected = 'all and none',
                         multiple = TRUE)

    modalDialog(
      UI,
      title = 'Create new single-cell dataset?',
      size = 's',
      footer = tagList(
        modalButton("Cancel"),
        actionButton(session$ns("confirm_quant"), "Quantify", class = 'pull-left btn-warning')
      )
    )
  }


  observe({
    updateSelectizeInput(session, 'selected_dataset', choices =  rbind(NA, datasets()), server = TRUE, options = options)
  })

  # show/hide integration/label-transfer forms
  show_integration <- reactive(input$show_integration %% 3 == 2)
  show_subset <- reactive(input$show_integration %% 3 == 1)
  show_label_transfer <- reactive(input$show_label_transfer %% 2 != 0)

  observe(toggleClass(id = "show_label_transfer", 'btn-primary', condition = show_label_transfer()))

  # distinguish between 1/2 toggles for integration vs subset
  observe({
    toggleClass(id = "show_integration", 'btn-primary', condition = show_integration() | show_subset())
  })

  observe({
    icon <- icon('object-ungroup', 'far fa-fw')
    if (show_integration()) icon <- icon('object-group', 'far fa-fw')

    updateActionButton(session, 'show_integration', icon = icon)
  })


  # hide integration/label-transfer buttons no dataset
  observe({
    toggle('show_label_transfer-parent', condition = dataset_exists())
  })


  return(list(
    dataset_name = dataset_name,
    scseq = scseq,
    annot = annot,
    annot_path = annot_path,
    datasets = datasets,
    show_integration = show_integration,
    show_subset = show_subset,
    show_label_transfer = show_label_transfer,
    is_integrated = is_integrated,
    dataset_exists = dataset_exists,
    species = species
  ))
}


#' Logic for selecting cluster to plot label origin for integrated dataset
#' @export
#' @keywords internal
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

#' Logic for labels plot for integrated dataset
#' @export
#' @keywords internal
scLabelsPlot <- function(input, output, session, sc_dir, selected_cluster, scseq) {

  output$labels_plot <- plotly::renderPlotly({
    scseq <- scseq()
    if (is.null(scseq) | is.null(scseq$batch)) return(NULL)

    cluster <- as.numeric(selected_cluster())
    cluster <- levels(scseq$cluster)[cluster]
    if (is.na(cluster)) return(NULL)

    plot_cluster_labels(scseq, cluster, sc_dir)
  })

}


#' Logic for single-cell sample comparison plots
#'
#' setup to allow for ggplot/plotly
#'
#' @export
#' @keywords internal
scSampleMarkerPlot <- function(input, output, session, selected_gene, plot_fun) {

  res <- reactive({
    gene <- selected_gene()
    req(gene)
    suppressMessages(plot_fun()(gene))
  })

  height <- reactive(ifelse(is_plotly(), 1, res()$height))

  is_plotly <- reactive('plotly' %in% class(res()))

  output$plot <- renderPlot(if (!is_plotly()) suppressMessages(print(res()$plot)) else NULL, height = height)
  output$plotly <- plotly::renderPlotly(if (is_plotly()) res() else NULL)
}

#' Logic for label transfer between datasets
#' @export
#' @keywords internal
labelTransferForm <- function(input, output, session, sc_dir, datasets, show_label_transfer, dataset_name, scseq, species) {
  label_transfer_inputs <- c('transfer_study', 'submit_transfer', 'overwrite_annot', 'ref_name')

  ref_preds <- reactiveVal()
  new_preds <- reactiveVal()
  new_annot <- reactiveVal()

  # show/hide label transfer forms
  observe({
    toggle(id = "label-transfer-form", anim = TRUE, condition = show_label_transfer())
  })

  # saved label transfer predictions
  preds <- reactive({
    new_preds()
    query_name <- dataset_name()
    req(query_name)

    # load previously saved reference preds
    preds_path <- scseq_part_path(sc_dir, query_name, 'preds')
    preds <- if (file.exists(preds_path)) readRDS(preds_path) else list()
    preds <- validate_preds(preds, sc_dir)
    return(preds)
  })

  observeEvent(dataset_name(), new_preds(NULL))

  # update annotation transfer choices
  observe({
    preds <- preds()

    datasets <- datasets()
    dataset_name <- dataset_name()
    species <- species()
    req(preds, datasets, species)

    choices <- get_label_transfer_choices(datasets, dataset_name, preds, species)
    updateSelectizeInput(session, 'ref_name', choices = choices, server = TRUE, selected = isolate(new_preds()), options = list(render = I('{option: transferLabelOption, item: scDatasetItem}')))
  })


  # submit annotation transfer
  observeEvent(input$ref_name, {

    query_name <- dataset_name()
    ref_name <- input$ref_name
    preds <- preds()

    req(ref_name != 'reset')
    req(query_name, ref_name, preds)
    req(!ref_name %in% names(preds))
    req(show_label_transfer())

    query <- scseq()
    req(query)

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
    senv <- loadNamespace('SingleR')
    updateProgress(1/n)

    if (ref_name %in% ls(senv)) {
      ref <- get(ref_name, envir = senv)()
      labels <- ref$label.main
      genes <- 'de'

    } else {
      ref_path <- scseq_part_path(sc_dir, ref_name, 'scseq')
      ref_date <- file.info(ref_path)$ctime
      ref <- readRDS(ref_path)

      # check if ref and query have the same founder
      rfound <- readRDS(scseq_part_path(sc_dir, ref_name, 'founder'))
      qfound <- readRDS(scseq_part_path(sc_dir, query_name, 'founder'))

      # use common cells to transfer labels if so
      cells <- intersect(colnames(ref), colnames(query))

      if (identical(qfound, rfound) && length(cells)) {

        ref_cluster <- ref[, cells]$cluster
        query_cluster <- query[, cells]$cluster
        tab <- table(assigned = ref_cluster, cluster = query_cluster)

      } else {
        markers_path <- scseq_part_path(sc_dir, ref_name, 'top_markers')
        genes <- readRDS(markers_path)

        # use aggregated reference for speed
        ref_path <- scseq_part_path(sc_dir, ref_name, 'aggr_ref')
        if (file.exists(ref_path)) {
          ref <- readRDS(ref_path)

        } else {
          set.seed(100)
          ref <- SingleR::aggregateReference(ref, labels=ref$cluster)
          saveRDS(ref, ref_path)
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

    # keep track of date that reference was used so that can invaludate if overwritten
    if (!missing(ref_date)) names(pred) <- ref_date

    preds_path <- scseq_part_path(sc_dir, query_name, 'preds')
    preds[[ref_name]] <- pred
    saveRDS(preds, preds_path)

    new_preds(ref_name)
    ref_preds(preds[[ref_name]])

    enableAll(label_transfer_inputs)
  })


  # show transfered labels immediately upon selection if have
  observe({
    query_name <- dataset_name()
    ref_name <- input$ref_name
    req(query_name)

    # load previously saved reference preds
    preds_path <- scseq_part_path(sc_dir, query_name, 'preds')
    preds <- if (file.exists(preds_path)) readRDS(preds_path) else list()

    ref_preds(preds[[ref_name]])
  })



  pred_annot <- reactive({
    # react to new annotation
    new_annot()
    ref_name <- input$ref_name
    ref_preds <- ref_preds()
    dataset_name <- dataset_name()
    req(dataset_name)

    # show saved annot if nothing selected or label transfer not open
    if (!isTruthy(ref_name) | !show_label_transfer()) {
      annot_path <- scseq_part_path(sc_dir, dataset_name, 'annot')
      annot <- readRDS(annot_path)

    } else {
      annot <- get_pred_annot(ref_preds, ref_name, dataset_name, sc_dir)
    }

    return(annot)
  })

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
    dataset_name <- dataset_name()

    req(dataset_name)

    showModal(transferModal())
  })

  observeEvent(input$confirm_overwrite, {
    removeModal()
    ref_name <- input$ref_name
    ref_preds <- ref_preds()
    dataset_name <- dataset_name()

    req(dataset_name)

    pred_annot <- get_pred_annot(ref_preds, ref_name, dataset_name, sc_dir)
    annot_path <- scseq_part_path(sc_dir, dataset_name, 'annot')
    saveRDS(pred_annot, annot_path)

    new_annot(pred_annot)
  })



  return(pred_annot)
}

#' Logic for subsetting a datatset
#' @export
#' @keywords internal
subsetForm <- function(input, output, session, sc_dir, scseq, datasets, show_subset, selected_dataset, dataset_dir, cluster_choices, is_integrated) {
  contrastOptions <- list(render = I('{option: contrastOptions, item: contrastItem}'))

  subset_name <- reactive(input$subset_name)
  new_dataset <- reactiveVal()

  subset_inputs <- c('subset_name',
                     'submit_subset',
                     'subset_clusters',
                     'toggle_exclude')

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
  observeEvent(input$submit_subset, {

    subset <- input$subset_clusters
    cluster_choices <- cluster_choices()
    metric_choices <- metric_choices()
    is_include <- is_include()

    exclude_clusters <- intersect(cluster_choices$value, subset)
    subset_metrics <- intersect(metric_choices$value, subset)

    if (is_include && length(exclude_clusters)) {
      exclude_clusters <- setdiff(cluster_choices$value, exclude_clusters)
    }

    from_dataset <- selected_dataset()
    founder <- get_founder(sc_dir, from_dataset)
    subset_name <- input$subset_name
    dataset_name <- paste(founder, subset_name, sep = '_')

    # clear error and disable button
    disableAll(subset_inputs)

    # Create a Progress object
    on.exit(progress$close())
    progress <- Progress$new(session, min=0, max = 9)

    progress$set(message = "Subsetting dataset", value = 0)
    subset_saved_scseq(sc_dir = sc_dir,
                       founder = founder,
                       from_dataset = from_dataset,
                       dataset_name = dataset_name,
                       exclude_clusters = exclude_clusters,
                       subset_metrics = subset_metrics,
                       is_include = is_include,
                       is_integrated = is_integrated(),
                       progress = progress)


    # re-enable, clear inputs, and trigger update of available datasets
    new_dataset(dataset_name)
    updateTextInput(session, 'subset_name', value = '')
    enableAll(subset_inputs)
  })

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
#' @export
#' @keywords internal
integrationForm <- function(input, output, session, sc_dir, datasets, show_integration, selected_dataset) {
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
    int  <- readRDS.safe(file.path(sc_dir, 'integrated.rds'))
    prev <- readRDS.safe(file.path(sc_dir, 'prev_dataset.rds'))
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

  # update exclude clusters
  observe({
    updateSelectizeInput(session, 'exclude_clusters', choices = exclude_choices(),selected = isolate(input$subset_clusters), options = excludeOptions, server = TRUE)
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
      write.csv(samples(), con, row.names = FALSE)
    }
  )

  observe({
    infile <- input$up_pairs
    req(infile)
    pairs(read.csv(infile$datapath, row.names = 'sample'))
  })

  # make upload green when have data
  observe(toggleClass(id = "click_up", 'btn-success', condition = isTruthy(pairs())))

  # change UI of exclude toggle
  observe({
    toggleClass(id = 'toggle_icon', 'fa-plus text-success', condition = is_include())
    toggleClass(id = 'toggle_icon', 'fa-minus text-warning', condition = !is_include())
  })

  # run integration
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

    error_msg <- validate_integration(test, ctrl, pairs)

    if (is.null(error_msg)) {
      # clear error and disable button
      removeClass('name-container', 'has-error')
      disableAll(integration_inputs)

      # Create a Progress object
      on.exit(progress$close())
      ntype <- length(integration_types)
      progress <- Progress$new(session, min=0, max = 9*ntype)

      for (i in seq_along(integration_types)) {
        vali <- (i-1)*9
        progress$set(message = paste(integration_types[i], 'integration:'), detail = '', value = vali)

        # run integration
        integrate_saved_scseqs(sc_dir,
                               test = test,
                               ctrl = ctrl,
                               integration_name = integration_name,
                               integration_type = integration_types[i],
                               exclude_clusters = exclude_clusters,
                               pairs = pairs,
                               progress = progress,
                               value = vali)
      }



      # re-enable, clear inputs, and trigger update of available datasets
      ctrl(NULL)
      test(NULL)
      pairs(NULL)
      new_dataset(integration_name)
      updateTextInput(session, 'integration_name', value = '')
      enableAll(integration_inputs)

    } else {
      # show error message
      pairs(NULL)
      html('error_msg', html = error_msg)
      addClass('name-container', class = 'has-error')
    }

  })

  return(new_dataset)
}


#' Logic for comparison type toggle for integrated datasets
#' @export
#' @keywords internal
comparisonType <- function(input, output, session, scseq, is_integrated) {

  # always show clusters if not integrated
  observe({
    if(!is_integrated())
      updateRadioGroupButtons(session, 'comparison_type', selected = 'clusters')
  })

  return(reactive(input$comparison_type))
}

#' Logic for cluster comparison input
#' @export
#' @keywords internal
clusterComparison <- function(input, output, session, dataset_dir, scseq, annot_path, ref_preds) {


  contrast_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))
  selected_cluster <- reactiveVal()
  markers <- reactiveVal(list())

  # things that return for plotting
  annot <- reactiveVal(NULL)
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
      choices <- get_cluster_choices(clusters, scseq = scseq())
    }

    return(choices)
  })


  observe({
    ref_preds <- ref_preds()
    if (!is.null(ref_preds)) annot(ref_preds)
    else annot(readRDS(annot_path()))
  })

  # set to first cluster if switch to showing contrasts
  observeEvent(input$show_contrasts, {
    if (show_contrasts()) {
      selected_cluster(NULL)
    } else {
      selected_cluster(test_cluster())
    }
  })

  observeEvent(input$selected_cluster, {
    selected_cluster(input$selected_cluster)
  })

  # reset if switch dataset
  observeEvent(dataset_dir(), {
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
    mod_annot <- readRDS(annot_path())
    mod_annot[sel_idx] <- ref_preds[sel_idx] <- input$new_cluster_name
    mod_annot <- make.unique(mod_annot, '_')

    # save on disc
    saveRDS(mod_annot, annot_path())

    # update annot and set selected cluster to new name
    annot(mod_annot)
  })


  # update UI for contrast/cluster choices
  observeEvent(choices(), {

    updateSelectizeInput(session, 'selected_cluster',
                         choices = rbind(NA, choices()),
                         selected = isolate(selected_cluster()),
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
    dataset_dir <- dataset_dir()
    markers <- markers()
    req(sel, dataset_dir)

    if (sel %in% names(markers)) return(NULL)


    if (show_contrasts()) {
      con_markers <- get_contrast_markers(sel, dataset_dir)
      markers <- c(markers, con_markers)

    } else {
      markers_path <- file.path(dataset_dir, paste0('markers_', sel, '.rds'))
      markers[[sel]] <- readRDS(markers_path)
    }

    markers(markers)
  })


  observe({
    sel <- selected_cluster()
    if (!isTruthy(sel)) selected_markers(NULL)
    else selected_markers(markers()[[sel]])
  })


  return(list(
    annot = annot,
    selected_markers = selected_markers,
    selected_cluster = selected_cluster
  ))
}


#' Logic for selected gene to show plots for
#' @export
#' @keywords internal
selectedGene <- function(input, output, session, dataset_name, dataset_dir, scseq, is_integrated, selected_markers, selected_cluster, type, qc_metrics = function()NULL, ambient = function()NULL) {
  gene_options <- list(render = I('{option: geneChoice, item: geneChoice}'))

  selected_gene <- reactiveVal(NULL)

  exclude_ambient <- reactive({
    if (is.null(input$exclude_ambient)) return(TRUE)
    input$exclude_ambient %% 2 != 1
  })

  # toggle for excluding ambient
  observe({
    toggleClass('exclude_ambient', class = 'btn-primary', condition = exclude_ambient())
  })

  # toggle for ridgeline
  gene_selected <- reactive({
    sel <- selected_gene()
    scseq <- scseq()
    if (!isTruthyAll(sel, scseq)) return(FALSE)
    sel %in% row.names(scseq)
  })

  show_ridge <- reactive(input$show_ridge %% 2 != 1 | !gene_selected())
  show_custom_metric <- reactive(type != 'samples' && (input$show_custom_metric %%2 != 0))

  observe(toggleClass(id = "show_ridge", 'btn-primary', condition = !show_ridge()))

  observe({
    toggle('custom_metric_panel', anim = TRUE, condition = show_custom_metric())
    if (show_custom_metric() & have_metric()) selected_gene(input$custom_metric)
  })

  observe(if (!show_custom_metric()) selected_gene(isolate(input$selected_gene)))
  observe(toggleClass('show_custom_metric', class = 'btn-primary', condition = show_custom_metric()))

  observe({
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

    if (class(res) == 'data.frame') {

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
  metrics_path <- reactive(file.path(dataset_dir(), 'saved_metrics.rds'))
  observe(saved_metrics(readRDS.safe(metrics_path())))

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

    saveRDS(res, metrics_path())
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
    get_gene_choices(markers, qc_metrics = qc_metrics, qc_first = qc_first, species = species())
  })

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


  # update choices
  observe({
    choices <- gene_choices()
    if(is.null(choices)) return(NULL)
    selected <- isolate(selected_gene())

    updateSelectizeInput(session, 'selected_gene',
                         choices = rbind(NA, choices, fill = TRUE),
                         selected = selected,
                         server = TRUE,
                         options = gene_options)
  })


  return(list(
    selected_gene = selected_gene,
    exclude_ambient = exclude_ambient,
    show_ridge = show_ridge,
    custom_metrics = custom_metrics,
    saved_metrics = saved_metrics
  ))


}


#' Logic for cluster plots
#' @export
#' @keywords internal
scClusterPlot <- function(input, output, session, scseq, selected_cluster, dataset_name, downloadable = FALSE) {


  fname_fun <- function() {
    paste0(dataset_name(), '_cluster_plot_data_', Sys.Date(), '.csv')
  }

  plot <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)

    label.highlight <- NULL
    cluster <- selected_cluster()
    cluster <- strsplit(cluster, '-vs-')[[1]]
    if (isTruthy(cluster)) label.highlight <- as.numeric(cluster)

    plot_tsne_cluster(scseq, label.highlight = label.highlight)
  })


  data_fun <- function() {
    plot <- plot()
    plot$data
  }

  callModule(downloadablePlot,
             "cluster_plot",
             plot_fun = plot,
             fname_fun = fname_fun,
             data_fun = data_fun,
             downloadable = downloadable)


  return(list(
    plot = plot
  ))
}

#' Logic for plot with downloadable data
#' @export
#' @keywords internal
downloadablePlot <- function(input, output, session, plot_fun, fname_fun, data_fun, downloadable) {

  # show download button only if plot visible and downloadable
  observe({
    toggleClass('download_container', class = 'visible-plot', condition = downloadable)
  })




  # click download
  output$download <- downloadHandler(
    filename = fname_fun,
    content = function(con) {
      write.csv(data_fun(), con)
    }
  )

  output$dl_plot <- renderPlot({
    plot_fun()
  }, height = 'auto')
}

#' Logic for marker feature plots
#' @export
#' @keywords internal
scMarkerPlot <- function(input, output, session, scseq, selected_feature, dataset_name, downloadable = TRUE, custom_metrics = function()NULL) {

  fname_fun <- function() {
    fname <- paste0(dataset_name(),
                    '_', selected_feature(),
                    '_marker_plot_data_', Sys.Date(), '.csv')
    return(fname)
  }

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
      pl <- plot_tsne_feature(scseq, feature)

    } else if (is_log) {
      ft <- cdata[[feature]]
      scseq$cluster <- factor(ft, levels = c(FALSE, TRUE))
      ncells <- sum(ft)
      pcells <- round(ncells / length(ft) * 100)

      title <- paste0(feature, ' (', format(ncells, big.mark=","), ' :: ', pcells, '%)')
      pl <- plot_tsne_cluster(scseq, label = FALSE, label.index = FALSE, order = TRUE, title = title)
    } else {
      pl <- NULL
    }

    return(pl)
  })


  data_fun <- function() {
    plot <- plot()
    plot$data
  }

  callModule(downloadablePlot,
             "marker_plot",
             plot_fun = plot,
             fname_fun = fname_fun,
             data_fun = data_fun,
             downloadable = downloadable)




  return(list(
    plot = plot
  ))
}

#' Logic for BioGPS plot
#' @export
#' @keywords internal
scBioGpsPlot <- function(input, output, session, selected_gene, species) {
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
#' @export
#' @keywords internal
scRidgePlot <- function(input, output, session, selected_gene, selected_cluster, scseq) {

  height <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) height <- 453
    else height <- length(levels(scseq$cluster))*50
    return(height)
  })

  output$ridge_plot <- renderPlot({
    gene <- selected_gene()
    cluster <- selected_cluster()
    req(gene)
    suppressMessages(print(plot_ridge(gene, scseq(), cluster)))
  }, height = height)
}


#' Logic for single cell cluster analyses for Single Cell, Drugs, and Pathways tabs
#' @export
#' @keywords internal
scSampleComparison <- function(input, output, session, dataset_dir, dataset_name, is_integrated = function()TRUE, input_scseq = function()NULL, is_sc = function()TRUE, exclude_ambient = function()FALSE, comparison_type = function()'samples') {
  contrast_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))
  input_ids <- c('download', 'selected_cluster')


  # use input scseq/annot if available
  scseq <- reactive({
    scseq <- input_scseq()
    if (is.null(scseq)) scseq <- load_scseq(dataset_dir())
    return(scseq)
  })


  annot <- reactive({
    req(is_sc())
    annot <- levels(input_scseq()$cluster)
    if (is.null(annot)) annot <- readRDS.safe(file.path(dataset_dir(), 'annot.rds'))
    return(annot)
  })

  has_replicates <- reactive(readRDS.safe(file.path(dataset_dir(), 'has_replicates.rds')))
  species <- reactive(readRDS.safe(file.path(dataset_dir(), 'species.rds')))


  # update cluster choices in UI
  cluster_choices <- reactive({
    req(is_integrated())
    dataset_dir <- dataset_dir()
    if (is.null(dataset_dir)) return(NULL)


    get_cluster_choices(clusters = annot(),
                        sample_comparison = TRUE,
                        dataset_dir = dataset_dir,
                        use_disk = TRUE,
                        top_tables = top_tables(),
                        has_replicates = has_replicates())
  })

  observe({
    updateSelectizeInput(session, 'selected_cluster',
                         choices = rbind(NA, cluster_choices()),
                         options = contrast_options, server = TRUE)
  })


  # path to lmfit, cluster markers, drug query results, and goanna pathway results
  clusters_str <- reactive(collapse_sorted(input$selected_cluster))
  drug_paths <- reactive(get_drug_paths(dataset_dir(), clusters_str()))
  go_path <- reactive(file.path(dataset_dir(), paste0('go_', clusters_str(), '.rds')))
  kegg_path <- reactive(file.path(dataset_dir(), paste0('kegg_', clusters_str(), '.rds')))
  top_tables_paths <- reactive(file.path(dataset_dir(), 'top_tables.rds'))


  # require cluster markers and fit result
  lm_fit <- reactive(readRDS.safe(file.path(dataset_dir(), 'lm_fit_0svs.rds')))
  cluster_markers <- reactive(readRDS.safe(file.path(dataset_dir(), paste0('markers_', clusters_str(), '.rds'))))

  # plot functions for left
  sel <- reactive(input$selected_cluster)
  pfun_left <- reactive({
    req(is_integrated())

    function(gene) {
      if(!isTruthy(gene)) return(NULL)

      if (has_replicates()) {
        if(is.null(top_tables())) return(NULL)
        pfun <- plot_scseq_dprimes(gene, annot(), sel(), top_tables(), exclude_ambient())

      } else {
        if(is.null(scseq())) return(NULL)
        pfun <- list(plot = plot_tsne_feature_sample(gene, scseq(), 'test'), height = 453)
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
      if (has_replicates()) pfun <- list(plot = plot_tsne_feature(scseq, gene), height = 453)
      else pfun <- list(plot = plot_tsne_feature_sample(gene, scseq, 'ctrl'), height = 453)
      return(pfun)
    }
  })
  pfun_right_bottom <- reactive({
    req(is_integrated())

    function(gene) {
      scseq <- scseq(); sel <- sel()
      if(!isTruthyAll(sel, scseq, gene)) return(list(plot=NULL, height=1))
      plot_ridge(gene, scseq, sel, by.sample = TRUE, with.height = TRUE)
    }
  })

  dataset_ambient <- reactive(readRDS.safe(file.path(dataset_dir(), 'ambient.rds')))

  # differential expression top tables for all clusters
  top_tables <- reactive({

    dataset_dir <- dataset_dir()
    if (is.null(dataset_dir)) return(NULL)
    tts_path <- file.path(dataset_dir, 'top_tables.rds')

    if (file.exists(tts_path)) {
      tts <- readRDS(tts_path)

    } else {
      fit <- lm_fit()
      dataset_ambient <- dataset_ambient()
      if (is.null(fit) | is.null(dataset_ambient)) return(NULL)

      tts <- list()
      for (i in seq_along(fit)) {

        cluster <- names(fit)[i]
        tt <- get_top_table(fit[[cluster]])
        markers <- get_cluster_markers(cluster, dataset_dir)
        ambient <- decide_ambient(dataset_ambient, tt, markers)
        tt$ambient <- row.names(tt) %in% ambient

        # add ambient-excluded adjusted pvals
        tt$adj.P.Val.Amb[!tt$ambient] <- p.adjust(tt$P.Value[!tt$ambient], method = 'BH')
        tts[[cluster]] <- tt
      }

      saveRDS(tts, file.path(dataset_dir(), 'top_tables.rds'))
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

  # drug query results
  drug_queries <- reactive({
    dpaths <- drug_paths()
    saved_drugs <- any(grepl('^cmap_res_', list.files(dataset_dir())))

    if (!isTruthy(input$selected_cluster)) {
      res <- NULL

    } else if (saved_drugs) {
      if (!file.exists(dpaths$cmap)) res <- NULL
      else res <- lapply(dpaths, readRDS)

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
        paths <- get_drug_paths(dataset_dir(), cluster)
        run_drug_queries(tt, paths, es, ambient, species)
      }

      progress$inc(1)
      enableAll(input_ids)
      if (!file.exists(dpaths$cmap)) res <- NULL
      else res <- lapply(dpaths, readRDS)

    }
    return(res)

  })


  # goana pathway result
  path_res <- reactive({
    go_path <- go_path()
    kegg_path <- kegg_path()

    if (file.exists(kegg_path)) {
      res <- list(go = readRDS(go_path),
                  kg = readRDS(kegg_path))

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

      # exclude ambient
      ambient <- ambient()
      lm_fit <- within(lm_fit, fit <- fit[!row.names(fit) %in% ambient, ])
      ebfit <- fit_ebayes(lm_fit, 'test-ctrl')
      res <- get_path_res(ebfit, go_path, kegg_path, species)
      progress$inc(1)
      enableAll(input_ids)
    }

    return(res)
  })

  pairs <- reactive(readRDS.safe(file.path(dataset_dir(), 'pairs.rds')))

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
    clusts <- paste0(annot[sort(clusts)], collapse = '_')
    return(clusts)
  })

  # name for  downloading
  dl_fname <- reactive({
    date <- paste0(Sys.Date(), '.zip')
    clusts <- annot_clusters()
    paste('single-cell', dataset_name(), clusts, date , sep='_')
  })

  data_fun <- function(file) {
    #go to a temp dir to avoid permission issues
    owd <- setwd(tempdir())
    on.exit(setwd(owd))

    tt_fname <- 'top_table.csv'
    go_fname <- 'go.csv'
    kg_fname <- 'kegg.csv'
    ab_fname <- 'abundances.csv'

    tt <- top_table()
    path_res <- path_res()
    abundances <- abundances()
    write.csv(tt, tt_fname)
    write.csv(path_res$go, go_fname)
    write.csv(path_res$kg, kg_fname)
    write.csv(abundances, ab_fname)

    #create the zip file
    zip(file, c(tt_fname, go_fname, kg_fname, ab_fname))
  }

  output$download <- downloadHandler(
    filename = function() {
      dl_fname()
    },
    content = data_fun
  )

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
