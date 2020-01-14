#' Logic for Single Cell Exploration page
#' @export
#' @keywords internal
scPage <- function(input, output, session, sc_dir, indices_dir) {

  # the analysis and options
  scForm <- callModule(scForm, 'form',
                       sc_dir = sc_dir,
                       indices_dir = indices_dir)


  cluster_data_fname <- function() {
    paste0(scForm$dataset_name(), '_cluster_plot_data_', Sys.Date(), '.csv')
  }

  scCluster <- callModule(scClusterPlot, 'cluster_plot',
                          scseq = scForm$scseq,
                          fname_fun = cluster_data_fname,
                          downloadable = TRUE,
                          cached_plot = reactive(NULL))

  # filename generator for marker plot data
  marker_data_fname <- function(type, selected_group = NULL) {
    return(function() {

      if (type == 'cluster') {
        fname <- paste0(scForm$dataset_name(),
                        '_', scForm$selected_gene_cluster(),
                        '_marker_plot_data_', Sys.Date(), '.csv')

      } else if (type =='sample') {
        fname <- paste0(scForm$dataset_name(),
                        '_', selected_group,
                        '_', scForm$selected_gene_sample(),
                        '_marker_plot_data_', Sys.Date(), '.csv')
      }

      return(fname)
    })
  }

  # showing cluster comparison
  scMarkerCluster <- callModule(scMarkerPlot, 'marker_plot_cluster',
                                scseq = scForm$scseq,
                                selected_gene = scForm$selected_gene_cluster,
                                selected_group = 'all',
                                fname_fun = marker_data_fname('cluster'),
                                cached_plot = reactive(NULL))

  callModule(scBioGpsPlot, 'biogps_plot',
             selected_gene = scForm$selected_gene_cluster)


  callModule(scRidgePlot, 'ridge_plot',
             selected_gene = scForm$selected_gene_cluster,
             selected_cluster = scForm$selected_cluster,
             scseq = scForm$scseq)

  # showing sample comparison
  scMarkerSample <- callModule(scMarkerPlot, 'marker_plot_test',
                               scseq = scForm$scseq,
                               selected_gene = scForm$selected_gene_sample,
                               selected_group = 'test',
                               fname_fun = marker_data_fname('sample', 'test'),
                               cached_plot = reactive(NULL))

  callModule(scMarkerPlot, 'marker_plot_ctrl',
             scseq = scForm$scseq,
             selected_gene = scForm$selected_gene_sample,
             selected_group = 'ctrl',
             fname_fun = marker_data_fname('sample', 'ctrl'),
             cached_plot = scMarkerSample$plot)

  # showing labels comparison

  label_plot1 <- reactive({
    anal <- scForm$label_anals()[1]
    req(anal)

    plot <- scCluster$plot()
    scseq <- scForm$scseq()
    annot <- readRDS(scseq_part_path(sc_dir, anal, 'annot'))

    get_label_plot(anal, scseq, annot, plot)
  })

  label_plot2 <- reactive({
    anal <- scForm$label_anals()[2]
    req(anal)

    plot <- scCluster$plot()
    scseq <- scForm$scseq()
    annot <- readRDS(scseq_part_path(sc_dir, anal, 'annot'))

    get_label_plot(anal, scseq, annot, plot)
  })

  callModule(scClusterPlot, 'label_plot1',
             cached_plot = label_plot1)

  callModule(scClusterPlot, 'label_plot2',
             cached_plot = label_plot2)

  observe({
    toggle(id = "sample_comparison_row",  condition = scForm$comparison_type() == 'samples')
    toggle(id = "cluster_comparison_row", condition = scForm$comparison_type() == 'clusters')
    toggle(id = "label_comparison_row", condition = scForm$comparison_type() == 'labels')
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

  # TODO: implement loading single cell dataset

  # updates if new integrated dataset
  new_dataset <- reactive(scIntegration())


  # the dataset and options
  scDataset <- callModule(scSelectedDataset, 'dataset',
                          sc_dir = sc_dir,
                          new_dataset = new_dataset,
                          indices_dir = indices_dir)

  observe({
    toggle('form_container', condition = scDataset$dataset_exists())
  })

  #TODO move label transfer and integration into dataset

  # label transfer between datasets
  scLabelTransfer <- callModule(labelTransferForm, 'transfer',
                                sc_dir = sc_dir,
                                datasets = scDataset$datasets,
                                show_label_transfer = scDataset$show_label_transfer,
                                dataset_name = scDataset$dataset_name,
                                scseq = scDataset$scseq)

  # dataset integration
  scIntegration <- callModule(integrationForm, 'integration',
                              sc_dir = sc_dir,
                              datasets = scDataset$datasets,
                              show_integration = scDataset$show_integration)


  # ---


  # comparison type
  comparisonType <- callModule(comparisonType, 'comparison',
                               scseq = scDataset$scseq,
                               is.integrated = scDataset$is.integrated)


  # the selected cluster/gene for cluster comparison
  dataset_dir <- reactive(file.path(sc_dir, scDataset$dataset_name()))

  scClusterComparison <- callModule(clusterComparison, 'cluster',
                                    dataset_dir = dataset_dir,
                                    scseq = scDataset$scseq,
                                    markers = scDataset$markers,
                                    annot_path = scDataset$annot_path,
                                    ref_preds = scLabelTransfer)

  scClusterGene <- callModule(selectedGene, 'gene_clusters',
                              dataset_name = scDataset$dataset_name,
                              scseq = scDataset$scseq,
                              selected_markers = scClusterComparison$selected_markers,
                              cluster_markers = NULL,
                              selected_cluster = scClusterComparison$selected_cluster,
                              annot_path = scDataset$annot_path,
                              comparison_type = comparisonType)

  # the selected clusters/gene for sample comparison

  scSampleComparison <- callModule(scSampleComparison, 'sample',
                                   dataset_dir = dataset_dir)

  scSampleGene <- callModule(selectedGene, 'gene_samples',
                             dataset_name = scDataset$dataset_name,
                             scseq = scDataset$scseq,
                             selected_markers = scSampleComparison$top_table,
                             cluster_markers = scSampleComparison$cluster_markers,
                             selected_cluster = scSampleComparison$clusters,
                             annot_path = scDataset$annot_path,
                             comparison_type = comparisonType)


  # the original labels for integrated datasets
  integrationAnnotAnals <- callModule(selectedAnnot, 'annot',
                                      scseq = scDataset$scseq,
                                      is.integrated = scDataset$is.integrated,
                                      sc_dir = sc_dir)



  # update scseq with annotation changes
  scseq <- reactive({
    scseq <- scDataset$scseq()
    annot <- scClusterComparison$annot()
    req(annot, scseq)
    levels(scseq$cluster) <- annot
    return(scseq)
  })

  # show the toggle if dataset is integrated
  observe({
    toggle(id = "comparison_toggle_container",  condition = scDataset$is.integrated())
  })



  # show appropriate inputs based on comparison type
  observe({
    toggle(id = "label_comparison_inputs",  condition = comparisonType() == 'labels')
    toggle(id = "sample_comparison_inputs",  condition = comparisonType() == 'samples')
    toggle(id = "cluster_comparison_inputs", condition = comparisonType() == 'clusters')
  })

  return(list(
    scseq = scseq,
    selected_cluster = scClusterComparison$selected_cluster,
    selected_gene_cluster = scClusterGene$selected_gene,
    selected_gene_sample = scSampleGene$selected_gene,
    label_anals = integrationAnnotAnals,
    comparison_type = comparisonType,
    dataset_name = scDataset$dataset_name,
    show_ridge = scClusterGene$show_ridge

  ))
}

#' Logic for selected dataset part of scForm
#' @export
#' @keywords internal
scSelectedDataset <- function(input, output, session, sc_dir, new_dataset, indices_dir) {

  # get directory with fastqs
  roots <- c('single-cell' = sc_dir)
  shinyFiles::shinyDirChoose(input, "new_dataset_dir", roots = roots)

  dataset_name <- reactive({
    req(input$selected_dataset)
    req(!is.create())
    input$selected_dataset
  })

  dataset_exists <- reactive(isTruthy(input$selected_dataset) & !is.create())

  # get's used for saving annotation to disc
  annot_path <- reactive({
    scseq_part_path(sc_dir, dataset_name(), 'annot')
  })

  # load annotation for clusters
  annot <- reactive(readRDS(annot_path()))

  # load scseq
  scseq <- reactive({
    scseq_path <- scseq_part_path(sc_dir, dataset_name(), 'scseq')
    scseq <- readRDS(scseq_path)
    return(scseq)
  })

  is.integrated <- reactive({
    dataset_name <- dataset_name()
    req(dataset_name)
    integrated <- readRDS(file.path(sc_dir, 'integrated.rds'))
    return(dataset_name %in% integrated)
  })



  # load markers and name using annot
  markers <- reactive({
    annot <- annot()
    req(annot)

    markers_path <- scseq_part_path(sc_dir, dataset_name(), 'markers')
    markers <- readRDS(markers_path)
    names(markers) <- annot
    return(markers)
  })

  # available single-cell datasets
  datasets <- reactive({
    # reactive to new single cell datasets
    new_dataset()
    get_sc_dataset_choices(sc_dir)
  })


  # are we creating a new dataset?
  is.create <- reactive({
    dataset_name <- input$selected_dataset
    datasets <- datasets()

    !dataset_name %in% unlist(datasets)
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

  # run single-cell quantification
  observeEvent(input$confirm_quant, {
    removeModal()
    disable('selected_dataset')

    # standardize cellranger files
    fastq_dir <- new_dataset_dir()
    dataset_name <- input$selected_dataset
    is.cellranger <- check_is_cellranger(fastq_dir)
    if (is.cellranger) standardize_cellranger(fastq_dir)

    # Create a Progress object
    progress <- Progress$new(session, min=0, max = 8)
    on.exit({progress$close(); enable('selected_dataset')})

    progress$set(message = "Quantifying files", value = 1)
    if (!is.cellranger) run_kallisto_scseq(indices_dir, fastq_dir)

    # TODO: figure out whitelist vs kneelist
    # previously subseted to kneelist in load_scseq then to whitelist just after

    progress$set(message = "Loading and QC", value = 2)
    type <- ifelse(is.cellranger, 'cellranger', 'kallisto')
    scseq <- load_scseq(fastq_dir, project = dataset_name, type = type)
    scseq <- scseq[, scseq$whitelist]
    gc()

    progress$set(message = "Normalizing", value = 3)
    scseq <- normalize_scseq(scseq)
    gc()

    progress$set(message = "Clustering", value = 4)
    scseq <- add_hvgs(scseq)
    scseq <- add_scseq_clusters(scseq)
    gc()

    progress$set(message = "Reducing dimensions", value = 5)
    scseq <- run_tsne(scseq)
    gc()

    progress$set(message = "Getting markers", value = 6)
    wilcox_tests <- pairwise_wilcox(scseq)
    markers <- get_scseq_markers(wilcox_tests)

    # top markers for SingleR
    top_markers <- scran::getTopMarkers(wilcox_tests$statistics, wilcox_tests$pairs)

    progress$set(message = "Saving", value = 7)
    anal <- list(scseq = scseq, markers = markers, annot = names(markers), top_markers = top_markers)
    save_scseq_data(anal, dataset_name, sc_dir)

    progress$set(value = 9)
  })

  # modal to confirm adding single-cell dataset
  quantModal <- function() {

    UI <- withTags(dl(dd('This will take a while.')))

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


  # update if options change
  observe({
    updateSelectizeInput(session, 'selected_dataset', choices = datasets())
  })

  # show/hide integration/label-transfer forms
  show_integration <- reactive(input$show_integration %% 2 != 0)
  show_label_transfer <- reactive(input$show_label_transfer %% 2 != 0)

  observe(toggleClass(id = "show_label_transfer", 'btn-primary', condition = show_label_transfer()))
  observe(toggleClass(id = "show_integration", 'btn-primary', condition = show_integration()))

  # hide integration/label-transfer buttons no dataset
  observe({
    toggleSelectizeButtons('selected_dataset', c('show_integration', 'show_label_transfer'), dataset_exists())
  })


  # return anal and options to app
  return(list(
    dataset_name = dataset_name,
    scseq = scseq,
    markers = markers,
    annot = annot,
    annot_path = annot_path,
    datasets = datasets,
    show_integration = show_integration,
    show_label_transfer = show_label_transfer,
    is.integrated = is.integrated,
    dataset_exists = dataset_exists
  ))
}

get_sc_dataset_choices <- function(sc_dir) {

  # make sure integrated rds exists
  int_path <- file.path(sc_dir, 'integrated.rds')
  if (!file.exists(int_path)) saveRDS(NULL, int_path)

  # use saved datasets as options
  integrated <- readRDS(file.path(sc_dir, 'integrated.rds'))
  individual <- setdiff(list.files(sc_dir), c(integrated, 'integrated.rds'))

  # exclude individual without scseq (e.g. folder with fastq.gz files only)
  # unlist for case when no individual scseqs
  has.scseq <- sapply(individual, function(ind) any(list.files(file.path(sc_dir, ind)) == 'scseq.rds'))
  individual <- individual[unlist(has.scseq)]

  # must be a list if length one for option groups to work
  if (length(integrated) == 1) integrated <- list(integrated)
  if (length(individual) == 1) individual <- list(individual)

  list(Individual = c('', individual), Integrated = integrated)
}




#' Logic for label transfer between datasets
#' @export
#' @keywords internal
labelTransferForm <- function(input, output, session, sc_dir, datasets, show_label_transfer, dataset_name, scseq) {
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

    # load previously saved reference preds
    preds_path <- scseq_part_path(sc_dir, query_name, 'preds')
    if (file.exists(preds_path)) readRDS(preds_path) else list()
  })

  # update annotation transfer choices
  observe({
    preds <- preds()


    datasets <- datasets()
    dataset_name <- dataset_name()
    req(preds, datasets)

    choices <- get_label_transfer_choices(datasets, dataset_name, preds)
    updateSelectizeInput(session, 'ref_name', choices = choices, server = TRUE, selected = isolate(new_preds()), options = list(render = I('{option: transferLabelOption}')))
  })

  # submit annotation transfer
  observe({

    query_name <- dataset_name()
    ref_name <- input$ref_name
    preds <- preds()

    req(query_name, ref_name, preds)
    req(!ref_name %in% names(preds))
    req(show_label_transfer())

    toggleAll(label_transfer_inputs)

    # Create a Progress object
    progress <- Progress$new()
    progress$set(message = "Transfering labels", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())

    # Create a callback function to update progress.
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 2
      }
      progress$set(value = value, detail = detail)
    }
    n = 3

    # get arguments for SingleR
    query <- scseq()
    senv <- loadNamespace('SingleR')

    if (ref_name %in% ls(senv)) {
      ref <- get(ref_name, envir = senv)()
      labels <- ref$label.main
      genes <- 'de'

    } else {
      markers_path <- scseq_part_path(sc_dir, ref_name, 'top_markers')
      genes <- readRDS(markers_path)
      ref <- load_saved_scseq(ref_name, sc_dir)

      # need until SingleR #77 fixed
      common <- intersect(row.names(ref), row.names(query))
      genes <- lapply(genes, function(x) {lapply(x, function(y) intersect(y, common))})

      labels <- as.character(ref$cluster)
    }

    updateProgress(1/n)

    # take best label for each cluster
    pred <- SingleR::SingleR(test = query, ref = ref, labels = labels, genes = genes)
    tab <- table(assigned = pred$pruned.labels, cluster = query$cluster)

    pred <- row.names(tab)[apply(tab, 2, which.max)]

    preds_path <- scseq_part_path(sc_dir, query_name, 'preds')
    preds[[ref_name]] <- pred
    saveRDS(preds, preds_path)

    new_preds(ref_name)
    ref_preds(preds[[ref_name]])

    toggleAll(label_transfer_inputs)
  })

  # disable submit label transfer when already have preds
  observe({
    ref_name <- input$ref_name
    shinyjs::toggleState('submit_transfer', condition = is.null(ref_preds()))
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
    anal_name <- dataset_name()

    # show saved annot if nothing selected or label transfer not open
    if (is.null(ref_preds) | !show_label_transfer()) {
      annot_path <- scseq_part_path(sc_dir, anal_name, 'annot')
      annot <- readRDS(annot_path)

    } else {
      annot <- get_pred_annot(ref_preds, ref_name, anal_name, sc_dir)
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
    anal_name <- dataset_name()

    req(anal_name)

    showModal(transferModal())
  })

  observeEvent(input$confirm_overwrite, {
    removeModal()
    ref_name <- input$ref_name
    ref_preds <- ref_preds()
    anal_name <- dataset_name()

    req(anal_name)

    pred_annot <- get_pred_annot(ref_preds, ref_name, anal_name, sc_dir)
    annot_path <- scseq_part_path(sc_dir, anal_name, 'annot')
    saveRDS(pred_annot, annot_path)

    new_annot(pred_annot)
  })



  return(pred_annot)
}


#' Logic for integration form toggled by showIntegration
#' @export
#' @keywords internal
integrationForm <- function(input, output, session, sc_dir, datasets, show_integration) {

  integration_inputs <- c('ctrl_integration', 'integration_name', 'submit_integration', 'test_integration', 'exclude_clusters')


  integration_name <- reactive(input$integration_name)
  integration_options <- reactive(datasets()$Individual)

  ctrl <- reactiveVal()
  test <- reactiveVal()
  new_anal <- reactiveVal()

  # show/hide integration forms
  observe({
    toggle(id = "integration-form", anim = TRUE, condition = show_integration())
  })

  observe(ctrl(input$ctrl_integration))
  observe(test(input$test_integration))





  # update test dataset choices
  observe({
    options <- integration_options()
    updateSelectizeInput(session, 'test_integration', choices = options[!options %in% ctrl()], selected = isolate(test()))
  })

  # update control dataset choices
  observe({
    options <- integration_options()
    updateSelectizeInput(session, 'ctrl_integration', choices = options[!options %in% test()], selected = isolate(ctrl()))
  })

  # update exclude clusters
  observe({
    anal_names <- c(test(), ctrl())
    anal_colors <- get_palette(anal_names)
    exclude_choices <- get_exclude_choices(anal_names, sc_dir, anal_colors)
    updateSelectizeInput(session, 'exclude_clusters',
                         choices = exclude_choices,
                         options = list(render = I('{option: excludeOptions, item: excludeOptions}')),
                         server = TRUE)
  })

  # run integration
  observeEvent(input$submit_integration, {

    test_anals <- test()
    ctrl_anals <- ctrl()
    exclude_clusters <- input$exclude_clusters
    anal_name <- input$integration_name
    datasets <- datasets()

    error_msg <- validate_integration(test_anals, ctrl_anals, anal_name, datasets)

    if (is.null(error_msg)) {
      # clear error and disable button
      removeClass('validate', 'has-error')
      toggleAll(integration_inputs)

      # Create a Progress object
      progress <- Progress$new()
      progress$set(message = "Integrating datasets", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())

      # Create a callback function to update progress.
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = detail)
      }


      # run integration
      integrate_saved_scseqs(sc_dir,
                             test = test_anals,
                             ctrl = ctrl_anals,
                             exclude_clusters = exclude_clusters,
                             anal_name = anal_name,
                             updateProgress = updateProgress)


      # re-enable, clear inputs, and trigger update of available anals
      ctrl(NULL)
      test(NULL)
      new_anal(anal_name)
      updateTextInput(session, 'integration_name', value = '')
      toggleAll(integration_inputs)

    } else {
      # show error message
      html('error_msg', html = error_msg)
      addClass('validate', class = 'has-error')
    }

  })

  return(new_anal)
}

#' Logic for selected anals in integrated dataset to show original annotation plots for.
#' @export
#' @keywords internal
selectedAnnot <- function(input, output, session, scseq, is.integrated, sc_dir) {


  orig_anals <- reactive({
    req(is.integrated())
    scseq <- scseq()
    return(unique(scseq$batch))
  })

  observe({
    updateSelectizeInput(session, 'integration_anals', choices = orig_anals())
  })

  anals <- reactive(input$integration_anals)
  return(anals)
}

#' Logic for comparison type toggle for integrated analyses
#' @export
#' @keywords internal
comparisonType <- function(input, output, session, scseq, is.integrated) {

  # always show clusters if not integrated
  observe({
    if(!is.integrated())
      updateRadioGroupButtons(session, 'comparison_type', selected = 'clusters')
  })

  return(reactive(input$comparison_type))
}

#' Logic for cluster comparison input
#' @export
#' @keywords internal
clusterComparison <- function(input, output, session, dataset_dir, scseq, markers, annot_path, ref_preds) {


  contrast_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))
  selected_cluster <- reactiveVal()
  show_contrasts <- reactive({ input$show_contrasts %% 2 != 0 })
  con_markers <- reactiveVal(list())

  # things that return for plotting
  annot <- reactiveVal(NULL)
  selected_markers <- reactiveVal(NULL)

  show_rename <- reactive({
    (input$rename_cluster + input$show_rename) %% 2 != 0
  })

  test_cluster <- reactive({
    test_cluster <- input$selected_cluster
    req(test_cluster)
    gsub(' vs .+?$', '', test_cluster)
  })

  # update data.frame for cluster/contrast choices
  choices <- reactive({
    clusters <- annot()
    req(clusters)

    if (show_contrasts()) {
      test <- isolate(test_cluster())
      choices <- get_contrast_choices(clusters, test)

    } else {
      choices <- get_cluster_choices(clusters, dataset_dir())
    }

    return(choices)
  })


  all_markers <- reactive({
    con_markers <- con_markers()
    markers <- markers()
    names(markers) <- seq(1, along.with = markers)
    return(c(markers, con_markers))
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

  # reset if switch analysis
  observeEvent(dataset_dir(), {
    con_markers(list())
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
    toggleState('show_rename', condition = !show_contrasts())
    toggleClass(id = "show_contrasts", 'btn-primary', condition = show_contrasts())
  })


  # modify/save annot if rename a cluster
  observeEvent(input$rename_cluster, {

    req(input$new_cluster_name)

    # update reactive annotation
    choices <- choices()
    sel_clust <- selected_cluster()
    sel_idx <- which(choices$value == sel_clust)

    mod_annot <- annot()
    mod_annot[sel_idx] <- input$new_cluster_name
    mod_annot <- make.unique(mod_annot, '_')

    # save on disc
    saveRDS(mod_annot, annot_path())

    # update annot and set selected cluster to new name
    annot(mod_annot)
  })


  # update UI for contrast/cluster choices
  observeEvent(choices(), {

    updateSelectizeInput(session, 'selected_cluster',
                         choices = choices(),
                         selected = selected_cluster(),
                         options = contrast_options, server = TRUE)
  })


  # update ui for renaming a cluster
  observe({
    choices <- choices()
    name <- choices[choices$value == input$selected_cluster, 'name']
    if (!show_rename())
      updateTextInput(session, 'new_cluster_name', value = '', placeholder = paste('Type new name for', name, '...'))
  })


  # get cluster if don't have (for comparing specific cluster)
  observeEvent(input$selected_cluster, {
    sel <- input$selected_cluster
    req(sel)

    if (!sel %in% names(all_markers())) {
      con <- strsplit(sel, ' vs ')[[1]]
      con_markers <- con_markers()

      # returns both directions
      tests <- pairwise_wilcox(scseq(), restrict = con)
      markers <- get_scseq_markers(tests)
      names(markers) <- paste(names(markers), 'vs', rev(names(markers)))

      con_markers <- c(con_markers, markers)
      con_markers(con_markers)
    }
  })


  observe({
    sel <- selected_cluster()
    req(sel)
    selected_markers <- all_markers()[[sel]]
    selected_markers(selected_markers)
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
selectedGene <- function(input, output, session, dataset_name, scseq, selected_markers, cluster_markers, selected_cluster, annot_path, comparison_type) {

  selected_gene <- reactiveVal(NULL)
  gene_options <- list(render = I('{option: geneChoice, item: geneChoice}'))

  exclude_ambient <- reactive({
    if (is.null(input$exclude_ambient)) return(FALSE)
    input$exclude_ambient %% 2 != 0
  })

  # toggle for excluding ambient
  observe({
    toggleClass('exclude_ambient', class = 'btn-primary', condition = exclude_ambient())
  })

  # toggle for ridgeline
  show_ridge <- reactive(input$show_ridge %% 2 != 0)
  observe(toggleClass(id = "show_ridge", 'btn-primary', condition = show_ridge()))

  filtered_markers <- reactive({

    ambient <- negative <- NULL
    markers <- selected_markers()
    if (is.null(markers)) return(NULL)

    if (exclude_ambient())
      ambient <- get_ambient(scseq(), markers = markers, cluster_markers = cluster_markers())

    if (comparison_type() == 'samples')
      negative <- row.names(markers)[markers$logFC < 0]

    supress <- unique(c(ambient, negative))
    markers <- supress.genes(markers, supress)
    return(markers)
  })

  # update marker genes based on cluster selection
  gene_choices <- reactive({
    scseq <- scseq()
    markers <- filtered_markers()
    selected_cluster <- selected_cluster()
    comparison_type <- isolate(comparison_type())

    # will error if labels
    req(comparison_type %in% c('samples', 'clusters'))
    if (is.null(markers) || is.null(selected_cluster)) return(NULL)

    get_gene_choices(scseq, markers)
  })

  # click genecards
  observeEvent(input$genecards, {
    gene_link <- paste0('https://www.genecards.org/cgi-bin/carddisp.pl?gene=', input$selected_gene)
    runjs(paste0("window.open('", gene_link, "')"))
  })


  # click download
  output$download <- downloadHandler(
    filename = function() {

      annot <- readRDS(annot_path())
      selected_idx <- as.numeric(selected_cluster()) + 1

      sc_dl_filename(cluster = annot[selected_idx],
                     anal = dataset_name(),
                     comparison_type = comparison_type())

    },
    content = function(con) {
      write.csv(filtered_markers(), con)
    }
  )





  # reset selected gene if analysis changes
  observeEvent(dataset_name(), {
    selected_gene(NULL)
  })

  # set return value based on input
  observeEvent(input$selected_gene, {
    selected_gene(input$selected_gene)
  })


  # update choices
  observe({
    updateSelectizeInput(session, 'selected_gene',
                         choices = gene_choices(), selected = NULL,
                         server = TRUE,
                         options = gene_options)
  })

  return(list(
    selected_gene = selected_gene,
    show_ridge = show_ridge
  ))

}


#' Logic for cluster plots
#' @export
#' @keywords internal
scClusterPlot <- function(input, output, session, scseq, cached_plot, fname_fun = function(){}, downloadable = FALSE) {

  plot <- reactive({
    cached_plot <- cached_plot()
    if (!is.null(cached_plot)) return(cached_plot)
    plot_tsne_cluster(scseq())

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
    toggleClass('download_container', class = 'visible-plot', condition = downloadable && isTruthy(plot_fun()))
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
  })
}

#' Logic for marker gene plots
#' @export
#' @keywords internal
scMarkerPlot <- function(input, output, session, scseq, selected_gene, selected_group, cached_plot, fname_fun = function(){}, downloadable = TRUE) {


  plot <- reactive({
    # cached plot if showing test samples
    cached <- cached_plot()
    if (!is.null(cached)) return(cached)

    req(selected_gene())
    plot_tsne_gene(scseq(), selected_gene())
  })

  ploted_plot <- reactive({
    pl <- plot()
    req(pl)

    if (selected_group != 'all')
      pl <- format_sample_gene_plot(pl, selected_group, selected_gene(), scseq())

    return(pl)
  })


  data_fun <- function() {
    plot <- ploted_plot()
    plot$data
  }

  callModule(downloadablePlot,
             "marker_plot",
             plot_fun = ploted_plot,
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
scBioGpsPlot <- function(input, output, session, selected_gene) {
  # plot BioGPS data
  output$biogps_plot <- renderPlot({
    plot_biogps(selected_gene())
  })
}

scRidgePlot <- function(input, output, session, selected_gene, selected_cluster, scseq) {
  # plot BioGPS data


  output$ridge_plot <- renderPlot({
    req(selected_gene())
    plot_ridge(selected_gene(), selected_cluster(), scseq())
  })
}

plot_ridge <- function(gene, selected_cluster, scseq) {

  cluster <- scseq$cluster
  logcounts <- SingleCellExperiment::logcounts(scseq)[gene, ]
  cols <- get_palette(levels(cluster))
  clus <- levels(cluster)[as.numeric(selected_cluster)]
  color <- cols[as.numeric(selected_cluster)]


  mean <- tapply(logcounts, cluster, mean)
  levs <- levels(cluster)[order(mean)]
  df <- data.frame(logcounts, cluster = factor(cluster, levels = levs))

  df$group <- 'out'
  df$group[cluster == clus] <- 'in'

  ggplot2::ggplot(df, ggplot2::aes(x = logcounts, y = cluster, fill = group, alpha = group, color = group)) +
    ggplot2::scale_fill_manual(values = c(color, 'grey')) +
    ggplot2::scale_color_manual(values = c('black', 'gray')) +
    ggplot2::scale_alpha_manual(values = c(0.95, 0.25)) +
    ggridges::geom_density_ridges(scale = 2.5, rel_min_height = 0.001) +
    ggridges::theme_ridges(center_axis_labels = TRUE) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_dimgray() +
    ggplot2::xlab('Expression') +
    ggplot2::theme(legend.position = 'none',
                   axis.text.y = ggplot2::element_text(color = '#333333', size = 14),
                   axis.title.x = ggplot2::element_text(color = '#333333', size = 14),
                   axis.title.y = ggplot2::element_blank()
    )

}


#' Logic for single cell cluster analyses for Single Cell, Drugs, and Pathways tabs
#' @export
#' @keywords internal
scSampleComparison <- function(input, output, session, dataset_dir, is_sc = function()TRUE) {
  contrast_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))
  input_ids <- c('run_comparison', 'selected_clusters')

  # single cell dataset
  scseq <- reactive(readRDS(file.path(dataset_dir(), 'scseq.rds')))

  # TODO: update if annotation change from Single Cell tab
  annot <- reactive(readRDS(file.path(dataset_dir(), 'annot.rds')))

  # update cluster choices in UI
  cluster_choices <- reactive(get_cluster_choices(annot(), dataset_dir(), TRUE))

  observe({
    req(is_sc())
    updateSelectizeInput(session, 'selected_clusters',
                         choices = cluster_choices(),
                         options = contrast_options, server = TRUE)
  })


  # path to lmfit, cluster markers, and drug query results
  clusters_str <- reactive(collapse_sorted(input$selected_clusters))

  lmfit_path <- reactive(file.path(dataset_dir(), paste0('lm_fit_', clusters_str(), '.rds')))
  markers_path <- reactive(file.path(dataset_dir(), paste0('markers_', clusters_str(), '.rds')))
  drug_paths <- reactive(get_drug_paths(dataset_dir(), clusters_str()))

  # do we have lm_fit and drug query results?
  saved_lmfit <- reactive(file.exists(lmfit_path()))
  saved_drugs <- reactive(file.exists(drug_paths()$cmap))

  # unlike bulk two-group comparisons
  # need user to indicate when all clusters selected
  # also require cluster markers and fit result
  lm_fit <- reactiveVal()
  cluster_markers <- reactiveVal()

  observeEvent(input$run_comparison, {

    if (saved_lmfit()) {
      resl <- list(
        fit = readRDS(lmfit_path()),
        cluster_markers = readRDS(markers_path())
      )

    } else {
      toggleAll(input_ids)
      resl <- run_limma_scseq(
        scseq(),
        input$selected_clusters,
        dataset_dir()
      )
      toggleAll(input_ids)
    }

    lm_fit(resl$fit)
    cluster_markers(resl$cluster_markers)
  })


  # differential expression top table
  top_table <- reactive({
    req(lm_fit())
    tt <- get_top_table(lm_fit())
    tt[order(tt$P.Value), ]
  })

  # need ambient for pathway and drug queries
  ambient <- reactive({
    get_ambient(scseq(),
                markers = top_table(),
                cluster_markers = cluster_markers())
  })

  # drug query results
  drug_queries <- reactive({
    if (is.null(lm_fit())) {
      res <- NULL

    } else if (saved_drugs()) {
      paths <- drug_paths()
      res <- lapply(paths, readRDS)

    } else {
      toggleAll(input_ids)
      top_table <- top_table()
      res <- run_drug_queries(top_table, drug_paths(), session, ambient())
      toggleAll(input_ids)
    }
    return(res)
  })

  path_res <- reactive({})

  # reset lm_fit if selected clusters change
  # also disable runing lm_fit if nothing selected
  observe({
    lm_fit(NULL)
    cluster_markers(NULL)
    toggleState('run_comparison', condition = length(input$selected_clusters))
  })

  # name for  downloading
  dl_name <- reactive({
    sel <- input$selected_clusters
    req(sel)
    annot <- gsub(' ', '-', annot())
    inds <- as.numeric(sel) + 1
    paste0(annot[sort(inds)], collapse = '_')
  })

  return(list(
    name = dl_name,
    lm_fit = lm_fit,
    ambient = ambient,
    top_table = top_table,
    cluster_markers = cluster_markers,
    drug_queries = drug_queries,
    path_res = path_res,
    clusters = reactive(input$selected_clusters)
  ))
}

