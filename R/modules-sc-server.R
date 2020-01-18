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
                          downloadable = TRUE)

  # cluster comparison plots ---

  # filename generator for marker plot data
  cluster_fname <- function() {
    fname <- paste0(scForm$dataset_name(),
                    '_', scForm$cluster_gene(),
                    '_marker_plot_data_', Sys.Date(), '.csv')
    return(fname)
  }

  scMarkerCluster <- callModule(scMarkerPlot, 'marker_plot_cluster',
                                scseq = scForm$scseq,
                                selected_gene = scForm$cluster_gene,
                                fname_fun = cluster_fname)

  callModule(scBioGpsPlot, 'biogps_plot',
             selected_gene = scForm$cluster_gene)


  callModule(scRidgePlot, 'ridge_plot',
             selected_gene = scForm$cluster_gene,
             selected_cluster = scForm$selected_cluster,
             scseq = scForm$scseq)


  # sample comparison plots ---

  callModule(scSampleMarkerPlot, 'left',
             selected_gene = scForm$sample_gene,
             plot_fun = scForm$samples_pfun_left,
             height = scForm$samples_plot_height_left)

  callModule(scSampleMarkerPlot, 'right',
             selected_gene = scForm$sample_gene,
             plot_fun = scForm$samples_pfun_right,
             height = scForm$samples_plot_height_right)


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

  # updates if new integrated dataset
  new_dataset <- reactive(scIntegration())

  # the dataset and options
  scDataset <- callModule(scSelectedDataset, 'dataset',
                          sc_dir = sc_dir,
                          new_dataset = new_dataset,
                          indices_dir = indices_dir)

  observe(toggle('form_container', condition = scDataset$dataset_exists()))

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


  # comparison type
  comparisonType <- callModule(comparisonType, 'comparison',
                               scseq = scDataset$scseq,
                               is.integrated = scDataset$is.integrated)


  # the selected cluster/gene for cluster comparison
  dataset_dir <- reactive(file.path(sc_dir, scDataset$dataset_name()))

  scClusterComparison <- callModule(clusterComparison, 'cluster',
                                    dataset_dir = dataset_dir,
                                    scseq = scDataset$scseq,
                                    annot_path = scDataset$annot_path,
                                    ref_preds = scLabelTransfer)

  scClusterGene <- callModule(selectedGene, 'gene_clusters',
                              dataset_name = scDataset$dataset_name,
                              selected_markers = scClusterComparison$selected_markers,
                              selected_cluster = scClusterComparison$selected_cluster,
                              comparison_type = comparisonType)

  # the selected clusters/gene for sample comparison
  scSampleComparison <- callModule(scSampleComparison, 'sample',
                                   dataset_dir = dataset_dir,
                                   dataset_name = scDataset$dataset_name,
                                   input_scseq = scseq)

  scSampleGene <- callModule(selectedGene, 'gene_samples',
                             dataset_name = scDataset$dataset_name,
                             selected_markers = scSampleComparison$top_table,
                             selected_cluster = scSampleComparison$selected_clusters,
                             comparison_type = comparisonType,
                             ambient = scSampleComparison$ambient)




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
    cluster_gene = scClusterGene$selected_gene,
    show_ridge = scClusterGene$show_ridge,
    sample_gene = scSampleGene$selected_gene,
    samples_pfun_left = scSampleComparison$pfun_left,
    samples_pfun_right = scSampleComparison$pfun_right,
    samples_plot_height_right = scSampleComparison$height_right,
    samples_plot_height_left = scSampleComparison$height_left,
    sample_clusters = scSampleComparison$selected_clusters,
    comparison_type = comparisonType,
    dataset_name = scDataset$dataset_name
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
    tests <- pairwise_wilcox(scseq)
    markers <- get_scseq_markers(tests)

    # top markers for SingleR
    top_markers <- scran::getTopMarkers(tests$statistics, tests$pairs)

    progress$set(message = "Saving", value = 7)
    anal <- list(scseq = scseq, markers = markers, tests = tests, annot = names(markers), top_markers = top_markers)
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


scSampleMarkerPlot <- function(input, output, session, selected_gene, plot_fun, height) {

  output$plot <- renderPlot({
    gene <- selected_gene()
    req(gene)

    plot_fun()(gene)
  }, height = height)
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

    if (show_contrasts()) {
      test <- isolate(test_cluster())
      choices <- get_contrast_choices(clusters, test)

    } else {
      choices <- get_cluster_choices(clusters, dataset_dir())
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

  # reset if switch analysis
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


  # get/load markers if don't have
  observeEvent(input$selected_cluster, {
    sel <- input$selected_cluster
    req(sel)
    if (sel %in% names(markers())) return(NULL)

    markers <- markers()
    dataset_dir <- dataset_dir()

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
    req(sel)
    selected_markers(markers()[[sel]])
  })


  return(list(
    annot = annot,
    selected_markers = selected_markers,
    selected_cluster = selected_cluster
  ))
}

#' Get markers for one cluster against one other cluster
#'
#' @param con String identifying clusters to compare e.g. \code{'1vs2'}
#' @param dataset_dir Path to folder for dataset
#'
#' @return Named list with data.frames for each comparison direction
#' @export
#' @keywords internal
#'
get_contrast_markers <- function(con, dataset_dir) {
  con <- strsplit(con, '-vs-')[[1]]
  tests_dir <- file.path(dataset_dir, 'tests')
  pairs_path <- file.path(tests_dir, 'pairs.rds')
  pairs <- readRDS(pairs_path)

  keep <- apply(pairs, 1, function(row) all(con %in% row))
  keep <- which(keep)

  stat_paths <- file.path(tests_dir, paste0('statistics_pair', keep, '.rds'))
  tests <- list(pairs = pairs[keep, ],
                statistics = lapply(stat_paths, readRDS))

  # returns both directions
  con_markers <- get_scseq_markers(tests)
  names(con_markers) <- paste0(names(con_markers), '-vs-', rev(names(con_markers)))

  return(con_markers)
}


#' Logic for selected gene to show plots for
#' @export
#' @keywords internal
selectedGene <- function(input, output, session, dataset_name, selected_markers, selected_cluster, comparison_type, ambient = function()NULL) {

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

    markers <- selected_markers()
    if (is.null(markers)) return(NULL)

    ambient <- NULL
    if (exclude_ambient()) ambient <- ambient()

    markers <- supress.genes(markers, ambient)
    return(markers)
  })

  # update marker genes based on cluster selection
  gene_choices <- reactive({
    markers <- filtered_markers()
    selected_cluster <- selected_cluster()
    type <- isolate(comparison_type())

    # will error if labels
    req(type %in% c('samples', 'clusters'))
    if (is.null(markers) || is.null(selected_cluster)) return(NULL)

    get_gene_choices(markers)
  })

  # click genecards
  observeEvent(input$genecards, {
    gene_link <- paste0('https://www.genecards.org/cgi-bin/carddisp.pl?gene=', input$selected_gene)
    runjs(paste0("window.open('", gene_link, "')"))
  })


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
scClusterPlot <- function(input, output, session, scseq, fname_fun = function(){}, downloadable = FALSE) {

  plot <- reactive({
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
scMarkerPlot <- function(input, output, session, scseq, selected_gene, fname_fun = function(){}, downloadable = TRUE) {


  plot <- reactive({
    req(selected_gene())
    plot_tsne_gene(scseq(), selected_gene())
  })

  ploted_plot <- reactive({
    pl <- plot()
    req(pl)
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

  height <- reactive(length(levels(scseq()$cluster))*50)

  output$ridge_plot <- renderPlot({
    req(selected_gene())
    plot_ridge(selected_gene(), scseq(), selected_cluster())
  }, height = height)
}


#' Logic for single cell cluster analyses for Single Cell, Drugs, and Pathways tabs
#' @export
#' @keywords internal
scSampleComparison <- function(input, output, session, dataset_dir, dataset_name, input_scseq = function()NULL, is_sc = function()TRUE) {
  contrast_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))
  input_ids <- c('run_comparison', 'selected_clusters')

  # use input scseq/annot if available
  scseq <- reactive({
    scseq <- input_scseq()
    if (is.null(scseq)) scseq <- readRDS(file.path(dataset_dir(), 'scseq.rds'))
    return(scseq)
  })

  annot <- reactive({
    annot <- levels(input_scseq()$cluster)
    if (is.null(annot)) annot <- readRDS(file.path(dataset_dir(), 'annot.rds'))
    return(annot)
  })

  has_replicates <- reactive(readRDS(file.path(dataset_dir(), 'has_replicates.rds')))
  summed <- reactive(readRDS(file.path(dataset_dir(), 'summed.rds')))

  # pseudo bulk eset
  pbulk <- reactiveVal()
  pbulk_path <- reactive(file.path(dataset_dir(), 'pbulk_eset.rds'))

  # update cluster choices in UI
  cluster_choices <- reactive(get_cluster_choices(annot(), dataset_dir(), TRUE))

  observe({
    req(is_sc())
    updateSelectizeInput(session, 'selected_clusters',
                         choices = cluster_choices(),
                         options = contrast_options, server = TRUE)
  })


  # path to lmfit, cluster markers, drug query results, and goanna pathway results
  clusters_str <- reactive(collapse_sorted(input$selected_clusters))
  lmfit_path <- reactive(file.path(dataset_dir(), 'lm_fit_0svs.rds'))
  markers_path <- reactive(file.path(dataset_dir(), paste0('markers_', clusters_str(), '.rds')))
  drug_paths <- reactive(get_drug_paths(dataset_dir(), clusters_str()))
  goana_path <- reactive(file.path(dataset_dir(), paste0('goana_', clusters_str(), '.rds')))

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
      reps <- has_replicates()
      obj <- if (reps) summed() else scseq()
      resl <- list(
        fit = run_limma_scseq(obj, dataset_dir(), is_summed = reps),
        cluster_markers = get_cluster_markers(input$selected_clusters, dataset_dir())
      )
      toggleAll(input_ids)
    }

    if (has_replicates()) {
      pbulk <- readRDS(pbulk_path())
      levels(pbulk$cluster) <- annot()
      pbulk(pbulk)
    }

    lm_fit(resl$fit)
    cluster_markers(resl$cluster_markers)
  })


  # plot functions for left
  pfun_left_reps <- reactive(function(gene) plot_scseq_gene_medians(gene, pbulk(), top_tables()))
  pfun_left_noreps <- reactive(function(gene) plot_tsne_gene_sample(gene, scseq(), 'test'))
  pfun_left <- reactive(if (has_replicates()) pfun_left_reps() else pfun_left_noreps())

  # plot function for right
  pfun_right_reps <- reactive(function(gene) plot_ridge(gene, scseq(), input$selected_clusters, by.sample = TRUE))
  pfun_right_noreps  <- reactive(function(gene) plot_tsne_gene_sample(gene, scseq(), 'ctrl'))
  pfun_right <- reactive(if (has_replicates()) pfun_right_reps() else pfun_right_noreps())

  # plot heights
  height_left <- reactive(if(has_replicates()) length(annot())*50 else 453)
  height_right <- reactive(if(has_replicates()) length(unique(scseq()$batch))*76 else 453)


  # differential expression top tables for all clusters
  top_tables <- reactive({
    req(lm_fit())

    fit <- lm_fit()
    groups <-  c('test', 'ctrl')
    clusts <- cluster_choices()$value

    tests <- paste0('test_', clusts)
    ctrls <- paste0('ctrl_', clusts)
    have.groups <- tests %in% colnames(fit$mod) & ctrls %in% colnames(fit$mod)

    contrasts <- paste(tests, ctrls, sep = '-')[have.groups]
    tts <- get_top_tables(fit, contrasts = contrasts)

    return(tts)

  })

  # top table for selected cluster only
  top_table <- reactive({
    sel <- input$selected_clusters
    contrast <- paste0('test_', sel, '-ctrl_', sel)
    tts <- top_tables()
    req(contrast %in% names(tts))
    tts[[contrast]]
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


  # goana pathway result
  path_res <- reactive({
    goana_path <- goana_path()

    if (file.exists(goana_path)) {
      goana_res <- readRDS(goana_path)

    } else {
      lm_fit <- lm_fit()

      selected_clusters <- input$selected_clusters
      groups <- paste(c('test', 'ctrl'), selected_clusters, sep = '_')

      # loses sync when groups selected and change dataset
      req(all(groups %in% colnames(lm_fit$mod)))

      # exclude ambient
      ambient <- ambient()
      lm_fit <- within(lm_fit, fit <- fit[!row.names(fit) %in% ambient, ])

      contrast <- paste0(groups[1], '-', groups[2])
      ebfit <- fit_ebayes(lm_fit, contrast)
      goana_res <- limma::goana(ebfit, species = 'Hs', geneid = 'ENTREZID')

      saveRDS(goana_res, goana_path)
    }

    return(goana_res)
  })


  # reset lm_fit if selected clusters change
  # also disable runing lm_fit if nothing selected
  observe({
    lm_fit(NULL)
    cluster_markers(NULL)
    toggleState('run_comparison', condition = length(input$selected_clusters))
  })


  # enable download
  observe({
    toggleState('download', condition = isTruthy(lm_fit()))
  })

  # name for  downloading
  dl_fname <- reactive({
    date <- paste0(Sys.Date(), '.zip')

    # annotated clusters
    clusts <- input$selected_clusters
    req(clusts)
    annot <- gsub(' ', '-', annot())
    clusts <- as.numeric(clusts)
    clusts <- paste0(annot[sort(clusts)], collapse = '_')

    paste('single-cell', dataset_name(), clusts, date , sep='_')
  })

  data_fun <- function(file) {
    #go to a temp dir to avoid permission issues
    owd <- setwd(tempdir())
    on.exit(setwd(owd))

    tt_fname <- 'top_table.csv'
    go_fname <- 'goana.csv'

    tt <- top_table()
    tt$ambient <- row.names(tt) %in% ambient()
    write.csv(tt, tt_fname)
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
    lm_fit = lm_fit,
    ambient = ambient,
    top_table = top_table,
    cluster_markers = cluster_markers,
    drug_queries = drug_queries,
    path_res = path_res,
    selected_clusters = reactive(input$selected_clusters),
    pfun_left = pfun_left,
    pfun_right = pfun_right,
    height_left = height_left,
    height_right = height_right
  ))
}


#' Format gene plots for sample comparison for drugseqr app
#'
#' @param plot Returned by \code{\link{plot_umap_gene}}
#' @param group Level in \code{scseq$orig.ident} to show cells for. Either \code{'ctrl'} or \code{'test'}
#' @param scseq \code{Seurat} object.
#'
#' @return \code{plot} formatted for drugseqr app
#' @export
#' @keywords internal
plot_tsne_gene_sample <- function(gene, scseq, group = 'test') {

  plot <- plot_tsne_gene(scseq, gene)
  gene <- make.names(gene)

  # the min and max gene expression value
  lims <- range(plot$data[[gene]])

  # show selected group only
  sel.cells <- colnames(scseq)[scseq$orig.ident == group]
  plot$data <- plot$data[row.names(plot$data) %in% sel.cells, ]

  # add selected group as title
  plot <- plot + ggplot2::ggtitle(toupper(group)) +
    ggplot2::theme(plot.title = ggplot2::element_text(color = 'black'))

  # use the same scale for each sample so that comparable
  suppressMessages(plot <- plot + ggplot2::scale_color_continuous(low ="lightgray", high = "blue", limits = lims))

  # remove control plot labels and legend
  if (group == 'ctrl')
    plot <- plot + ggplot2::xlab('') + ggplot2::ylab('') +
    ggplot2::theme(legend.position = 'none')

  return(plot)
}


plot_ridge <- function(gene, scseq, selected_cluster, by.sample = FALSE) {


  # get color for selected cluster
  clus_levs <- levels(scseq$cluster)
  seli <- as.numeric(selected_cluster)
  sel <- clus_levs[seli]
  colors <- get_palette(clus_levs)
  color  <- colors[seli]

  # either highlight test group or selected cluster
  if (by.sample) {
    scseq <- scseq[, scseq$cluster %in% sel]
    y <- factor(scseq$batch)
    hl <- scseq$orig.ident == 'test'
    title <- paste('Expression by Sample for', sel)

  } else {
    y <- scseq$cluster
    hl <- y == sel
    title <- 'Expression by Cluster'
  }

  # order y-axis by mean logcounts
  x <- SingleCellExperiment::logcounts(scseq)[gene, ]
  m <- tapply(x, y, mean)
  y <- factor(y, levels = levels(y)[order(m)])
  df <- data.frame(x, hl, y)

  (pl <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = hl, alpha = hl, color = hl)) +
      ggplot2::scale_fill_manual(values = c('grey', color)) +
      ggplot2::scale_alpha_manual(values = c(0.25, 0.95)) +
      ggplot2::scale_color_manual(values = c('gray', 'black')) +
      ggridges::geom_density_ridges(scale = 3, rel_min_height = 0.001) +
      ggridges::theme_ridges(center_axis_labels = TRUE) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::coord_cartesian(clip = "off") +
      theme_dimgray() +
      ggplot2::xlab('') +
      ggplot2::ggtitle(title) +
      ggplot2::theme(legend.position = 'none',
                     plot.title = ggplot2::element_text(color = '#333333', size = 16, face = 'plain', margin = ggplot2::margin(b = 25) ),
                     axis.text.y = ggplot2::element_text(color = '#333333', size = 14),
                     axis.text.x = ggplot2::element_text(color = '#333333', size = 14),
                     axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank()
      ))

  ncells <- tapply(x, y, length)
  xlim <- ggplot2::ggplot_build(pl)$layout$panel_scales_x[[1]]$range$range[2]

  for (i in seq_along(ncells)) {
    (pl <- pl + ggplot2::annotation_custom(
      grob = grid::textGrob(label = paste(ncells[i], 'cells'), vjust = -0.3, just = 'right',
                            gp = grid::gpar(fontsize = 14, col = '#333333')),
      ymax = i,
      ymin = i,
      xmax = xlim,
      xmin = xlim
    ))
  }



  return(pl)
}


plot_scseq_gene_medians <- function(gene, pbulk, tts) {

  logcounts <- Biobase::assayDataElement(pbulk, 'vsd')

  group <- pbulk$orig.ident
  clust <- pbulk$cluster
  batch <- pbulk$batch

  # highlight clusters where this gene is significant
  pvals <- sapply(tts, function(tt) tt[gene, 'adj.P.Val'])
  have.clus <- gsub('^.+?_([0-9]+)$', '\\1', names(pvals))
  have.clus <- levels(clust)[as.numeric(have.clus)]
  miss.clus <- setdiff(levels(clust), have.clus)
  names(pvals) <- have.clus
  pvals[miss.clus] <- 1
  pvals <- pvals[levels(clust)]

  reds <- RColorBrewer::brewer.pal(9, 'YlOrRd')[c(3, 6, 9)]

  setup_color <- function(clust, group, sel, pvals) {

    res <- rep('white', length(clust))
    res[group == 'test'] <- '#87CEEB'
    res[group == 'test' & pvals[clust] < 0.15] <- reds[1]
    res[group == 'test' & pvals[clust] < 0.10] <- reds[2]
    res[group == 'test' & pvals[clust] < 0.05] <- reds[3]
    return(factor(res, levels = c(rev(reds), '#87CEEB', 'white')))
  }

  tb <- tibble::tibble(logcounts = logcounts[gene, ],
                       clust,
                       group,
                       batch,
                       pt.color = setup_color(clust, group, sel, pvals))

  tb <- tb %>%
    group_by(batch, clust) %>%
    add_tally() %>%
    summarise(meds = median(logcounts),
              n = unique(n),
              mads = stats::mad(logcounts),
              group = unique(group),
              pt.color = unique(pt.color))

  # order by significance
  tb$clust <- factor(tb$clust, levels(clust)[order(-pvals)])

  labels <- c('p<0.05', 'p<0.1', 'p<0.15', 'p>0.15', 'Control')
  labels <- labels[table(tb$pt.color) > 0]


  ggplot2::ggplot(tb, ggplot2::aes(x = meds, y = clust, fill = pt.color)) +
    ggplot2::scale_fill_identity('', guide = 'legend', labels = labels) +
    ggplot2::geom_point(shape = 21, color = '#333333', size = 4, alpha = 0.8) +
    ggplot2::xlab('') +
    ggplot2::ylab('') +
    ggplot2::ggtitle('Pseudobulk Expression by Cluster') +
    ggpubr::theme_pubr(legend = 'top')  +
    theme_dimgray() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 14, color = '#333333'),
                   plot.title = ggplot2::element_text(size = 16, color = '#333333', face = 'plain', margin = ggplot2::margin(b = 25)),
                   axis.title.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_line(size = 0),
                   legend.text = ggplot2::element_text(size = 16, color = '#333333'),
                   panel.grid.major.y = ggplot2::element_line(linetype = 'longdash', size = 0.3, color = 'black'))




}





