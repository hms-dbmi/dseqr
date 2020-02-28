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
                          selected_cluster = scForm$selected_cluster,
                          fname_fun = cluster_data_fname,
                          downloadable = TRUE)

  # cluster comparison plots ---

  # filename generator for marker plot data
  cluster_fname <- function() {
    fname <- paste0(scForm$dataset_name(),
                    '_', scForm$clusters_gene(),
                    '_marker_plot_data_', Sys.Date(), '.csv')
    return(fname)
  }

  scMarkerCluster <- callModule(scMarkerPlot, 'marker_plot_cluster',
                                scseq = scForm$scseq,
                                selected_gene = scForm$clusters_gene,
                                fname_fun = cluster_fname)

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

scLabelsPlot <- function(input, output, session, sc_dir, selected_cluster, scseq) {

  output$labels_plot <- plotly::renderPlotly({
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)

    cluster <- as.numeric(selected_cluster())
    cluster <- levels(scseq$cluster)[cluster]
    if (is.na(cluster)) return(NULL)

    plot_cluster_labels(scseq, cluster, sc_dir)
  })

}

plot_cluster_labels <- function(scseq, clust, sc_dir) {

  df <- tibble::as_tibble(scseq@colData) %>%
    dplyr::filter(cluster == clust) %>%
    dplyr::group_by(batch) %>%
    dplyr::mutate(nsample = n()) %>%
    dplyr::group_by(batch, orig.cluster) %>%
    dplyr::summarize(ncell = n(),
                     group = orig.ident[1],
                     nsample = nsample[1],
                     pcell = n() / nsample[1] * 100,
                     boc = paste(batch[1], orig.cluster[1], sep='_')) %>%
    dplyr::arrange(group, batch, pcell) %>%
    dplyr::select(-orig.cluster)

  # get mapping between original cluster index and name
  annots <- c()
  batches <- unique(df$batch)

  for (i in seq_along(batches)) {
    batch <- batches[i]
    annot_path <- file.path(sc_dir, batch, 'annot.rds')
    annot <- readRDS(annot_path)
    names(annot) <- paste(batch, seq_along(annot), sep='_')
    annots <- c(annots, annot)
  }

  df$oc <- annots[df$boc]
  df$customdata <- paste(df$ncell, 'of', df$nsample, 'cells')

  (pl <- plotly::plot_ly(data = df,
                         y = ~boc,
                         x = ~pcell,
                         color = ~batch,
                         customdata = ~customdata,
                         height = (nrow(df)*30) + 140,
                         text = ~nsample,
                         type = 'scatter',
                         mode = 'markers',
                         marker = list(size = 9, line = list(width = 1, color = '#333333')),
                         hoverlabel = list(bgcolor = '#000000', align = 'left'),
                         hovertemplate = paste0(
                           '<span style="color: crimson; font-weight: bold; text-align: left;">From Sample</span>: %{customdata}<br>',
                           '<extra></extra>')
  ) %>%
      plotly::layout(
        hovermode= 'closest',
        margin = list(t = 65, r = 20, l = 0, pad = 10),
        title = list(text = 'For Each Sample: Percent of Cells From Original Clusters', x = 0, y = .99, yanchor = 'top', font = list(size = 16, color = '#333333')),
        legend = list(traceorder = 'reversed', itemclick = FALSE, font = list(size = 14, color = '#333333')),
        xaxis = list(title = '', range = c(-4, 104), fixedrange=TRUE, side = 'top', zeroline = FALSE, tickfont = list(size = 14, color = '#333333')),
        yaxis = list(title = '', fixedrange=TRUE, tickvals = df$boc, ticktext = df$oc, tickmode = 2, gridwidth = 1, gridcolor = 'gray', ticklen = 10, tickcolor = 'white', tickfont = list(size = 14, color = '#333333'))
      ) %>%
      plotly::config(displayModeBar = 'hover',
                     displaylogo = FALSE,
                     doubleClick = 0,
                     showAxisDragHandles = FALSE,
                     showAxisRangeEntryBoxes = FALSE,
                     showTips = FALSE,
                     modeBarButtonsToRemove = c('zoom2d',
                                                'pan2d',
                                                'autoScale2d',
                                                'resetScale2d',
                                                'hoverClosestCartesian',
                                                'hoverCompareCartesian',
                                                'select2d',
                                                'lasso2d',
                                                'zoomIn2d',
                                                'zoomOut2d',
                                                'toggleSpikelines'),
                     toImageButtonOptions = list(format = "png")
      ))

  return(pl)
}


#' Logic for form on Single Cell Exploration page
#' @export
#' @keywords internal
scForm <- function(input, output, session, sc_dir, indices_dir) {

  # updates if new integrated dataset
  new_dataset <- reactiveVal()
  observe(new_dataset(scIntegration()))

  # the dataset and options
  scDataset <- callModule(scSelectedDataset, 'dataset',
                          sc_dir = sc_dir,
                          new_dataset = new_dataset,
                          indices_dir = indices_dir)

  observe(toggle('form_container', condition = scDataset$dataset_exists()))

  # update scseq with annotation changes
  scseq <- reactive({
    scseq <- scDataset$scseq()
    annot <- scClusterComparison$annot()
    if (!isTruthy(annot) | !isTruthy(scseq)) return(NULL)
    levels(scseq$cluster) <- annot
    return(scseq)
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
                              dataset_name = scDataset$dataset_name,
                              show_integration = scDataset$show_integration)


  # comparison type
  comparisonType <- callModule(comparisonType, 'comparison',
                               scseq = scseq,
                               is.integrated = scDataset$is.integrated)


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
                              dataset_name = scDataset$dataset_name,
                              is.integrated = scDataset$is.integrated,
                              selected_markers = scClusterComparison$selected_markers,
                              selected_cluster = scClusterComparison$selected_cluster,
                              qc_metrics = scClusterComparison$qc_metrics,
                              type = 'clusters')

  # the selected clusters/gene for sample comparison
  scSampleComparison <- callModule(scSampleComparison, 'sample',
                                   dataset_dir = dataset_dir,
                                   dataset_name = scDataset$dataset_name,
                                   is.integrated = scDataset$is.integrated,
                                   input_scseq = scseq,
                                   comparison_type = comparisonType,
                                   exclude_ambient = scSampleGene$exclude_ambient)

  scSampleGene <- callModule(selectedGene, 'gene_samples',
                             dataset_name = scDataset$dataset_name,
                             is.integrated = scDataset$is.integrated,
                             selected_markers = scSampleComparison$top_table,
                             selected_cluster = scSampleComparison$selected_cluster,
                             type = 'samples',
                             ambient = scSampleComparison$ambient)



  scLabelsComparison <- callModule(scLabelsComparison, 'labels',
                                   cluster_choices = scSampleComparison$cluster_choices)


  # show the toggle if dataset is integrated
  observe({
    toggle(id = "comparison_toggle_container",  condition = scDataset$is.integrated())
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
    show_ridge = scClusterGene$show_ridge,
    samples_pfun_left = scSampleComparison$pfun_left,
    samples_pfun_right = scSampleComparison$pfun_right,
    clusters_cluster = scClusterComparison$selected_cluster,
    samples_cluster = scSampleComparison$selected_cluster,
    labels_cluster = scLabelsComparison$selected_cluster,
    selected_cluster = selected_cluster,
    comparison_type = comparisonType,
    dataset_name = scDataset$dataset_name,
    species = scDataset$species
  ))
}

scLabelsComparison <- function(input, output, session, cluster_choices) {
  contrast_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))

  observe({
    updateSelectizeInput(session, 'selected_cluster',
                         choices = cluster_choices(),
                         options = contrast_options, server = TRUE)
  })

  return(list(
    selected_cluster = reactive(input$selected_cluster)
  ))
}

#' Logic for selected dataset part of scForm
#' @export
#' @keywords internal
scSelectedDataset <- function(input, output, session, sc_dir, new_dataset, indices_dir) {
  dataset_inputs <- c('selected_dataset', 'show_integration', 'show_label_transfer')
  options <- list(render = I('{option: scDatasetOptions, item: scDatasetItem}'))

  metric_choices <- c('low_lib_size',
                      'low_n_features',
                      'high_subsets_mito_percent',
                      'low_subsets_ribo_percent',
                      'high_doublet_score',
                      'high_outlyingness')

  # get directory with fastqs
  roots <- c('single-cell' = sc_dir)
  shinyFiles::shinyDirChoose(input, "new_dataset_dir", roots = roots)

  dataset_exists <- reactive(isTruthy(input$selected_dataset) & !is.create())

  dataset_name <- reactive({
    if (!dataset_exists()) return(NULL)
    input$selected_dataset
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
    toggleAll(dataset_inputs)

    dataset_dir <- file.path(sc_dir, dataset_name())
    scseq <- load_scseq(dataset_dir)

    toggleAll(dataset_inputs)
    return(scseq)
  })

  is.integrated <- reactive({
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


  # update last selected dataset on-file if changes
  prev_path <- file.path(sc_dir, 'prev_dataset.rds')

  observeEvent(input$selected_dataset, {
    sel <- input$selected_dataset
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

  # run single-cell quantification
  observeEvent(input$confirm_quant, {

    metrics <- input$qc_metrics
    # none, all, all and none: can't combine
    if (length(metrics) > 1 && !all(metrics %in% metric_choices)) return(NULL)

    if (!isTruthy(metrics)) metrics <- 'none'
    if (metrics == 'none') metrics <- NULL
    if (metrics == 'all') metrics <- metric_choices

    removeModal()
    toggleAll(dataset_inputs)

    # Create a Progress object
    progress <- Progress$new(session, min=0, max = 9)
    on.exit(progress$close())

    fastq_dir <- new_dataset_dir()
    dataset_name <- input$selected_dataset

    if (metrics == 'all and none') {
      load_raw_scseq(paste0(dataset_name, '_QC0'), fastq_dir, sc_dir, indices_dir, progress, metrics = NULL)
      load_raw_scseq(paste0(dataset_name, '_QC1'), fastq_dir, sc_dir, indices_dir, progress, metrics = metric_choices, founder = paste0(dataset_name, '_QC0'))
      saveRDS(paste0(dataset_name, '_QC1'), prev_path)

    } else {
      load_raw_scseq(dataset_name, fastq_dir, sc_dir, indices_dir, progress, metrics = metrics)
      saveRDS(dataset_name, prev_path)
    }

    toggleAll(dataset_inputs)
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
    selected <- readRDS.safe(prev_path)
    updateSelectizeInput(session, 'selected_dataset', choices = rbind(NA, datasets()), selected = selected, server = TRUE, options = options)
  })

  # show/hide integration/label-transfer forms
  show_integration <- reactive(input$show_integration %% 2 != 0)
  show_label_transfer <- reactive(input$show_label_transfer %% 2 != 0)

  observe(toggleClass(id = "show_label_transfer", 'btn-primary', condition = show_label_transfer()))
  observe(toggleClass(id = "show_integration", 'btn-primary', condition = show_integration()))

  # hide integration/label-transfer buttons no dataset
  observe({
    toggle('show_label_transfer-parent', condition = dataset_exists())
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
    dataset_exists = dataset_exists,
    species = species
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

  # founder for subsets as type
  ind_type <- sapply(individual,
                     function(ind) readRDS.safe(file.path(sc_dir, ind, 'founder.rds'), .nofile = 'Individual', .nullfile = 'Individual'), USE.NAMES = FALSE)
  sub <- ind_type != 'Individual'

  # exclude founder name from option label
  opt_label <- individual
  opt_label[sub] <- stringr::str_replace(opt_label[sub], paste0(ind_type[sub], '_'), '')

  choices <- data.frame(value = c(integrated, individual),
                        type = c(rep('Integrated', length(integrated)), ind_type),
                        itemLabel = stringr::str_trunc(c(integrated, individual), 35),
                        optionLabel = stringr::str_trunc(c(integrated, opt_label), 35),
                        stringsAsFactors = FALSE)

  return(choices)
}



scSampleMarkerPlot <- function(input, output, session, selected_gene, plot_fun) {

  res <- reactive({
    gene <- selected_gene()
    req(gene)
    suppressMessages(plot_fun()(gene))
  })

  height <- reactive(ifelse(is_plotly(), 1, res()$height))

  is_plotly <- output$is_plotly <- reactive('plotly' %in% class(res()))

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
    if (file.exists(preds_path)) readRDS(preds_path) else list()
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
      ref <- readRDS(scseq_part_path(sc_dir, ref_name, 'scseq'))

      # check if ref and query have the same founder
      rfound <- readRDS(scseq_part_path(sc_dir, ref_name, 'founder'))
      qfound <- readRDS(scseq_part_path(sc_dir, query_name, 'founder'))

      # use common cells to transfer labels if so
      cells <- intersect(colnames(ref), colnames(query))

      if (identical(qfound, rfound) && length(cells)) {

        ref_cluster <- ref[, cells]$cluster
        query_cluster <- query[, cells]$cluster
        tab <-  table(assigned = ref_cluster, cluster = query_cluster)

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

    preds_path <- scseq_part_path(sc_dir, query_name, 'preds')
    preds[[ref_name]] <- pred
    saveRDS(preds, preds_path)

    new_preds(ref_name)
    ref_preds(preds[[ref_name]])

    toggleAll(label_transfer_inputs)
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
integrationForm <- function(input, output, session, sc_dir, datasets, show_integration, dataset_name) {

  integration_inputs <- c('ctrl_integration',
                          'integration_name',
                          'submit_integration',
                          'test_integration',
                          'integration_type',
                          'exclude_clusters',
                          'click_up',
                          'click_dl',
                          'toggle_exclude')


  integration_name <- reactive(input$integration_name)

  # datasets() with server side selectize causes bug
  integration_choices <- reactive({
    ds <- datasets()
    ds <- ds[!ds$type %in% 'Integrated', ]
    sel <- dataset_name()

    if (isTruthy(sel)) {
      is.sel <- ds$value == sel
      type   <- ds$type[is.sel]
      is.sub <- type != 'Individual'

      # move within-founder datasets to top of choices
      if (is.sub) ds$type <- factor(ds$type, levels = unique(c(type, ds$type)))
    }

    choices <- ds %>%
      dplyr::group_by(type) %>%
      dplyr::summarise(values = list(value))

    names(choices$values) <- choices$type
    choices$values
  })

  ctrl <- reactiveVal()
  test <- reactiveVal()
  new_anal <- reactiveVal()
  selected_datasets <- reactive(c(test(), ctrl()))


  is_include <- reactive({ input$toggle_exclude %% 2 != 0 })
  is_subset <- reactive(length(test()) == 1 && is.null(ctrl()))
  is_integration <- reactive(length(test()) && length(ctrl()))
  allow_pairs <- reactive(length(selected_datasets()) > 2 & is_integration())

  # show/hide integration forms
  observe({
    toggle(id = "integration-form", anim = TRUE, condition = show_integration())
  })

  observe(ctrl(input$ctrl_integration))
  observe(test(input$test_integration))

  # update test dataset choices
  observe({
    choices <- integration_choices()
    choices <- lapply(choices, function(x) x[!x %in% ctrl()])
    updateSelectizeInput(session, 'test_integration', choices = choices, selected = isolate(test()))
  })

  # update control dataset choices
  observe({
    choices <- integration_choices()
    choices <- lapply(choices, function(x) x[!x %in% test()])
    updateSelectizeInput(session, 'ctrl_integration', choices = choices, selected = isolate(ctrl()))
  })



  # hide pairing if not enough datasets
  observe({
    toggle(id = "click_dl", condition = allow_pairs())
    toggle(id = "click_up", condition = allow_pairs())
  })

  # show cluster type choices if enough datasets
  observe(toggle(id = 'integration_type', condition = is_integration()))


  # show exclude/exclude and new dataset only if something selected
  observe(toggle(id = 'exclude-container', condition = selected_datasets()))
  observe(toggle(id = 'name-container', condition = selected_datasets()))

  # update placeholder for dataset name based on if subsettting
  observe({
    placeholder <- ifelse(is_subset(), 'eg: QC2 (appended to founder dataset name)', '')
    updateTextInput(session, 'integration_name', placeholder = placeholder)
  })

  excludeOptions <- list(render = I('{option: excludeOptions, item: excludeOptions}'))
  contrastOptions <- list(render = I('{option: contrastOptions, item: contrastItem}'))
  options <- reactiveVal(excludeOptions)

  exclude_choices <- reactive({
    selected <- selected_datasets()

    if (is_subset()) {
      dataset_dir <- file.path(sc_dir, selected)
      clusters <- readRDS(file.path(dataset_dir, 'annot.rds'))
      choices <- get_cluster_choices(clusters, dataset_dir)
      choices$value <- paste(selected, choices$value, sep = '_')
      options(contrastOptions)

    } else {
      colors <- get_palette(selected)
      choices <- get_exclude_choices(selected, sc_dir, colors)
      options(excludeOptions)
    }

    return(choices)
  })

  # update exclude clusters
  observe({
    updateSelectizeInput(session, 'exclude_clusters', choices = exclude_choices(),selected = isolate(input$exclude_clusters), options = options(), server = TRUE)
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
    anals <- c(test(), ctrl())
    req(anals)
    data.frame(sample = anals, pair = NA)
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

    exclude_clusters <- input$exclude_clusters

    if (is_include()) {
      choices <- exclude_choices()
      exclude_clusters <- setdiff(choices$value, exclude_clusters)
    }

    test_anals <- test()
    ctrl_anals <- ctrl()
    datasets <- datasets()
    pairs <- pairs()
    is_subset <- is_subset()

    type <- input$integration_type
    anal_name <- input$integration_name

    error_msg <- validate_integration(test_anals, ctrl_anals, anal_name, datasets, pairs)

    if (is.null(error_msg)) {
      # clear error and disable button
      removeClass('name-container', 'has-error')
      toggleAll(integration_inputs)

      # Create a Progress object
      on.exit(progress$close())
      if (is_subset) {
        progress <- Progress$new(session, min=0, max = 8)
        progress$set(message = "Subsetting dataset", value = 0)

        founder <- get_founder(sc_dir, test_anals)

        subset_saved_scseq(sc_dir,
                           dataset_name = test_anals,
                           exclude_clusters = exclude_clusters,
                           save_name = paste(founder, anal_name, sep = '_'),
                           progress = progress)

      } else {
        progress <- Progress$new(session, min=0, max = 9)
        progress$set(message = "Integrating datasets", value = 0)

        # run integration
        integrate_saved_scseqs(sc_dir,
                               test = test_anals,
                               ctrl = ctrl_anals,
                               exclude_clusters = exclude_clusters,
                               anal_name = anal_name,
                               type = type,
                               pairs = pairs,
                               progress = progress)

      }

      # re-enable, clear inputs, and trigger update of available anals
      ctrl(NULL)
      test(NULL)
      pairs(NULL)
      new_anal(anal_name)
      updateTextInput(session, 'integration_name', value = '')
      toggleAll(integration_inputs)

    } else {
      # show error message
      pairs(NULL)
      html('error_msg', html = error_msg)
      addClass('name-container', class = 'has-error')
    }

  })

  return(new_anal)
}

subset_saved_scseq <- function(sc_dir, dataset_name, save_name, exclude_clusters, progress = NULL) {
  if (is.null(progress)) {
    progress <- list(set = function(value, message = '', detail = '') {
      cat(value, message, detail, '...\n')
    })
  }

  progress$set(1, detail = 'loading')
  scseq <- load_scseqs_for_integration(dataset_name, exclude_clusters, sc_dir)[[1]]
  founder <- get_founder(sc_dir, dataset_name)

  process_raw_scseq(scseq, save_name, sc_dir, progress = progress, value = 1, founder = founder)
}

get_founder <- function(sc_dir, dataset_name) {

  # check for founder of parent
  fpath <- file.path(sc_dir, dataset_name, 'founder.rds')
  founder <- readRDS(fpath)
  if (is.null(founder)) founder <- dataset_name
  return(founder)
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
    req(dataset_dir())

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
    annot(ref_preds)
  })


  # update UI for contrast/cluster choices
  observeEvent(choices(), {

    updateSelectizeInput(session, 'selected_cluster',
                         choices = rbind(NA, choices()),
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

  qc_metrics <- reactive({
    scseq <- scseq()
    if(is.null(scseq)) return(NULL)
    qc <- colnames(scseq@colData)
    qc[match(c('outlier_any',
               'low_lib_size',
               'low_n_features',
               'high_subsets_mito_percent',
               'low_subsets_ribo_percent',
               'high_doublet_score',
               'high_outlyingness',
               'log10_sum',
               'log10_detected',
               'mito_percent',
               'ribo_percent',
               'doublet_score',
               'outlingness'), qc, nomatch = 0)]
  })


  observe({
    sel <- selected_cluster()
    if (!isTruthy(sel)) selected_markers(NULL)
    else selected_markers(markers()[[sel]])
  })


  return(list(
    annot = annot,
    selected_markers = selected_markers,
    selected_cluster = selected_cluster,
    qc_metrics = qc_metrics
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
selectedGene <- function(input, output, session, dataset_name, is.integrated, selected_markers, selected_cluster, type, qc_metrics = function()NULL, ambient = function()NULL) {

  selected_gene <- reactiveVal(NULL)
  gene_options <- list(render = I('{option: geneChoice, item: geneChoice}'))

  exclude_ambient <- reactive({
    if (is.null(input$exclude_ambient)) return(TRUE)
    input$exclude_ambient %% 2 != 1
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
    qc_metrics <- qc_metrics()

    # will error if labels
    # also prevents intermediate redraws
    if ((is.null(markers) & is.null(qc_metrics)) ||
        (is.null(markers) & isTruthy(selected_cluster))) return(NULL)


    get_gene_choices(markers, qc_metrics = qc_metrics)
  })

  # click genecards
  observeEvent(input$genecards, {
    gene_link <- paste0('https://www.genecards.org/cgi-bin/carddisp.pl?gene=', input$selected_gene)
    runjs(paste0("window.open('", gene_link, "')"))
  })


  # reset selected gene if analysis changes
  observe({
    sel <- input$selected_gene
    if (!isTruthy(sel)| !isTruthy(dataset_name())) selected_gene(NULL)
    else selected_gene(sel)
  })


  # update choices
  observe({
    prev <- isolate(input$selected_gene)
    choices <- gene_choices()

    updateSelectizeInput(session, 'selected_gene',
                         choices = choices,
                         server = TRUE,
                         options = gene_options)
  })

  return(list(
    selected_gene = selected_gene,
    exclude_ambient = exclude_ambient,
    show_ridge = show_ridge
  ))

}


#' Logic for cluster plots
#' @export
#' @keywords internal
scClusterPlot <- function(input, output, session, scseq, selected_cluster, dataset_name, fname_fun = function(){}, downloadable = FALSE) {

  plot <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)

    label.highlight <- NULL
    cluster <- selected_cluster()
    cluster <- strsplit(cluster, '-vs-')[[1]]
    if (isTruthy(cluster)) label.highlight <- as.numeric(cluster)

    plot_tsne_cluster(scseq, label.highlight)
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

  observe({
    toggle('plot_container', condition = isTruthy(plot_fun()))
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

#' Logic for marker gene plots
#' @export
#' @keywords internal
scMarkerPlot <- function(input, output, session, scseq, selected_gene, fname_fun = function(){}, downloadable = TRUE) {


  plot <- reactive({
    gene <- selected_gene()
    scseq <- scseq()
    if (!isTruthy(gene) || !isTruthy(scseq)) return(NULL)
    if (!gene %in% row.names(scseq) && !gene %in% colnames(scseq@colData)) return(NULL)
    plot_tsne_feature(scseq, gene)
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

scRidgePlot <- function(input, output, session, selected_gene, selected_cluster, scseq) {

  height <- reactive(length(levels(scseq()$cluster))*50)

  output$ridge_plot <- renderPlot({
    gene <- selected_gene()
    cluster <- selected_cluster()
    req(gene)
    suppressMessages(print(plot_ridge(gene, scseq(), cluster)))
  }, height = height)
}

readRDS.safe <- function(path, .nofile = NULL, .nullfile = NULL) {
  res <- .nofile
  if (isTruthy(path) && file.exists(path))
    res <- readRDS(path)

  if (is.null(res)) return(.nullfile)
  return(res)
}


#' Logic for single cell cluster analyses for Single Cell, Drugs, and Pathways tabs
#' @export
#' @keywords internal
scSampleComparison <- function(input, output, session, dataset_dir, dataset_name, is.integrated = function()TRUE, input_scseq = function()NULL, is_sc = function()TRUE, exclude_ambient = function()FALSE, comparison_type = function()'samples') {
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
    req(is.integrated())
    dir <- dataset_dir()
    if (is.null(dir)) return(NULL)
    get_cluster_choices(annot(),
                        dir,
                        scseq(),
                        sample_comparison = TRUE,
                        top_tables = top_tables(),
                        has_replicates = has_replicates())
  })

  observe({
    updateSelectizeInput(session, 'selected_cluster',
                         choices = cluster_choices(),
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
  sel <- reactive({sel <- input$selected_cluster; req(sel); sel})
  pfun_left_reps <- reactive(function(gene) plot_scseq_gene_medians(gene, annot(), sel(), top_tables(), exclude_ambient()))

  pfun_left_noreps <- reactive(function(gene) list(plot = plot_tsne_feature_sample(gene, scseq(), 'test'), height = 453))
  pfun_left <- reactive(if (has_replicates()) pfun_left_reps() else pfun_left_noreps())

  # plot function for right
  pfun_right_reps <- reactive(function(gene) plot_ridge(gene, scseq(), sel(), by.sample = TRUE, with.height = TRUE))
  pfun_right_noreps  <- reactive(function(gene) list(plot = plot_tsne_feature_sample(gene, scseq(), 'ctrl'), height = 453))
  pfun_right <- reactive(if (has_replicates()) pfun_right_reps() else pfun_right_noreps())

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
      toggleAll(input_ids)

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
      toggleAll(input_ids)
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
      toggleAll(input_ids)
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
      toggleAll(input_ids)
    }

    return(res)
  })

  pairs <- reactive(readRDS.safe(file.path(dataset_dir(), 'pairs.rds')))

  abundances <- reactive(diff_abundance(scseq(), annot(), pairs()))


  # reset lm_fit if selected clusters change
  # also disable runing lm_fit if nothing selected or don't have ctrl and test cells
  comparison_valid <- reactive({
    req(is_sc())
    sel <- as.numeric(input$selected_cluster)
    if (!isTruthy(sel) | !isTruthy(dataset_dir())) return(FALSE)

    stats <- get_cluster_stats(dataset_dir())
    stats$ntest[sel] > 0 & stats$nctrl[sel] > 0
  })


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

  return(list(
    ambient = ambient,
    top_table = top_table,
    cluster_markers = cluster_markers,
    cluster_choices = cluster_choices,
    drug_queries = drug_queries,
    path_res = path_res,
    selected_cluster = reactive(input$selected_cluster),
    annot_clusters = annot_clusters,
    pfun_left = pfun_left,
    pfun_right = pfun_right
  ))
}




#' Format gene plots for sample comparison for drugseqr app
#'
#' @param group Level in \code{scseq$orig.ident} to show cells for. Either \code{'ctrl'} or \code{'test'}
#' @param scseq \code{SingleCellExperiment} object.
#'
#' @return \code{plot} formatted for drugseqr app
#' @export
#' @keywords internal
plot_tsne_feature_sample <- function(gene, scseq, group = 'test') {

  plot <- plot_tsne_feature(scseq, gene)
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


plot_ridge <- function(feature, scseq, selected_cluster, by.sample = FALSE, with.height = FALSE, decreasing = feature %in% c('ribo_percent', 'log10_sum', 'log10_detected')) {
  if (is.null(selected_cluster)) selected_cluster <- ''

  # for one vs one comparisons
  selected_cluster <- strsplit(selected_cluster, '-vs-')[[1]]

  # get color for selected cluster
  clus_levs <- levels(scseq$cluster)
  seli <- as.numeric(selected_cluster)
  sel <- clus_levs[seli]
  nsel <- length(sel)
  colors <- get_palette(clus_levs)
  color  <- colors[seli]

  # either highlight test group or selected cluster
  if (by.sample) {
    scseq <- scseq[, scseq$cluster %in% sel]
    y <- factor(scseq$batch)
    hl <- scseq$orig.ident
    title <- paste('Expression by Sample for', sel)

  } else {
    y <- scseq$cluster
    hl <- as.character(y)
    hl[!hl %in% sel] <- 'out'
    hl <- factor(hl, levels = c(sel, 'out'))
    title <- 'Expression by Cluster'
  }


  if (feature %in% row.names(scseq))
    x <- as.numeric(SingleCellExperiment::logcounts(scseq[feature, ]))
  else
    x <- scseq[[feature]]

  # errors if boolean
  if (is.logical(x)) return(NULL)

  m <- tapply(x, y, mean)
  y <- factor(y, levels = levels(y)[order(m, decreasing = decreasing)])
  df <- data.frame(x, hl, y) %>%
    dplyr::add_count(y) %>%
    dplyr::filter(n > 2)

  pl <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = hl, alpha = hl, color = hl)) +
    ggplot2::scale_fill_manual(values = c(color, 'gray')) +
    ggplot2::scale_alpha_manual(values = c(rep(0.95, nsel), 0.25)) +
    ggplot2::scale_color_manual(values = c(rep('black', nsel), 'gray')) +
    ggridges::geom_density_ridges(scale = 3, rel_min_height = 0.001) +
    ggridges::theme_ridges(center_axis_labels = TRUE) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_dimgray() +
    ggplot2::xlab('') +
    ggplot2::ggtitle(title) +
    ggplot2::theme(legend.position = 'none',
                   plot.title.position = "plot", panel.grid.major.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(color = '#333333', size = 16, face = 'plain', margin = ggplot2::margin(b = 25) ),
                   axis.text.y = ggplot2::element_text(color = '#333333', size = 14),
                   axis.text.x = ggplot2::element_text(color = '#333333', size = 14),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()
    )


  ncells <- tapply(df$x, as.character(df$y), length)
  ncells <- ncells[levels(df$y)]
  xlim <- ggplot2::ggplot_build(pl)$layout$panel_scales_x[[1]]$range$range[2]



  for (i in seq_along(ncells)) {
    pl <- pl + ggplot2::annotation_custom(
      grob = grid::textGrob(label = paste(ncells[i], 'cells'), vjust = -0.3, just = 'right',
                            gp = grid::gpar(fontsize = 14, col = '#333333')),
      ymax = i,
      ymin = i,
      xmax = xlim,
      xmin = xlim
    )
  }

  if (with.height) pl <- list(plot = pl, height = max(453, length(ncells) * 50))
  return(pl)
}


plot_scseq_gene_medians <- function(gene, annot, selected_cluster, tts, exclude_ambient) {

  tt <- lapply(tts, function(x) x[gene, ])
  tt <- do.call(rbind, tt)
  tt <- tt[!is.na(tt$t), ]
  path_df <- get_path_df(tt, path_id = '')

  seli <- which(row.names(tt) == selected_cluster)
  clusters <- as.numeric(row.names(tt))

  # highlight clusters where this gene is significant
  row.names(tt) <- row.names(path_df) <- path_df$Gene <- annot[clusters]

  path_df$ambient <- tt$ambient
  link <- as.character(path_df$Gene)
  path_df$Link <- paste0('<span style="color: dimgray">', link, '</span>')
  path_df$Link[seli] <- gsub('dimgray', 'black', path_df$Link[seli])
  path_df$color <- ifelse(tt$adj.P.Val < 0.05, 'black', 'gray')

  if (exclude_ambient) path_df$color[path_df$ambient] <- 'gray'

  # plot trips up if numbered clusters
  is.number <-  suppressWarnings(!is.na(as.numeric(link)))
  path_df$Gene[is.number] <- paste('Cluster', link[is.number])

  path_df <- path_df[order(abs(tt$dprime), decreasing = TRUE), ]
  dprimesPlotly(path_df, drugs = FALSE)
}


#' Get gslist for pathway analysis
#'
#' @param species Species identifier
#' @param universe NULL
#' @param type either 'go' or 'kegg'
#' @param gs_dir Directory to save gslist to
#'
#' @return gslist
#' @export
get_gslist <- function(species = 'Hs', universe = NULL, type = 'go', gs_dir = '/srv/drugseqr/gs_dir') {

  if (!dir.exists(gs_dir)) dir.create(gs_dir)

  fname <- paste('gslist', type, species, 'rds', sep = '.')
  gslist_path <- file.path(gs_dir, fname)

  if (file.exists(gslist_path)) {
    gslist <- readRDS(gslist_path)

  } else if (type == 'go') {

    #	Get access to package of GO terms
    suppressPackageStartupMessages(OK <- requireNamespace("GO.db",quietly=TRUE))
    if(!OK) stop("GO.db package required but is not installed (or can't be loaded)")

    #	Get access to required annotation functions
    suppressPackageStartupMessages(OK <- requireNamespace("AnnotationDbi",quietly=TRUE))
    if(!OK) stop("AnnotationDbi package required but is not installed (or can't be loaded)")

    #	Load appropriate organism package
    orgPkg <- paste0("org.",species,".eg.db")
    suppressPackageStartupMessages(OK <- requireNamespace(orgPkg,quietly=TRUE))
    if(!OK) stop(orgPkg," package required but is not installed (or can't be loaded)")

    #	Get GO to Entrez Gene mappings
    obj <- paste0("org.",species,".egGO2ALLEGS")
    egGO2ALLEGS <- tryCatch(getFromNamespace(obj,orgPkg), error=function(e) FALSE)
    if(is.logical(egGO2ALLEGS)) stop("Can't find gene ontology mappings in package ",orgPkg)

    gslist <- as.list(egGO2ALLEGS)
    saveRDS(gslist, gslist_path)

  } else if (type == 'kegg') {

    kegg_species <- get_kegg_species(species)
    gkl <- limma::getGeneKEGGLinks(kegg_species, convert = TRUE)
    gkl <- gkl %>%
      dplyr::group_by(PathwayID) %>%
      dplyr::summarise(gslist = list(GeneID))

    gslist <- gkl$gslist
    names(gslist) <- gkl$PathwayID
    saveRDS(gslist, gslist_path)
  }

  return(gslist)
}

#' Convert from species identifier to KEGG species identifier
#'
#' @param species species identifier
#'
#' @return KEGG species identifier
#' @export
get_kegg_species <- function(species) {
  species <- match.arg(species, c("Ag", "At", "Bt", "Ce", "Dm", "Dr", "EcK12", "EcSakai", "Gg", "Hs", "Mm", "Mmu", "Pf", "Pt", "Rn", "Ss", "Xl"))
  #	Convert from Bioconductor to KEGG species codes
  species.KEGG <- switch(species, "Ag"="aga", "At"="ath", "Bt"="bta", "Ce"="cel", "Cf"="cfa", "Dm"="dme", "Dr"="dre", "EcK12"="eco", "EcSakai"="ecs", "Gg"="gga", "Hs"="hsa", "Mm"="mmu", "Mmu"="mcc", "Pf"="pfa", "Pt"="ptr", "Rn"="rno", "Ss"="ssc", "Xl"="xla")

  return(species.KEGG)
}


#' Get names of gene set
#'
#' @param gslist result of \code{\link{get_gslist}}
#' @param type either 'go' or 'kegg'
#' @param species species identifier
#' @param gs_dir Directory to save results to
#'
#' @return Description of \code{gslist} gene sets
#' @export
#'
get_gs.names <- function(gslist, type = 'go', species = 'Hs', gs_dir = '/srv/drugseqr/gs_dir') {
  if (!dir.exists(gs_dir)) dir.create(gs_dir)

  fname <- paste('gs.names', type, species, 'rds', sep = '.')
  gs.names_path <- file.path(gs_dir, fname)

  if (file.exists(gs.names_path)) {
    gs.names <- readRDS(gs.names_path)

  } else if (type == 'go') {
    GOID <- names(gslist)
    TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db,keys=GOID,columns="TERM"))
    gs.names <- TERM$TERM
    names(gs.names) <- TERM$GOID
    saveRDS(gs.names, gs.names_path)

  } else if (type == 'kegg') {
    kegg_species <- get_kegg_species(species)
    gs.names <- limma::getKEGGPathwayNames(kegg_species, remove=TRUE)
    row.names(gs.names) <- gs.names$PathwayID
    gs.names <- gs.names[names(gslist), 'Description']
    names(gs.names) <- names(gslist)
    saveRDS(gs.names, gs.names_path)
  }

  return(gs.names)
}



