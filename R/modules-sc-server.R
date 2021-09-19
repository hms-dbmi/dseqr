
#' Logic for Single Cell Tab
#'
#' @inheritParams bulkPage
#' @export
#'
#' @return Called with \link[shiny]{callModule} to generate logic for
#'   single-cell tab.
#'
scPage <- function(input, output, session, sc_dir, indices_dir, tx2gene_dir, gs_dir, is_mobile, add_sc, remove_sc) {

  # the analysis and options
  scForm <- callModule(
    scForm, 'form',
    sc_dir = sc_dir,
    indices_dir = indices_dir,
    tx2gene_dir = tx2gene_dir,
    gs_dir = gs_dir,
    is_mobile = is_mobile,
    add_sc = add_sc,
    remove_sc = remove_sc)

  # prevent grid differential expression on contrast change
  observeEvent(scForm$compare_groups(), scForm$show_pbulk(FALSE))

  # cluster plot in top right
  clusters_view <- callModule(
    scClusterPlot, 'cluster_plot',
    scseq = scForm$scseq,
    annot = scForm$annot,
    is_mobile = is_mobile,
    clusters_marker_view = clusters_marker_view,
    grid_abundance = grid_abundance,
    grid_expression_fun = scForm$grid_expression_fun,
    selected_gene = scForm$samples_gene,
    show_pbulk = scForm$show_pbulk)

  # cluster comparison plots ---


  clusters_marker_view <- callModule(
    scMarkerPlot, 'marker_plot_cluster',
    scseq = scForm$scseq,
    custom_metrics = scForm$custom_metrics,
    selected_feature = scForm$clusters_gene,
    h5logs = scForm$h5logs,
    show_controls = TRUE,
    is_mobile = is_mobile,
    clusters_view = clusters_view)

  callModule(scBioGpsPlot, 'biogps_plot',
             selected_gene = scForm$clusters_gene,
             species = scForm$species)


  callModule(scViolinPlot, 'violin_plot',
             selected_gene = scForm$clusters_gene,
             selected_cluster = scForm$clusters_cluster,
             scseq = scForm$scseq,
             annot = scForm$annot,
             plots_dir =scForm$plots_dir,
             h5logs = scForm$h5logs,
             is_mobile = is_mobile)


  # sample comparison plots ---

  have_comparison <- reactive(length(scForm$compare_groups()) == 2)

  # grid abundance layer data
  grid_abundance <- callModule(
    scGridAbundance, 'grid_abundance',
    scseq = scForm$scseq,
    dataset_dir = scForm$dataset_dir,
    dplots_dir = scForm$dplots_dir,
    compare_groups = scForm$compare_groups,
    meta = scForm$meta,
    sc_dir = sc_dir)

  test_markers_view <- callModule(
    scMarkerPlot, 'expr_test',
    scseq = scForm$scseq_meta,
    custom_metrics = scForm$custom_metrics,
    selected_feature = scForm$samples_gene,
    h5logs = scForm$h5logs,
    group = 'test',
    is_mobile = is_mobile,
    show_plot = have_comparison,
    clusters_view = clusters_view,
    markers_view = ctrl_markers_view)

  ctrl_markers_view <- callModule(
    scMarkerPlot, 'expr_ctrl',
    scseq = scForm$scseq_meta,
    custom_metrics = scForm$custom_metrics,
    selected_feature = scForm$samples_gene,
    h5logs = scForm$h5logs,
    group = 'ctrl',
    is_mobile = is_mobile,
    show_plot = have_comparison,
    clusters_view = clusters_view,
    markers_view = test_markers_view)


  callModule(scSamplePlot, 'expr_sample_violin',
             selected_gene = scForm$samples_gene,
             plot_fun = scForm$samples_violin_pfun)




  observe({
    toggle(id = "sample_comparison_row",  condition = scForm$comparison_type() == 'samples')
    toggle(id = "cluster_comparison_row", condition = scForm$comparison_type() == 'clusters')
  })


  observe({
    toggle(id = 'biogps_container', condition = !scForm$show_biogps())
    toggle(id = 'violin_container', condition = scForm$show_biogps())
  })

  return(NULL)
}


#' Logic for form on Single Cell Exploration page
#'
#' @keywords internal
#' @noRd
scForm <- function(input, output, session, sc_dir, indices_dir, tx2gene_dir, gs_dir, is_mobile, add_sc, remove_sc) {

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
    resoln_dir <- file.path(dataset_dir, get_resoln_dir(resoln))
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

  dplots_dir <- reactive({
    dataset_dir <- dataset_dir()
    if (is.null(dataset_dir)) return(NULL)
    dplots_dir <- file.path(dataset_dir, 'plots')
    if (!dir.exists(dplots_dir)) dir.create(dplots_dir)
    return(dplots_dir)
  })

  annot <- reactiveVal()
  annot_path <- reactive(file.path(resoln_dir(), 'annot.qs'))

  observe(annot(qread.safe(annot_path())))

  observe(toggle('form_container', condition = scDataset$dataset_exists()))

  # update scseq with cluster changes (from resolution)
  scseq_clusts <- reactive({
    scseq <- scDataset$scseq()
    if (is.null(scseq)) return(NULL)

    clusters <- scResolution$clusters()
    if (!is.null(clusters)) scseq$cluster <- clusters
    return(scseq)
  })


  # read transposed hdf5 logcounts for fast row indexing
  h5logs <- reactive({
    dataset_dir <- dataset_dir()
    if (is.null(dataset_dir)) return(NULL)
    fpath <- file.path(dataset_dir, 'tlogs.tenx')
    res <- HDF5Array::TENxMatrix(fpath, group="mm10")
    return(t(res))
  })

  # dgCMatrix logcounts other
  dgclogs <- reactive({
    dataset_dir <- dataset_dir()
    if (is.null(dataset_dir)) return(NULL)
    progress <- Progress$new(session, min = 0, max = 2)
    progress$set(message = "Loading logcounts", value = 1)
    on.exit(progress$close())
    res <- qs::qread(file.path(dataset_dir, 'dgclogs.qs'))
    progress$set(value = 2)
    return(res)
  })

  counts <- reactive({
    dataset_dir <- dataset_dir()
    if (is.null(dataset_dir)) return(NULL)
    qs::qread(file.path(dataset_dir, 'counts.qs'))
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

  # update scseq with group names
  # gets used for sample comparison marker plots
  scseq_meta <- reactive({
    scseq <- scseq()
    meta <- scSampleGroups$meta()
    groups <- scSampleGroups$groups()

    if (!isTruthyAll(scseq, meta, groups)) return(NULL)
    if (!all(row.names(meta) %in% scseq$batch)) return(NULL)

    attach_meta(scseq, meta = meta, groups = groups)
  })

  qc_metrics <- reactive({
    scseq <- scseq()
    if(is.null(scseq)) return(NULL)
    metrics <- scseq@colData
    qc <- colnames(metrics)

    names(qc) <- sapply(metrics, class)

    qc <- qc[names(qc) %in% c('numeric', 'logical')]
    qc <- qc[!grepl('^sum$|^total$|^subsets|^percent|^sizeFactor|^mapping|^predicted', qc)]
  }, )

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

  # hide sample cluster/feature inputs if not two groups
  have_contrast <- reactive(length(scSampleGroups$groups()) == 2)
  have_cluster <- reactive(isTruthy(scSampleClusters$top_table()))
  observe(toggle(id = "sample_cluster_input",condition = have_contrast()))
  observe(toggle(id = "sample_gene_input",condition = have_contrast()))

  selected_cluster <- reactiveVal('')
  observe({
    type <- comparisonType()
    req(type)

    old <- isolate(selected_cluster())
    new <- switch(type,
                  'clusters' = scClusterComparison$selected_cluster(),
                  'samples' = scSampleClusters$selected_cluster(),
                  'labels' = scLabelsComparison$selected_cluster())

    if (is.null(new)) new <- ''
    if (new != old) selected_cluster(new)
  })


  # the dataset and options
  scDataset <- callModule(scSelectedDataset, 'dataset',
                          sc_dir = sc_dir,
                          new_dataset = new_dataset,
                          indices_dir = indices_dir,
                          tx2gene_dir = tx2gene_dir,
                          add_sc = add_sc,
                          remove_sc = remove_sc)


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
                             counts = counts,
                             dgclogs = dgclogs,
                             snn_graph = scDataset$snn_graph,
                             annot_path = annot_path,
                             show_label_resoln = scDataset$show_label_resoln,
                             compare_groups = scSampleGroups$groups)

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
                                    clusters = scResolution$clusters,
                                    dgclogs = dgclogs)

  scClusterGene <- callModule(selectedGene, 'gene_clusters',
                              scseq = scseq_clusts,
                              h5logs = h5logs,
                              dataset_name = scDataset$dataset_name,
                              resoln_name = scResolution$resoln_name,
                              resoln_dir = resoln_dir,
                              tx2gene_dir = tx2gene_dir,
                              is_integrated = scDataset$is_integrated,
                              cluster_markers = scClusterComparison$cluster_markers,
                              selected_markers = scClusterComparison$selected_markers,
                              selected_cluster = scClusterComparison$selected_cluster,
                              qc_metrics = qc_metrics,
                              type = 'clusters')



  # the selected clusters/gene for sample comparison

  scSampleGroups <- callModule(scSampleGroups, 'sample_groups',
                               dataset_dir = dataset_dir,
                               resoln_dir = resoln_dir,
                               input_scseq = scseq,
                               counts = counts,
                               dataset_name = scDataset$dataset_name,
                               show_pbulk = scSampleGene$show_pbulk)

  scSampleClusters <- callModule(scSampleClusters, 'sample_clusters',
                                 input_scseq = scseq,
                                 meta = scSampleGroups$meta,
                                 h5logs = h5logs,
                                 lm_fit = scSampleGroups$lm_fit,
                                 lm_fit_grid = scSampleGroups$lm_fit_grid,
                                 groups = scSampleGroups$groups,
                                 input_annot = annot,
                                 dataset_dir = dataset_dir,
                                 resoln_dir = resoln_dir,
                                 resoln = resoln,
                                 plots_dir = plots_dir,
                                 dataset_name = scDataset$dataset_name,
                                 sc_dir = sc_dir,
                                 gs_dir = gs_dir,
                                 is_integrated = scDataset$is_integrated,
                                 comparison_type = comparisonType,
                                 exclude_ambient = scSampleGene$exclude_ambient,
                                 applied = scResolution$applied,
                                 is_mobile = is_mobile)

  scSampleGene <- callModule(selectedGene, 'gene_samples',
                             scseq = scDataset$scseq,
                             h5logs = h5logs,
                             dataset_name = scDataset$dataset_name,
                             resoln_name = scResolution$resoln_name,
                             resoln_dir = resoln_dir,
                             tx2gene_dir = tx2gene_dir,
                             is_integrated = scDataset$is_integrated,
                             selected_markers = scSampleClusters$top_table,
                             selected_cluster = scSampleClusters$selected_cluster,
                             type = 'samples',
                             ambient = scSampleClusters$ambient)



  scLabelsComparison <- callModule(scLabelsComparison, 'labels',
                                   cluster_choices = scSampleClusters$cluster_choices)



  return(list(
    scseq = scseq,
    scseq_meta = scseq_meta,
    samples_gene = scSampleGene$selected_gene,
    clusters_gene = scClusterGene$selected_gene,
    custom_metrics = scClusterGene$custom_metrics,
    show_biogps = scClusterGene$show_biogps,
    show_pbulk = scSampleGene$show_pbulk,
    samples_violin_pfun = scSampleClusters$violin_pfun,
    grid_expression_fun = scSampleClusters$grid_expression_fun,
    clusters_cluster = scClusterComparison$selected_cluster,
    samples_cluster = scSampleClusters$selected_cluster,
    labels_cluster = scLabelsComparison$selected_cluster,
    selected_cluster = selected_cluster,
    comparison_type = comparisonType,
    dataset_name = scDataset$dataset_name,
    species = scDataset$species,
    plots_dir = plots_dir,
    dplots_dir = dplots_dir,
    dataset_dir = dataset_dir,
    annot = annot,
    meta = scSampleGroups$meta,
    compare_groups = scSampleGroups$groups,
    h5logs = h5logs
  ))
}


#' Logic for single cell sample comparison groups for Single Cell and Drugs tabs
#'
#' IMPORTANT! USED IN DRUGS TAB:
#' As a result changes here can lead to cryptic bugs in drugs tab.
#'
#' @keywords internal
#' @noRd
#'
scSampleGroups <- function(input, output, session, dataset_dir, resoln_dir, dataset_name, input_scseq = function()NULL, show_pbulk = function()FALSE, counts = function()NULL) {
  group_options <- list(render = I('{option: bulkContrastOptions, item: bulkContrastItem}'))
  input_ids <- c('click_dl_meta', 'click_up_meta', 'compare_groups')


  # need for drugs tab
  scseq <- reactive({
    scseq <- input_scseq()
    if (!is.null(scseq)) return(scseq)

    dataset_dir <- dataset_dir()
    if (!isTruthy(dataset_dir)) return(NULL)

    scseq <- load_scseq_qs(dataset_dir)
    attach_clusters(scseq, resoln_dir())
  })

  prev_path <- reactive({
    if (!isTruthy(dataset_dir())) return(NULL)
    file.path(dataset_dir(), 'prev_groups.qs')
  })

  meta_path <- reactive(file.path(dataset_dir(), 'meta.qs'))

  # download/upload metadata
  observeEvent(input$click_dl_meta, {
    req(scseq())
    shinyjs::click('dl_meta')
  })

  observeEvent(input$click_up_meta, {
    req(scseq())
    shinyjs::click('up_meta')
  })

  ref_meta <- reactive({
    scseq <- scseq()
    samples <- unique(scseq$batch)
    data.frame('Group name' = rep(NA, length(samples)),
               'Pair' = NA,
               row.names = samples,
               check.names = FALSE, stringsAsFactors = FALSE)
  })

  fname <- reactive(paste0(dataset_name(), '_meta.csv'))

  output$dl_meta <- downloadHandler(
    filename = fname,
    content = function(con) {
      utils::write.csv(ref_meta(), con, row.names = TRUE)
    }
  )

  # uploaded annotation
  up_meta <- reactiveVal()
  error_msg <- reactiveVal()

  # reset when change dataset
  observe({
    up_meta(qread.safe(meta_path()))
  })

  observe({
    msg <- error_msg()
    html('error_msg', html = msg)
    toggleClass('validate-up', 'has-error', condition = isTruthy(msg))
  })


  #TODO: convert contrasts to fixed ids based on samples
  observeEvent(input$up_meta, {
    ref <- ref_meta()
    req(ref)

    infile <- input$up_meta
    if (!isTruthy(infile)){
      res <- msg <- NULL

    } else {
      res <- utils::read.csv(infile$datapath, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
      msg <- validate_up_meta(res, ref)

      if (is.null(msg)) {
        colnames(res) <- c('group', 'pair')
        qs::qsave(res, meta_path())

      } else {
        res <- NULL
      }
    }

    error_msg(msg)
    up_meta(res)
  })


  group_choices <- reactive({
    meta <- up_meta()
    if (is.null(meta)) return(NULL)
    groups <- unique(na.exclude(meta$group))
    data.frame(name = groups,
               value = groups, stringsAsFactors = FALSE)
  })

  prev_choices <- reactive({
    qread.safe(prev_path())
  })

  groups <- reactiveVal()
  observeEvent(dataset_name(), groups(NULL))
  observe(groups(input$compare_groups))

  # reset when change resolution
  first_set <- reactiveVal(TRUE)
  observe({
    groups <- groups()
    if (is.null(groups) || groups != 'reset') return(NULL)
    first <- isolate(first_set())
    if (!first)
      updateSelectizeInput(session, 'compare_groups', selected = '')


    first_set(FALSE)
  })

  observe({
    # group_choices may not change with dataset_name change
    dataset_name()

    updateSelectizeInput(session,
                         'compare_groups',
                         choices = group_choices(),
                         selected = prev_choices(),
                         server = TRUE,
                         options = group_options)
  })

  summed <- reactive(qs::qread(file.path(resoln_dir(), 'summed.qs')))
  species <- reactive(qread.safe(file.path(dataset_dir(), 'species.qs')))


  # save groups as previous
  observe({
    groups <- groups()
    prev_path <- isolate(prev_path())
    if (is.null(prev_path)) return(NULL)
    qs::qsave(groups, prev_path)
  })


  summed_grid <- reactive({
    summed_path <- file.path(resoln_dir(), 'summed_grid.qs')

    if (!file.exists(summed_path)){
      # need counts to aggregate
      scseq <- scseq()
      SingleCellExperiment::counts(scseq) <- counts()
      grid <- get_grid(scseq)

      scseq$cluster <- factor(grid$cluster)
      summed <- aggregate_across_cells(scseq)
      qs::qsave(summed, summed_path)
    } else {
      summed <- qs::qread(summed_path)
    }

    return(summed)
  })

  lm_fit <- reactive({
    resoln_dir <- resoln_dir()
    if (is.null(resoln_dir)) return(NULL)

    groups <- input$compare_groups
    if (is.null(groups)) return(NULL)
    if (length(groups) != 2) return(NULL)

    # make sure meta is current
    meta <- up_meta()
    if (!all(groups %in% meta$group)) return(NULL)

    # add hash using uploaded metadata to detect changes
    tohash <- list(meta = meta)
    meta_hash <- digest::digest(tohash, algo = 'murmur32')
    fit_file <- paste0('lm_fit_0svs_', meta_hash, '.qs')
    fit_path <- file.path(resoln_dir, fit_file)

    if (file.exists(fit_path)) {
      fit <- qs::qread(fit_path)

    } else {
      # make sure summed is current
      summed <- summed()
      if (!all(row.names(meta) %in% summed$batch)) return(NULL)

      disableAll(input_ids)
      progress <- Progress$new(session, min = 0, max = 4)
      progress$set(message = "Pseudobulking:", detail = 'clusters', value = 1)
      on.exit(progress$close())

      fit <- run_limma_scseq(summed = summed,
                             meta = meta,
                             species = species(),
                             progress = progress,
                             trend = FALSE,
                             with_fdata = TRUE,
                             min.total.count = 15,
                             min.count = 10)

      progress$set(message = "Saving fits", detail = "", value = 5)
      qs::qsave(fit, fit_path)
      enableAll(input_ids)
    }
    return(fit)
  })


  lm_fit_grid <- reactive({
    if (!show_pbulk()) return(NULL)
    meta <- up_meta()
    if (max(table(meta$group)) < 2) return(NULL)

    # add hash using uploaded metadata to detect changes
    tohash <- list(meta = meta)
    meta_hash <- digest::digest(tohash, algo = 'murmur32')

    fit_file <- paste0('lm_fit_grid_0svs_', meta_hash, '.qs')
    fit_path <- file.path(dataset_dir(), fit_file)

    if (file.exists(fit_path)) {
      lm_fit <- qs::qread(fit_path)

    } else {
      dataset_name <- dataset_name()

      disableAll(input_ids)
      progress <- Progress$new(session, min = 0, max = 5)
      on.exit(progress$close())
      progress$set(message = "Pseudobulking:", detail = 'grid', value = 1)
      summed <- summed_grid()

      lm_fit <- run_limma_scseq(
        summed = summed,
        meta = meta,
        species = species(),
        trend = TRUE,
        method = 'RLE',
        with_fdata = FALSE,
        progress = progress,
        min.total.count = 3,
        min.count = 1)

      progress$set(message = "Saving fits", detail = "", value = 5)
      qs::qsave(lm_fit, fit_path)

      enableAll(input_ids)
    }
    return(lm_fit)
  })

  return(list(
    lm_fit = lm_fit,
    lm_fit_grid = lm_fit_grid,
    groups = groups,
    meta = up_meta,
    scseq = scseq
  ))
}


#' Logic for single cell sample comparison cluster for Single Cell and Drugs tabs
#'
#' IMPORTANT! USED IN DRUGS TAB:
#' As a result changes here can lead to cryptic bugs in drugs tab.
#'
#' @keywords internal
#' @noRd
#'
scSampleClusters <- function(input, output, session, input_scseq, meta, lm_fit, groups, dataset_dir, resoln_dir, resoln, plots_dir, dataset_name, sc_dir, gs_dir = NULL, lm_fit_grid = function()NULL, input_annot = function()NULL, is_integrated = function()TRUE, is_sc = function()TRUE, exclude_ambient = function()FALSE, comparison_type = function()'samples', applied = function()TRUE, is_mobile = function()FALSE, h5logs = function()NULL, page = 'single-cell') {
  cluster_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))
  input_ids <- c('click_dl_anal', 'selected_cluster')

  contrast_dir <- reactiveVal()

  # update contrast_dir if groups or resoln changes
  observe({

    groups <- groups()
    dataset_dir <- isolate(dataset_dir())
    if (length(groups) != 2) {
      cdir <- NULL
    } else {
      contrast <- paste0(groups, collapse = '_vs_')
      resoln <- resoln()
      cdir <- file.path(dataset_dir, get_resoln_dir(resoln), contrast)
    }
    contrast_dir(cdir)
  })

  # set contrast_dir to saved if dataset_dir changes
  observeEvent(dataset_dir(), {
    dataset_dir <- dataset_dir()
    if (is.null(dataset_dir)) {
      contrast_dir(NULL)
      return()
    }

    groups_path <- file.path(dataset_dir, 'prev_groups.qs')
    groups <- qread.safe(groups_path)
    if (is.null(groups)) {
      contrast_dir(NULL)
      return()
    }

    contrast <- paste0(groups, collapse = '_vs_')
    resoln_name <- load_resoln(dataset_dir)
    contrast_dir(file.path(dataset_dir, resoln_name, contrast))
  })

  scseq <- reactive({
    meta <- meta()
    scseq <- input_scseq()
    groups <- groups()
    if (!isTruthyAll(meta, scseq, groups)) return(NULL)
    if (length(groups) != 2) return(NULL)
    if (!all(groups %in% meta$group)) return(NULL)
    if (!all(row.names(meta) %in% scseq$batch)) return(NULL)
    scseq <- attach_meta(scseq, meta=meta, groups=groups)
    scseq <- subset_contrast(scseq)
    return(scseq)
  })


  annot <- reactive({
    annot <- input_annot()
    if (!is.null(annot)) return(annot)

    annot_path <- file.path(resoln_dir(), 'annot.qs')
    if (!file.exists(annot_path)) return(NULL)
    qs::qread(annot_path)
  })

  # update cluster choices in UI
  cluster_choices <- reactive({

    integrated <- is_integrated()
    resoln_dir <- resoln_dir()
    contrast_dir <- contrast_dir()

    if (!isTruthyAll(resoln_dir, integrated, contrast_dir)) return(NULL)
    # need to be in sync (don't take from elsewhere)
    annot_path  <- file.path(resoln_dir, 'annot.qs')
    annot <- qread.safe(annot_path)

    if (is.null(annot)) return(NULL)

    top_tables <- top_tables()
    if (is.null(top_tables)) return(NULL)

    tryCatch({
      choices <- get_cluster_choices(
        clusters = c(annot, 'All Clusters'),
        sample_comparison = TRUE,
        resoln_dir = resoln_dir,
        contrast_dir = contrast_dir,
        use_disk = TRUE,
        top_tables = top_tables)

      choices$disabled <- !choices$value %in% names(top_tables)
      return(choices)

    },
    error = function(e) return(NULL))
  })

  observe({
    updateSelectizeInput(session, 'selected_cluster',
                         choices = rbind(NA, cluster_choices()),
                         options = cluster_options, server = TRUE)
  })


  # path to lmfit, cluster markers, drug query results, and goanna pathway results
  drug_paths <- reactive({
    cdir <- contrast_dir()
    if (is.null(cdir)) return(NULL)
    get_drug_paths(cdir, clusters_str())
  })

  clusters_str <- reactive(collapse_sorted(input$selected_cluster))
  pathways_dir <- reactive(file.path(contrast_dir(), 'pathways'))
  goana_path <- reactive(file.path(pathways_dir(), paste0('goana_', clusters_str(), '.qs')))
  top_tables_paths <- reactive(file.path(pathways_dir(), 'top_tables.qs'))


  # plot functions
  sel <- reactive(input$selected_cluster)

  violin_pfun <- reactive({
    pfun <- function(gene) {
      sel <- sel()
      scseq <- scseq()
      if(!isTruthyAll(sel, gene, scseq)) return(NULL)

      violin_data <- get_violin_data(gene, scseq, sel, by.sample = TRUE, with_all = TRUE, h5logs=h5logs())

      if (all(violin_data$df$x == 0)) return(NULL)
      plot <- plot_violin(violin_data = violin_data, with.height = TRUE, is_mobile = is_mobile())
      return(plot)
    }
    return(pfun)
  })



  grid_expression_fun <- reactive({
    req(is_integrated())
    scseq <- scseq()

    fun <- function(gene) {
      top_tables <- top_tables_grid()
      if (!isTruthyAll(scseq, gene, top_tables)) return(NULL)

      grid <- get_grid(scseq)
      grid_expression <- get_grid_expression(gene, top_tables, grid)
      return(grid_expression)
    }

    return(fun)
  })

  summed <- reactive(qs::qread(file.path(resoln_dir(), 'summed.qs')))

  cluster_ambient <- reactive({
    if (length(groups()) != 2) return(NULL)

    sel <- sel()
    resoln_dir <- resoln_dir()
    if (is.null(resoln_dir)) return(NULL)
    ambience_path <- file.path(resoln_dir, paste0('ambience_', sel, '.qs'))

    if (file.exists(ambience_path)) {
      amb <- qs::qread(ambience_path)

    } else {
      # no cluster e.g. for 'All Clusters'
      summed <- summed()
      if (!sel %in% levels(summed$cluster)) return(NULL)

      # check that not disabled
      disabled <- cluster_choices()[sel, 'disabled']
      if (is.null(disabled) || disabled) return(NULL)

      amb <- get_ambience(scseq())
      if (ncol(amb) != ncol(summed))

        disableAll(input_ids)
      progress <- Progress$new(session, min = 0, max = 2)
      on.exit(progress$close())
      progress$set(message = paste("Detecting ambient genes: cluster", sel), value = 1)

      amb <- calc_cluster_ambience(summed, amb, sel)
      qs::qsave(amb, ambience_path)

      enableAll(input_ids)
    }
    return(amb)
  })

  # differential expression top tables for all clusters
  top_tables <- reactive({

    contrast_dir <- contrast_dir()
    resoln_dir <- resoln_dir()
    if (is.null(contrast_dir)) return(NULL)

    tts_path <- file.path(contrast_dir, 'top_tables.qs')

    if (file.exists(tts_path)) {
      tts <- qs::qread(tts_path)

    } else {
      fit <- lm_fit()
      groups <- groups()
      if (is.null(fit) | is.null(groups)) return(NULL)

      # prevent error when change dataset
      if (!all(groups %in% colnames(fit[[1]]$mod))) return(NULL)

      disableAll(input_ids)
      nmax <- length(fit)+1
      progress <- Progress$new(session, min = 0, max = nmax)
      on.exit(progress$close())
      progress$set(message = "Differential Expression:", value = 1)

      tts <- list()
      for (i in seq_along(fit)) {
        cluster <- names(fit)[i]
        progress$set(detail=paste('cluster', cluster), value = i)

        tt <- tryCatch({
          crossmeta::get_top_table(
            fit[[cluster]],
            groups,
            robust = TRUE,
            allow.no.resid = TRUE)
        },
        error = function(e) NULL)
        if (is.null(tt)) next()


        # need dprime for drug plots and queries
        if (is.null(tt$dprime)) tt$dprime <- tt$logFC
        tts[[cluster]] <- tt
      }

      # add 'All Clusters' result
      progress$set(detail='all clusters', value = nmax)
      annot <-  qs::qread(file.path(resoln_dir, 'annot.qs'))
      all <- as.character(length(annot)+1)
      es <- run_esmeta(tts)

      if (is.null(es)) {
        tts[[all]] <- run_esmeta_logc(tts)

      } else {
        # perform p-val meta analysis for pvals
        # effect size too conservative
        es$pval <- es$fdr <- NA
        pvals <- run_pmeta(tts)
        common <- intersect(row.names(pvals), row.names(es))
        es[common, 'pval'] <- pvals[common, 'p.meta']
        es[common, 'fdr'] <- pvals[common, 'fdr']

        enids <- extract_enids(tts)
        cols <- colnames(tts[[1]])
        tts[[all]] <- es_to_tt(es, enids, cols)
      }

      unlink(contrast_dir, recursive = TRUE)
      dir.create(contrast_dir)
      qs::qsave(tts, tts_path)
      enableAll(input_ids)
    }
    return(tts)
  })

  # differential expression top tables for all 'grid' clusters
  top_tables_grid <- reactive({
    meta <- meta()
    groups <- groups()
    if (is.null(groups) | is.null(meta)) return(NULL)
    meta <- meta[meta$group %in% groups, ]
    if (nrow(meta) < 3) return(NULL)

    contrast_dir <- contrast_dir()
    if (is.null(contrast_dir)) return(NULL)

    tts_path <- file.path(contrast_dir, 'top_tables_grid.qs')

    if (file.exists(tts_path)) {
      tts <- qs::qread(tts_path)

    } else {
      fit <- lm_fit_grid()
      if (is.null(fit)) return(NULL)

      disableAll(input_ids)
      nmax <- length(fit)+1
      progress <- Progress$new(session, min = 0, max = nmax)
      on.exit(progress$close())
      progress$set(message = "Differential Expression:", value = 1)

      tts <- list()
      for (i in seq_along(fit)) {
        cluster <- names(fit)[i]
        progress$set(detail=paste('grid', cluster), value = i)

        tt <- tryCatch({
          crossmeta::get_top_table(
            fit[[cluster]],
            groups,
            trend = TRUE,
            allow.no.resid = TRUE)
        },
        error = function(e) NULL)
        if (is.null(tt)) next()

        tt <- tt[, colnames(tt) %in% c('logFC', 'P.Value'), drop = FALSE]
        tts[[cluster]] <- tt
      }

      qs::qsave(tts, tts_path)
      enableAll(input_ids)
    }
    return(tts)
  })


  # top table for selected cluster only
  top_table <- reactive({
    sel <- input$selected_cluster
    if (!isTruthy(sel)) return(NULL)

    tt <- top_tables()[[sel]]
    if (is.null(tt)) return(NULL)
    if (page == 'drugs') return(tt)

    ambient <- cluster_ambient()
    tt$ambient <- row.names(tt) %in% ambient

    # add ambient-excluded adjusted pvals
    tt$adj.P.Val.Amb[!tt$ambient] <- stats::p.adjust(tt$P.Value[!tt$ambient], method = 'BH')
    if (all(tt$ambient)) tt$adj.P.Val.Amb <- NA

    return(tt)
  })


  # need ambient for pathway
  ambient <- reactive({tt <- top_table() ; row.names(tt)[tt$ambient]})
  species <- reactive(qread.safe(file.path(dataset_dir(), 'species.qs')))

  # indicating if 'All Clusters' selected
  is.meta <- reactive(input$selected_cluster == tail(names(top_tables()), 1))

  # drug query results
  drug_queries <- reactive({
    dpaths <- drug_paths()
    contrast_dir <- contrast_dir()
    if (is.null(contrast_dir)) return(NULL)
    drugs_dir <- file.path(contrast_dir, 'drugs')

    saved_drugs <- any(grepl('^cmap_res_', list.files(drugs_dir)))

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
        paths <- get_drug_paths(contrast_dir(), cluster)
        run_drug_queries(tt, paths, es, species)
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
    goana_path <- goana_path()

    if (file.exists(goana_path)) {
      res <- qs::qread(goana_path)

    } else {
      disableAll(input_ids)
      progress <- Progress$new(session, min = 0, max = 2)
      progress$set(message = "Running pathway analysis", value = 1)
      on.exit(progress$close())

      lm_fit <- lm_fit()
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
        contrast <- paste0(groups(), collapse = '-')

        # no df.residual if 1v1
        have.df <- max(lm_fit$fit$df.residual) > 0
        if (have.df) {
          de <- crossmeta::fit_ebayes(lm_fit, contrast)
        } else {
          de <- top_table()
          de <- de[!de$ambient, ]
        }
      }

      pathways_dir <- pathways_dir()
      dir.create(pathways_dir, showWarnings = FALSE)

      # TODO: run for all clusters at same time
      res <- get_path_res(de, goana_path, gs_dir, species)
      qs::qsave(res, goana_path)

      progress$inc(1)
      enableAll(input_ids)
    }

    return(res)
  })

  pairs <- reactive(qread.safe(file.path(dataset_dir(), 'pairs.qs')))

  abundances <- reactive(diff_abundance(scseq(), annot(), pairs()))



  # enable download
  observe({
    toggleState('click_dl_anal', condition = isTruthy(top_table()))
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
    ab_fname <- 'abundances.csv'
    goup_fname <- 'go_up.csv'
    godn_fname <- 'go_down.csv'

    tt <- top_table()
    if (is.meta()) tt <- tt_to_es(tt)

    pres <- path_res()
    abundances <- abundances()
    tozip <- c()
    tozip <- write.csv.safe(tt, tt_fname, tozip)
    tozip <- write.csv.safe(pres$up, goup_fname, tozip)
    tozip <- write.csv.safe(pres$dn, godn_fname, tozip)
    tozip <- write.csv.safe(abundances, ab_fname, tozip)

    #create the zip file
    utils::zip(file, tozip)
  }

  output$dl_anal <- downloadHandler(
    filename = function() {
      dl_fname()
    },
    content = data_fun
  )

  # download can timeout so get objects before clicking
  observeEvent(input$click_dl_anal, {
    tt <- top_table()
    pres <- path_res()
    shinyjs::click("dl_anal")
  })


  selected_cluster <- reactiveVal()
  observe(selected_cluster(input$selected_cluster))

  return(list(
    ambient = ambient,
    top_table = top_table,
    cluster_choices = cluster_choices,
    drug_queries = drug_queries,
    path_res = path_res,
    selected_cluster = selected_cluster,
    annot_clusters = annot_clusters,
    grid_expression_fun = grid_expression_fun,
    violin_pfun = violin_pfun,
    is_integrated = is_integrated
  ))
}


#' Logic for selected dataset part of scForm
#'
#' @keywords internal
#' @noRd
scSelectedDataset <- function(input, output, session, sc_dir, new_dataset, indices_dir, tx2gene_dir, add_sc, remove_sc) {
  dataset_inputs <- c('selected_dataset', 'show_integration', 'show_label_resoln')
  options <- list(render = I('{option: scDatasetOptions, item: scDatasetItem}'),
                  searchField = c('optgroup', 'label'))

  dataset_name <- reactive({
    ds <- datasets()
    ds$name[ds$value == input$selected_dataset]
  })

  dataset_dir <- reactive(file.path(sc_dir, dataset_name()))
  snn_path <- reactive(file.path(dataset_dir(), 'snn_graph.qs'))

  dataset_exists <- reactive(isTruthy(dataset_name()))

  scseq <- reactive({
    dataset_name <- dataset_name()
    if (!isTruthy(dataset_name)) return(NULL)
    dataset_dir <- dataset_dir()
    load_scseq_qs(dataset_dir)
  })

  # load snn graph
  snn_graph <- reactive({
    snn_path <- snn_path()

    if (file.exists(snn_path)) {
      snn_graph <- qs::qread(snn_path)

    } else {
      scseq <- scseq()
      if (!isTruthy(scseq)) return(NULL)
      is.azi <- file.exists(file.path(dataset_dir(), 'azimuth_ref.qs'))
      if (is.azi) return(NULL)

      disableAll(dataset_inputs)
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
  curr_selected <- reactiveVal()
  datasets <- reactive({
    # reactive to new single cell datasets
    new_dataset()
    datasets <- get_sc_dataset_choices(sc_dir)
    prev <- isolate(prev_datasets())
    curr <- isolate(input$selected_dataset)

    if (isTruthy(prev) && isTruthy(curr)) {
      datasets <- keep_curr_selected(datasets, prev, curr)

      # set currently selected name
      curr_selected(prev[as.numeric(curr), 'name'])
    }

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

  # logic for upload table modal
  up_all <- reactiveVal()
  up_samples <- reactiveVal()

  observeEvent(input$up_raw, {
    prev <- up_all()
    new <- input$up_raw
    up_all(rbind.data.frame(prev, new))
  })

  up_table <- reactive({
    df <- up_all()
    if (is.null(df)) return(NULL)

    df <- df[, c('name', 'size')]
    df$size <- sapply(df$size, utils:::format.object_size, units = 'auto')
    colnames(df) <- c('File', 'Size')

    df <- dplyr::mutate(df, Sample = NA, .before = 1)
    samples <- isolate(up_samples())
    if (!is.null(samples)) df$Sample <- samples
    return(df)
  })

  output$up_table <- DT::renderDataTable({

    table <- up_table()
    if (is.null(table)) return(NULL)
    samples <- up_samples()
    if (!is.null(samples)) table$Sample <- samples

    DT::datatable(table,
                  class = 'cell-border',
                  rownames = FALSE,
                  escape = FALSE, # to allow HTML in table
                  selection = 'multiple',
                  extensions = 'Scroller',
                  options = list(
                    scrollX = TRUE,
                    ordering = FALSE,
                    dom = 't',
                    scroller = TRUE,
                    scrollY = min(250, length(samples)*38)
                  )) %>%  DT::formatStyle('Size', `text-align` = 'right')
  })


  allow_delete <- reactive(isTruthy(input$remove_datasets) & input$confirm_delete == 'delete')

  observe({
    toggleState('delete_dataset', condition = allow_delete())
    toggleClass('delete_dataset', class = 'btn-danger', condition = allow_delete())
  })

  observe({
    toggle('confirm_delete_container', condition = isTruthy(input$remove_datasets))
  })

  observeEvent(input$delete_dataset, {
    remove_datasets <- input$remove_datasets
    unlink(file.path(sc_dir, remove_datasets), recursive = TRUE)
    updateTextInput(session, 'confirm_delete', value = '')
    removeModal()
    new_dataset(remove_datasets)
  })

  observe({
    toggle('sample_name_container', condition = isTruthy(up_table()))
  })

  validate_add_sample <- function(sample, rows) {
    # TODO
    msg <- NULL
    return(msg)
  }

  observeEvent(input$add_sample, {
    sample <- input$sample_name
    rows <- input$up_table_rows_selected
    msg <- validate_add_sample(sample, rows)

    html('error_msg', html = msg)
    toggleClass('validate-up', 'has-error', condition = isTruthy(msg))

    if (is.null(msg)) {
      df <- up_all()
      up_all(df)

      samples <- up_samples()
      samples[rows] <- sample
      up_samples(samples)

      updateTextInput(session, 'sample_name', value = '')
    }
  })


  # open modal selectors
  observeEvent(add_sc(), {
    showModal(uploadModal(session, isTruthy(up_table())))
  })

  observeEvent(remove_sc(), {
    ds <- datasets()
    ds <- ds[!ds$type %in% 'Previous Session', ]
    ds <- tibble::as_tibble(ds)
    names(ds$name) <- ds$name

    choices <- ds %>%
      dplyr::group_by(type) %>%
      dplyr::summarise(names = list(name))

    names(choices$names) <- choices$type
    choices <- choices$names

    showModal(deleteModal(session, choices))
  })


  observeEvent(input$up_raw_progress, {
    # TODO: show table
  })

  # move uploaded to destination
  observeEvent(input$up_raw, {
    new <- input$up_raw
    prev <- up_samples()

    # initialize names using file prefixes
    pat <- paste0(
      '([_ -]+)?',
      c('barcodes[.]tsv(.+)?$',
        'features[.]tsv(.+)?$',
        'genes[.]tsv(.+)?$',
        'matrix[.]mtx(.+)?$',
        'filtered_feature_bc_matrix[.]h5$'), collapse = '|')

    fnames <- new$name
    new <- gsub(pat, '', fnames)
    new[new == ''] <- NA
    new[new == fnames] <- NA

    up_samples(c(prev, new))
  })

  observeEvent(input$click_existing, {
    removeModal()
    Sys.sleep(1)
    shinyjs::click('new_dataset_dir')
  })

  validate_import_samples <- function(up_df, samples) {
    msg <- NULL
    nsamp <- sum(!is.na(samples))
    neach <- table(samples)

    if (nsamp == 0) {
      msg <- 'Specify samples for some files'
    }

    return(msg)
  }

  detect_import_species <- function(up_df) {

    gene.file <- grep('features.tsv|genes.tsv', up_df$name)[1]
    h5.file <- grep('[.]h5$', up_df$name)[1]

    # support only human if no genes.tsv file
    if (is.na(gene.file) & is.na(h5.file)) return("Homo sapiens")

    if (length(h5.file)) {
      infile <- hdf5r::H5File$new(up_df$datapath[h5.file], 'r')
      genes <- infile[['matrix/features/id']][]
      genes <- data.frame(row.names = genes)
    } else {
      genes <- read.table(up_df$datapath[gene.file], row.names = 1)
    }

    get_species(genes)
  }

  # ask for confirmation
  observeEvent(input$import_samples, {

    up_df <- up_all()
    samples <- up_samples()
    msg <- validate_import_samples(up_df, samples)
    species <- detect_import_species(up_df)

    html('error_msg', html = msg)
    toggleClass('validate-up', 'has-error', condition = isTruthy(msg))

    if (is.null(msg)) {
      showModal(confirmModal(session, 'quant', metric_choices, species))
    }
  })


  metric_choices <- c('low_lib_size',
                      'low_n_features',
                      'high_subsets_mito_percent',
                      'low_subsets_ribo_percent',
                      'high_doublet_score')

  # run single-cell quantification
  qargs <- reactiveValues()
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

    up <- up_all()
    samples <- up_samples()
    uniq_samples <- unique(na.omit(samples))

    azimuth_ref <- input$azimuth_ref
    if (!isTruthy(azimuth_ref)) azimuth_ref <- NULL

    for (dataset_name in uniq_samples) {
      upi <- up[samples %in% dataset_name,, drop = FALSE]

      fastq_dir <- file.path(sc_dir, dataset_name)
      unlink(fastq_dir, recursive = TRUE)
      dir.create(fastq_dir)
      file.move(upi$datapath, file.path(fastq_dir, upi$name))

      if (metrics[1] == 'all and none') {
        opts <- list(
          list(dataset_name = paste0(dataset_name, '_QC0'),
               metrics = NULL,
               founder = dataset_name),
          list(dataset_name = paste0(dataset_name, '_QC1'),
               metrics = metric_choices,
               founder = dataset_name))


      } else {
        opts <- list(
          list(dataset_name = dataset_name,
               metrics = metrics,
               founder = NULL))


      }

      # add function that initiates quantification
      # allows to run n at a time
      qargs[[dataset_name]] <- list(
        opts = opts,
        fastq_dir = fastq_dir,
        sc_dir = sc_dir,
        indices_dir = indices_dir,
        tx2gene_dir = tx2gene_dir,
        azimuth_ref = azimuth_ref
      )
    }
  })

  # restrict to two imports at a time
  observe({
    invalidateLater(5000, session)
    todo <- reactiveValuesToList(qargs)
    todo <- names(todo)[!sapply(todo, is.null)]
    if (!length(todo)) return(NULL)

    doing <- reactiveValuesToList(quants)
    doing <- names(doing)[!sapply(doing, is.null)]
    nmax <- 2

    if (length(doing) >= nmax) return(NULL)

    while (length(doing) < nmax && length(todo)) {
      dataset_name <- todo[1]
      todo <- todo[-1]
      doing <- c(doing, dataset_name)

      args <- qargs[[dataset_name]]
      qargs[[dataset_name]] <- NULL

      quants[[dataset_name]] <- callr::r_bg(
        func = run_load_raw_scseq,
        package = 'dseqr',
        args = args
      )

      progress <- Progress$new(max=10*length(args$opts))
      msg <- paste(stringr::str_trunc(dataset_name, 33), "import:")
      progress$set(message = msg, value = 0)
      pquants[[dataset_name]] <- progress
    }
  })

  observe({
    invalidateLater(5000, session)
    handle_sc_progress(quants, pquants, new_dataset)
  })


  observe({
    datasets <- datasets()
    sel <- isolate(input$selected_dataset)

    # for when delete current dataset
    removed.curr <- !is.null(curr_selected()) && !curr_selected() %in% datasets$name
    if (removed.curr) sel <- ''

    datasets <- datasets_to_list(datasets)
    updateSelectizeInput(session, 'selected_dataset', selected = sel, choices = datasets, options = options)
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

  # color label-resoln button when toggled
  observe({
    toggleClass('show_label_resoln', 'btn-primary', condition = show_label_resoln())
  })

  observe({
    toggleClass('show_integration', 'btn-primary', condition = show_integration() | show_subset())
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


#' Logic for single-cell sample comparison plots
#'
#' setup to allow for ggplot/plotly
#'
#'
#' @keywords internal
#' @noRd
scSamplePlot <- function(input, output, session, selected_gene, plot_fun) {

  res <- reactive({
    gene <- selected_gene()
    req(gene)
    suppressMessages(plot_fun()(gene))
  })

  height <- reactive({
    h <- res()$height
    if (is.null(h)) return(1)
    else return(h)
  })


  plot <- reactive({
    res <- res()
    if (is.null(res)) return(NULL)
    return(res$plot)
  })

  output$plot <- shiny::renderPlot(plot(), height = height)
}


#' Logic for label transfer between datasets
#'
#' @keywords internal
#' @noRd
labelTransferForm <- function(input, output, session, sc_dir, dataset_dir, resoln_dir, resoln_name, annot_path, datasets, dataset_name, scseq, species, clusters, show_label_resoln) {
  label_transfer_inputs <- c('transfer_study', 'submit_transfer', 'overwrite_annot', 'ref_name', 'resoln')
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
        # use aggregated reference for speed
        ref_path <- scseq_part_path(sc_dir, ref_subname, 'aggr_ref')
        if (file.exists(ref_path)) {
          ref <- qs::qread(ref_path)

        } else {
          set.seed(100)
          ref <- SingleR::aggregateReference(ref, labels=ref$cluster)
          qs::qsave(ref, ref_path)
        }
        labels <- ref$label
      }
    }

    updateProgress(2/n)
    if (is.null(tab)) {
      # take best label for each cluster
      pred <- SingleR::SingleR(test = query, ref = ref, labels = labels)
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


  # show transferred labels immediately upon selection if have
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
    req(ref_name)

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
resolutionForm <- function(input, output, session, sc_dir, resoln_dir, dataset_dir, dataset_name, scseq, counts, dgclogs, snn_graph, annot_path, show_label_resoln, compare_groups) {
  resolution_inputs <- c('resoln')

  prev_resoln <- reactiveVal()
  resoln_path <- reactiveVal()
  resoln <- reactiveVal()

  observe({
    resoln <- qread.safe(resoln_path())
    prev_resoln(resoln)
  })

  # updateNumericInput removes focus preventing keyboard interaction
  observe({
    clusters <- clusters()
    nclus <- length(levels(clusters))
    type <- ifelse(is_azimuth(), 'nclus_azi', 'nclus')
    shinyjs::html(type, nclus)
  })

  first_set <- reactiveVal(TRUE)
  observeEvent(input[[rname()]], {
    set <- input[[rname()]]
    if (!first_set() && (!is.numstring(set) || set >= 0.1 & set <= 5.1)) {
      resoln(set)

      # prevent update to DE results after change resolution
      compare_groups('reset')
    }

    first_set(FALSE)
  }, ignoreInit = TRUE)

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
    isazi <- file.exists(apath)
    is_azimuth(isazi)

    rpath <- file.path(dataset_dir(), 'resoln.qs')
    resoln_path(rpath)
    init <- qread.safe(rpath, 1)
    resoln(init)

    if (isazi) {
      rname('resoln_azi')
      cols <- colnames(scseq()@colData)
      choices <- get_azimuth_cols(cols, 'cluster')
      updateSelectizeInput(session, 'resoln_azi', choices = choices, selected = init)

    } else {
      rname('resoln')
      updateNumericInput(session, 'resoln', value = init)
    }

  }, priority = 1)

  resoln_name <- reactive(file.path(dataset_name(), get_resoln_dir(resoln())))

  # clusters after change resolution
  clusters_path <- reactive(file.path(resoln_dir(), 'clusters.qs'))

  clusters <- reactive({
    resoln <- resoln()
    clusters_path <- clusters_path()
    clusters <- qread.safe(clusters_path)

    if (!is.null(clusters)) {
      qs::qsave(resoln, resoln_path())

      resoln_name <-  get_resoln_dir(resoln)
      applied_path <- file.path(sc_dir, dataset_name(),resoln_name, 'applied.qs')
      if (file.exists(applied_path)) {
        prev_resoln(resoln)
        return(clusters)
      }

      disableAll(resolution_inputs)
      progress <- Progress$new(session, min = 0, max = 2)
      progress$set(message = "Updating:", detail = 'clusters', value = 1)

    } else {

      g <- snn_graph()
      if (is.null(g)) return(NULL)

      # stop resolution calc when change to dataset with different resolution
      prev_resoln <- prev_resoln()
      if (prev_resoln == resoln) return(NULL)

      qs::qsave(resoln, resoln_path())
      disableAll(resolution_inputs)
      progress <- Progress$new(session, min = 0, max = 2)
      progress$set(message = "Updating:", detail = 'clusters', value = 1)

      clusters <- get_clusters(g, resolution = resoln)
      qs::qsave(clusters, clusters_path)

      # transfer annotation from prev clusters to new
      qs::qsave(levels(clusters), annot_path())
      transfer_prev_annot(resoln, prev_resoln, dataset_name(), sc_dir)
    }


    on.exit(progress$close())
    scseq <- scseq()
    # need counts for pseudobulk
    # need dgclogs for scseq sample (for label transfer)
    SingleCellExperiment::counts(scseq) <- counts()
    SingleCellExperiment::logcounts(scseq) <- dgclogs()

    # add new clusters and run post clustering steps
    scseq$cluster <- clusters
    run_post_cluster(scseq, dataset_name(), sc_dir, resoln, progress, 1, reset_annot = FALSE)

    # mark as previously applied
    prev_resoln(resoln)
    enableAll(resolution_inputs)
    return(clusters)
  })


  return(list(
    clusters = clusters,
    resoln = resoln,
    resoln_name = resoln_name
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
    get_cluster_choices(annot, scseq = scseq)
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

  # set azimuth refs

  species_refs <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)
    species <- scseq@metadata$species
    unname(azimuth_refs[names(azimuth_refs) == species])
  })

  observe({
    species_refs <- species_refs()
    updateSelectizeInput(session, 'azimuth_ref', choices = c('', species_refs))
  })

  observe({
    toggle('azimuth_ref', condition = length(species_refs()) > 0)
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
    error_msg <- validate_subset(selected_dataset(),
                                 input$subset_name,
                                 input$subset_clusters,
                                 is_include(),
                                 hvgs())

    if (is.null(error_msg)) {
      removeClass('name-container', 'has-error')
      showModal(confirmModal(session, 'subset'))

    } else {
      # show error message
      html('error_msg', html = error_msg)
      addClass('name-container', class = 'has-error')
    }
  })

  observeEvent(input$confirm_subset, {
    removeModal()

    subset_clusters <- input$subset_clusters
    subset_name <- input$subset_name
    cluster_choices <- cluster_choices()
    metric_choices <- metric_choices()
    is_include <- is_include()
    from_dataset <- selected_dataset()
    is_integrated <- is_integrated()
    hvgs <- hvgs()

    azimuth_ref <- input$azimuth_ref
    if (!isTruthy(azimuth_ref)) azimuth_ref <- NULL

    exclude_clusters <- intersect(cluster_choices$value, subset_clusters)
    subset_metrics <- intersect(metric_choices$value, subset_clusters)

    if (is_include && length(exclude_clusters)) {
      exclude_clusters <- setdiff(cluster_choices$value, exclude_clusters)
    }

    founder <- get_founder(sc_dir, from_dataset)
    dataset_name <- subsets_name <- paste(founder, subset_name, sep = '_')

    subsets[[dataset_name]] <- callr::r_bg(
      func = subset_saved_scseq,
      package = 'dseqr',
      args = list(
        sc_dir = sc_dir,
        founder = founder,
        from_dataset = from_dataset,
        dataset_name = dataset_name,
        exclude_clusters = exclude_clusters,
        subset_metrics = subset_metrics,
        is_integrated = is_integrated,
        is_include = is_include,
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

  integration_inputs <- c('integration_datasets',
                          'integration_name',
                          'submit_integration',
                          'integration_types',
                          'azimuth_ref')


  integration_name <- reactive(input$integration_name)

  # datasets() with server side selectize causes bug
  integration_choices <- reactive({
    ds <- datasets()
    if (!nrow(ds)) return(NULL)
    int  <- qread.safe(file.path(sc_dir, 'integrated.qs'))
    prev <- qread.safe(file.path(sc_dir, 'prev_dataset.qs'))
    ds <- ds[!ds$name %in% int & !ds$type %in% 'Previous Session', ]

    ds <- tibble::as_tibble(ds)
    names(ds$name) <- ds$name

    choices <- ds %>%
      dplyr::group_by(type) %>%
      dplyr::summarise(names = list(name))

    names(choices$names) <- choices$type
    choices$names
  })

  new_dataset <- reactiveVal()


  is_include <- reactive({ input$toggle_exclude %% 2 != 0 })

  allow_integration <- reactive(length(input$integration_datasets) > 1)


  # show/hide integration forms
  observe({
    toggle(id = "integration-form", anim = TRUE, condition = show_integration())
  })

  # update selected datasets
  observe({
    choices <- integration_choices()
    shinyWidgets::updatePickerInput(session, 'integration_datasets', choices = choices)
  })


  # disable integration if not enough selected
  observe({
    toggleState(id = 'submit_integration', condition = allow_integration())
  })


  # show cluster type choices if enough datasets
  observe(toggle(id = 'integration_types', condition = allow_integration()))


  # show name box only if something selected
  # observe(toggle(id = 'name-container', condition = allow_integration()))

  # set azimuth refs based on species
  species <- reactive({
    get_integration_species(input$integration_datasets, sc_dir)
  })

  species_refs <- reactive({
    species <- species()
    unname(azimuth_refs[names(azimuth_refs) == species])
  })

  enable_azi <- reactive(length(species_refs()) > 0)

  observe(toggle('azimuth_ref', condition = enable_azi()))
  observe({
    disabledChoices <- NULL
    if (!enable_azi()) disabledChoices <- 'Azimuth'

    updateCheckboxGroupButtons(session,
                               'integration_types',
                               disabledChoices = disabledChoices)
  })

  observe({
    updateSelectizeInput(session, 'azimuth_ref', choices = c('', species_refs()))
  })


  # run integration
  pintegs <- reactiveValues()
  integs <- reactiveValues()

  use_azimuth <- reactive({
    'Azimuth' %in% input$integration_types &&
      enable_azi() &&
      allow_integration()
  })

  observe({
    toggle('azimuth_ref_container', condition = use_azimuth())
  })

  observeEvent(input$submit_integration, {

    dataset_names <- input$integration_datasets
    types <- input$integration_types
    name <- input$integration_name

    azimuth_ref <- input$azimuth_ref
    if (!use_azimuth()) azimuth_ref <- NULL

    error_msg <- validate_integration(types, name, azimuth_ref, dataset_names, sc_dir)
    if (is.null(error_msg)) {
      removeClass('name-container', 'has-error')

      integs[[name]] <- callr::r_bg(
        func = run_integrate_saved_scseqs,
        package = 'dseqr',
        args = list(
          sc_dir = sc_dir,
          dataset_names = dataset_names,
          integration_name = name,
          integration_types = types,
          azimuth_ref = azimuth_ref
        )
      )

      progress <- Progress$new(max=8*length(types))
      msg <- paste(stringr::str_trunc(name, 33), "integration:")
      progress$set(message = msg, value = 0)
      pintegs[[name]] <- progress

    } else {
      # show error message
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
clusterComparison <- function(input, output, session, sc_dir, dataset_dir, dataset_name, resoln_dir, resoln, scseq, annot_path, annot, ref_preds, clusters, dgclogs) {
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
    if (!file.exists(annot_path)) return(NULL)
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

  cluster_markers <- reactive({
    resoln_dir <- resoln_dir()
    markers_path <- file.path(resoln_dir, 'markers.qs')

    if (!file.exists(markers_path)) {
      scseq <- scseq()
      if (is.null(scseq)) return(NULL)

      dataset_name <- dataset_name()

      # need dgclogs for presto
      logcounts(scseq) <- dgclogs()

      disableAll(cluster_inputs)
      progress <- Progress$new(session, min = 0, max = 3)
      progress$set(message = "Getting markers", value = 1)
      on.exit(progress$close())

      levels(scseq$cluster) <- seq_along(levels(scseq$cluster))
      markers <- get_presto_markers(scseq)
      resoln_name <- paste0(dataset_name, '/', get_resoln_dir(resoln()))

      progress$set(message = "Saving", value = 2)
      save_scseq_data(list(markers = markers), resoln_name, sc_dir, overwrite = FALSE)
      progress$set(value = 3)
      enableAll(cluster_inputs)

    } else {
      markers <- qs::qread(markers_path)
    }
    return(markers)
  })


  # get/load markers if don't have
  observeEvent(input$selected_cluster, {
    sel <- input$selected_cluster
    resoln_dir <- resoln_dir()
    markers <- markers()
    req(sel, resoln_dir)

    if (sel %in% names(markers)) return(NULL)
    if (!length(markers)) markers <- cluster_markers()

    if (grepl('-vs-', sel)) {
      con_markers <- get_contrast_markers(sel, markers)
      markers[[sel]] <- con_markers
    }

    markers(markers)
  })


  observe({
    sel <- debounce(selected_cluster, 20)
    sel <- sel()

    if (!is.null(sel)) {
      new <- markers()[[sel]]

      prev <- isolate(selected_markers())
      if (is.null(prev) || !identical(new, prev))
        selected_markers(new)
    }
  })


  return(list(
    annot = annot,
    cluster_markers = cluster_markers,
    selected_markers = selected_markers,
    selected_cluster = selected_cluster
  ))
}



#' Logic for selected gene to show plots for
#'
#' @keywords internal
#' @noRd
selectedGene <- function(input, output, session, dataset_name, resoln_name, resoln_dir, tx2gene_dir, scseq, h5logs, is_integrated, selected_markers, selected_cluster, type, cluster_markers = function()NULL, qc_metrics = function()NULL, ambient = function()NULL) {
  gene_options <- list(render = I('{option: geneOption, item: geneItem}'))

  selected_gene <- reactiveVal(NULL)

  # toggle for excluding ambient
  exclude_ambient <- reactive({
    if (is.null(input$exclude_ambient)) return(TRUE)
    input$exclude_ambient %% 2 != 1
  })

  observe({
    toggleClass('exclude_ambient', class = 'btn-primary', condition = exclude_ambient())
  })

  feature <- reactive({
    row <- input$gene_table_rows_selected
    gt <- isolate(gene_table())
    if (is.null(gt) | is.null(row)) return('')
    gt[row]$feature
  })


  gene_selected <- reactive({
    sel <- feature()
    scseq <- scseq()
    if (!isTruthyAll(sel, scseq)) return(FALSE)
    sel %in% row.names(scseq)
  })

  # toggle for violin plot
  have_biogps <- reactive({
    toupper(feature()) %in% biogps[, SYMBOL]
  })

  sel_violin <- reactive(input$show_biogps %% 2 != 1)
  show_biogps <- reactive(sel_violin() | !have_biogps())
  observe(toggleClass(id = "show_biogps", 'btn-primary', condition = !sel_violin()))

  # toggle for showing custom metric
  show_custom_metric <- reactive(type != 'samples' && (input$show_custom_metric %%2 != 0))

  observe({
    toggle('custom_metric_panel', anim = TRUE, condition = show_custom_metric())
    if (show_custom_metric() & have_metric()) selected_gene(input$custom_metric)
  })

  observe(if (!show_custom_metric()) selected_gene(isolate(feature())))
  observe(toggleClass('show_custom_metric', class = 'btn-primary', condition = show_custom_metric()))

  # toggle for showing pseudobulk grid layer
  show_pbulk <- reactiveVal(FALSE)
  observeEvent(input$show_pbulk, show_pbulk(type == 'samples' && !show_pbulk()))

  # prevent grid differential expression on dataset change
  observeEvent(dataset_name(), show_pbulk(FALSE))

  observe(toggleClass(id = "show_pbulk", 'btn-primary', condition = show_pbulk()))


  # disable buttons when not valid
  observe({
    toggleState('show_pbulk', condition = is_integrated())
    toggleState('show_biogps', condition = have_biogps())
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
    # need logcounts
    scseq <- scseq()
    req(metric, scseq, show_custom_metric())

    SingleCellExperiment::logcounts(scseq) <- h5logs()

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


  scseq_genes <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)
    bio <- SummarizedExperiment::rowData(scseq)$bio
    ord <- order(bio, decreasing = TRUE)

    data.frame(row.names = row.names(scseq)[ord])
  })


  # click genecards
  observeEvent(input$genecards, {
    gene_link <- paste0('https://www.genecards.org/cgi-bin/carddisp.pl?gene=', feature())
    runjs(paste0("window.open('", gene_link, "')"))
  })


  # reset selected gene if dataset or cluster changes
  observe({
    sel <- feature()
    if (!isTruthy(sel)| !isTruthy(dataset_name())) return(NULL)
    else selected_gene(sel)
  })

  # reset custom metric if dataset changes
  observeEvent(resoln_name(), custom_metrics(NULL))

  species <- reactive(scseq()@metadata$species)
  tx2gene <- reactive(load_tx2gene(species(), tx2gene_dir))

  qc_first <- reactive(is.null(selected_markers()))

  # update marker genes based on cluster selection
  gene_table <- reactive({

    markers <- selected_markers()
    selected_cluster <- selected_cluster()
    qc_metrics <- qc_metrics()

    # will error if labels
    # also prevents intermediate redraws
    if (is.null(markers) & isTruthy(selected_cluster)) return(NULL)

    if (is.null(markers)) {
      markers <- scseq_genes()
    }

    scseq <- scseq()
    if (is.null(markers) || is.null(scseq)) return(NULL)

    tables <- list()
    tables[[1]] <- get_gene_table(markers, species(), tx2gene())
    tables[[2]] <- get_qc_table(qc_metrics)

    # add other genes
    other <- setdiff(row.names(scseq), tables[[1]]$feature)
    tables[[3]] <- get_leftover_table(other, species(), tx2gene())

    cols <- colnames(tables[[1]])
    if (qc_first()) tables <- tables[c(2,1,3)]

    data.table::rbindlist(tables, fill = TRUE)[, ..cols]
  })

  output$gene_table <- DT::renderDataTable({

    gene_table <- gene_table()
    if (is.null(gene_table)) return(NULL)


    # non-html feature column is hidden and used for search
    # different ncol if contrast
    cols <- colnames(gene_table)
    ncols <- length(cols)
    pct_targs <- grep('%', cols)
    frac_targs <- grep('AUC|logFC', cols)
    pval_targs <- grep('FDR', cols)

    vis_targ <- (length(cols)-1)
    search_targs <- 0

    # prevent sort/filter when qc_first
    sort_targs <- 0
    filter <- list(position='top', clear = FALSE, vertical = TRUE, opacity = 0.85)
    if (qc_first()) {
      sort_targs <- '_all'
      filter = list(position='none')
    }

    if (length(pct_targs)) gene_table[, (pct_targs) := lapply(.SD, as.integer), .SDcols = pct_targs]
    if (length(frac_targs)) gene_table[, (frac_targs) := round(.SD, 2), .SDcols = frac_targs]
    if (length(pval_targs)) gene_table[, (pval_targs) := round(.SD, 3), .SDcols = pval_targs]

    regex <- input$gene_search
    if (grepl(', ', regex)) regex <- format_comma_regex(regex)

    dt <- DT::datatable(
      gene_table,
      class = 'cell-border',
      rownames = FALSE,
      escape = FALSE, # to allow HTML in table
      selection = 'single',
      extensions = 'Scroller',
      filter = filter,
      options = list(
        deferRender = TRUE,
        scroller = TRUE,
        dom = '<"hidden"f>t',
        bInfo = 0,
        scrollY=250,
        search = list(regex = TRUE, search = regex),
        language = list(search = 'Select feature to plot:'),
        columnDefs = list(
          list(visible = FALSE, targets = vis_targ),
          list(searchable = FALSE, targets = search_targs),
          list(sortable = FALSE, targets = sort_targs)
        ),
        headerCallback = DT::JS("function(thead, data, start, end, display) {$('th', thead).css('font-weight', 'normal');}")
      )
    )

    return(dt)

  }, server = TRUE)

  DTproxy <- DT::dataTableProxy("gene_table")

  observe({
    # keep search gene table changes
    regex <- input$gene_search
    if (grepl(', ', regex)) regex <- format_comma_regex(regex)

    DT::updateSearch(DTproxy, keywords = list('global' = regex))
  })

  observe({
    toggle('gene_search_input', condition = !is.null(gene_table()))
  })



  return(list(
    selected_gene = selected_gene,
    exclude_ambient = exclude_ambient,
    show_biogps = show_biogps,
    show_pbulk = show_pbulk,
    custom_metrics = custom_metrics,
    saved_metrics = saved_metrics
  ))


}


#' Logic for cluster plots
#'
#' @keywords internal
#' @noRd
scClusterPlot <- function(input, output, session, scseq, annot, is_mobile, clusters_marker_view, grid_abundance, grid_expression_fun, selected_gene, show_pbulk) {

  show_plot <- reactive(!is.null(scseq()))
  observe(toggle('cluster_plot_container', condition = show_plot()))

  coords <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)

    reds <- SingleCellExperiment::reducedDimNames(scseq)
    reds <- reds[reds %in% c('UMAP', 'TSNE')]
    red <- ifelse('UMAP' %in% reds, 'UMAP', reds[1])

    coords <- SingleCellExperiment::reducedDim(scseq, red)
    as.data.frame(coords)
  })

  labels <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)

    labels <- scseq$cluster
    annot <- levels(labels)
    annot <- format_violin_annot(annot)
    levels(labels) <- annot
    return(labels)
  })

  deck_props <- reactive({
    deck_props <- list()

    if (is_mobile()) {
      deck_props <- list(
        '_typedArrayManagerProps' = list(overAlloc = 1, poolSize = 0)
      )
    }

    return(deck_props)
  })


  label_repels <- reactive({
    coords <- coords()
    scseq <- scseq()
    labels <- as.character(labels())
    if (!isTruthyAll(coords, scseq, labels)) return(NULL)

    # show nums if too many labels/mobile
    label_coords <- get_label_coords(coords, labels)

    if (is_mobile() | length(annot) > 30)
      label_coords$label <- gsub('^(\\d+):.+?$', '\\1', label_coords$label)

    label_repels <- repel::repel_text(
      label_coords,
      xrange = range(coords[,1]),
      yrange = range(coords[,2]),
      direction = 'y')

    return(label_repels)
  })

  label_coords <- reactive({
    label_repels <- label_repels()
    title <- title()
    coords <- coords()
    show_grid <- show_grid()

    if (!isTruthyAll(label_repels, title, coords)) return(NULL)

    nlab <- nrow(label_repels)
    title <- ifelse(show_grid, title, '')

    label_repels <- rbind(
      c(min(coords[,1]), max(coords[,2]), title),
      label_repels)

    label_repels$anchor <- c('start', rep('middle', nlab))
    label_repels$baseline <- c('top', rep('center', nlab))
    label_repels$size <- c(18, rep(14, nlab))

    # only title label if showing grid
    if (show_grid)
      label_repels <- label_repels[1, ]

    return(label_repels)
  })

  show_grid <- reactive({
    show_grid <- input$cluster_plot_show_grid
    ifelse(is.null(show_grid), FALSE, show_grid)
  })

  text_props <- list(getSize=htmlwidgets::JS("d => d.size"),
                     getTextAnchor = htmlwidgets::JS("d => d.anchor"),
                     getAlignmentBaseline = htmlwidgets::JS("d => d.baseline"))

  title <- reactiveVal()

  polygons <- reactive({
    gene <- selected_gene()

    if (show_pbulk() & isTruthy(gene)) {
      grid_expression_fun <- grid_expression_fun()
      polygons <- grid_expression_fun(gene)
      if (is.null(polygons)) return(NULL)

      polygons <- add_grid_colors(polygons)
      title(paste0('Δ EXPRESSION: ', gene))

    } else {
      polygons <- grid_abundance()
      title('Δ CELLS')
    }

    return(polygons)
  })


  output$cluster_plot <- picker::renderPicker({
    coords <- coords()
    labels <- labels()
    deck_props <- deck_props()
    label_coords <- isolate(label_coords())
    polygons <- isolate(polygons())

    if (!isTruthyAll(coords, labels, deck_props)) return(NULL)

    annot <- levels(labels)
    pal <- get_palette(annot)
    names(pal) <- annot
    colors <- unname(pal[labels])

    picker::picker(coords,
                   colors,
                   labels,
                   label_coords,
                   polygons = polygons,
                   point_color_polygons = "white",
                   show_controls = FALSE,
                   deck_props = deck_props,
                   text_props = text_props)
  })

  proxy <- picker::picker_proxy('cluster_plot')
  observe(picker::update_picker(proxy, clusters_marker_view()))
  observe(picker::update_picker(proxy, label_coords = label_coords()))
  observe(picker::update_picker(proxy, polygons = polygons()))

  have_pbulk <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) return(FALSE)
    grid_expression_fun <- grid_expression_fun()
    polygons <- grid_expression_fun(row.names(scseq)[1])
    isTruthy(polygons)
  })

  observe(picker::update_picker(proxy, show_grid = show_pbulk() && have_pbulk()))

  return(reactive(input$cluster_plot_view_state))
}


# get grid abundance data
scGridAbundance <- function(input, output, session, scseq, dataset_dir, sc_dir, compare_groups, dplots_dir, meta) {

  grid_abundance <- reactive({

    dataset_dir <- dataset_dir()
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)

    groups <- compare_groups()
    if (length(groups) != 2) return(NULL)

    meta <- meta()
    if (sum(meta$group %in% groups) < 3) return(NULL)

    scseq <- attach_meta(scseq, meta=meta, groups=groups)
    scseq <- subset_contrast(scseq)

    # add hash for if change groups
    tohash <- list(meta = meta, groups = groups)
    meta_hash <- digest::digest(tohash, algo = 'murmur32')

    apath <- file.path(dplots_dir(), paste0('grid_abundance_', meta_hash, '.qs'))
    if (file.exists(apath)) {
      grid_abundance <- qs::qread(apath)

    } else {
      grid_abundance <- get_grid_abundance(scseq)
      qs::qsave(grid_abundance, apath)
    }

    add_grid_colors(grid_abundance)
  })

  return(grid_abundance)


}


#' Logic for marker feature plots
#'
#' @keywords internal
#' @noRd
scMarkerPlot <- function(input, output, session, scseq, selected_feature, h5logs, clusters_view, is_mobile, show_plot = function()TRUE, markers_view = function()NULL, group = NULL, show_controls = FALSE, deck_props = NULL, custom_metrics = function()NULL) {


  have_colors <- reactive(length(colors()))
  observe(toggle('marker_plot', condition = show_plot() && have_colors()))

  output$marker_plot <- picker::renderPicker({
    if (!show_plot()) return(NULL)

    scseq <- scseq()
    cells <- cells()
    title <- title()
    if (!isTruthy(cells) || !isTruthy(scseq)) return(NULL)

    reds <- SingleCellExperiment::reducedDimNames(scseq)
    reds <- reds[reds %in% c('UMAP', 'TSNE')]
    red <- ifelse('UMAP' %in% reds, 'UMAP', reds[1])

    coords <- SingleCellExperiment::reducedDim(scseq, red)
    coords <- data.frame(coords)

    # use xrange and yrange for all data
    xrange <- range(coords[,1])
    yrange <- range(coords[,2])

    # now subset
    cell.idx <- match(cells, colnames(scseq))
    scseq <- scseq[, cell.idx]
    coords <- coords[cell.idx, ]

    labels <- scseq$cluster
    levels(labels) <- format_violin_annot(levels(labels))

    # show group name as plot label
    label_coords <- NULL

    if (isTruthy(title)) {
      label_coords <- data.frame(x = xrange[1],
                                 y = yrange[2],
                                 label = title)
    }

    deck_props <- list()
    if (is_mobile()) {
      deck_props <- list(
        '_pickable' = FALSE,
        '_typedArrayManagerProps' = list(overAlloc = 1, poolSize = 0)
      )
    }

    picker::picker(coords,
                   isolate(colors()),
                   labels,
                   xrange = xrange,
                   yrange = yrange,
                   # show_controls = show_controls,
                   show_controls = FALSE,
                   label_coords = label_coords,
                   deck_props = deck_props)
  })


  # column data (custom metric + stored)
  cdata <- reactive({
    scseq <- scseq()
    if (!isTruthy(scseq)) return(NULL)

    metrics <- custom_metrics()
    cdata <- scseq@colData

    if (!is.null(metrics) && nrow(cdata) != nrow(metrics)) return(NULL)
    if (!is.null(metrics)) {
      cdata <- cbind(cdata, metrics)
    }

    return(cdata)
  })



  group_title <- if (is.null(group)) '' else toupper(group)
  title <- reactiveVal(group_title)
  cells <- reactiveVal()
  update_colors_proxy <- reactiveVal(TRUE)

  colors <- reactive({
    scseq <- scseq()
    feature <- selected_feature()
    cdata <- cdata()
    if (!isTruthyAll(feature, scseq, cdata)) return(NULL)

    is_gene <- feature %in% row.names(scseq)
    is_feature <- feature %in% colnames(cdata)
    if (!is_gene && !is_feature) return(NULL)

    # get feature
    if (is_gene) {
      ft <- h5logs()[feature, ]

    } else {
      ft <- cdata[[feature]]
      names(ft) <- colnames(scseq)
    }

    # e.g. group or sample columns
    if (is.character(ft) | is.factor(ft)) return(NULL)

    # order cell ids if logical
    ids <- colnames(scseq)
    if (!is.null(group))
      ids <- ids[scseq$orig.ident == group]

    bool.ft <- is.logical(ft)

    set.seed(0)
    if (bool.ft) ids <- ids[order(ft)]
    else ids <- sample(ids)

    prev <- isolate(cells())
    changed.ids <- !identical(prev, ids)
    if (changed.ids) cells(ids)

    update_colors_proxy(!changed.ids)

    # get colors
    ft.ids <- ft[ids]
    all.zero <- all(ft.ids == 0)

    if (bool.ft || all.zero) {
      cols <- const$colors$ft
      ntot <- length(ft.ids)

      colors <- rep(cols[1], ntot)
      colors[ft.ids] <- cols[2]

      # title is info
      ncells <- sum(ft.ids)
      pcells <- round(ncells/ntot*100, 1)

      title(sprintf("%s (%s cells :: %s%%)", feature, ncells, pcells))

    } else {

      # scale before subsetting to group
      ft.scaled <- scales::rescale(ft, c(0, 1))
      ft.scaled <- ft.scaled[ids]

      is_qc <- feature %in% const$features$qc
      reverse_scale <- feature %in% const$features$reverse

      cols <- if (is_qc) const$colors$qc else const$colors$ft
      if (reverse_scale) cols <- rev(cols)
      colors <- scales::seq_gradient_pal(cols[1], cols[2])(ft.scaled)

      # title is group
      prev <- isolate(title())
      if (prev != group_title) title(group_title)
    }

    return(colors)
  })


  title_coords <- reactive({

  })

  observe(picker::update_picker(proxy, label_coords = title_coords()))

  proxy <- picker::picker_proxy('marker_plot')
  observe(picker::update_picker(proxy, clusters_view()))
  observe(picker::update_picker(proxy, markers_view()))

  observe({
    if (!update_colors_proxy()) return(NULL)
    picker::update_picker(proxy, colors = colors())
  })

  return(reactive(input$marker_plot_view_state))
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
    gene <- toupper(selected_gene())
    if (!length(gene)) return(NULL)
    if (!length(species)) return(NULL)
    if (!gene %in% biogps[, SYMBOL]) return(NULL)

    plot_biogps(gene)
  })
}


#' Logic for violin plot for clusters
#'
#' @keywords internal
#' @noRd
scViolinPlot <- function(input, output, session, selected_gene, selected_cluster, scseq, annot, plots_dir, h5logs, is_mobile) {

  show_plot <- reactive(!is.null(plot()))
  observe(toggle('violin_plot', condition = show_plot()))

  height <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) height <- 453
    else height <- max(length(levels(scseq$cluster))*38, 420)
    return(height)
  })


  gene_d <- selected_gene %>% debounce(20)
  clus_d <- selected_cluster %>% debounce(20)

  violin_data <- reactive({
    gene <- gene_d()
    cluster <- clus_d()
    if (!isTruthy(gene)) return(NULL)
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)
    is.gene <- gene %in% row.names(scseq)
    is.num <- is.gene || is.numeric(scseq@colData[[gene]])
    if (!is.num) return(NULL)

    scseq <- scseq()
    h5logs <- if (is.gene) h5logs() else NULL
    if (is.null(scseq)) return(NULL)
    vdat <- get_violin_data(gene, scseq, cluster, with_all = TRUE, h5logs=h5logs)

    return(vdat)
  }) %>% debounce(20)

  plot <- reactive({
    violin_data <- violin_data()
    if (is.null(violin_data)) return(NULL)
    if (all(violin_data$df$x == 0)) return(NULL)
    plot_violin(violin_data = violin_data, is_mobile = is_mobile())
  }) %>% debounce(20)


  output$violin_plot <- renderPlot(plot(), height=height)
}


