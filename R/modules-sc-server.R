
#' Logic for Single Cell Tab
#'
#' @inheritParams bulkPage
#' @export
#'
#' @return Called with \link[shiny]{callModule} to generate logic for
#'   single-cell tab.
#'
scPage <- function(input, output, session, sc_dir, indices_dir, tx2gene_dir, gs_dir, is_mobile, add_sc, remove_sc, integrate_sc, export_sc) {

  # the analysis and options
  scForm <- callModule(
    scForm, 'form',
    sc_dir = sc_dir,
    indices_dir = indices_dir,
    tx2gene_dir = tx2gene_dir,
    gs_dir = gs_dir,
    is_mobile = is_mobile,
    add_sc = add_sc,
    remove_sc = remove_sc,
    integrate_sc = integrate_sc,
    export_sc = export_sc)

  # prevent grid differential expression on contrast change
  observeEvent(scForm$groups(), scForm$show_pbulk(FALSE))

  # cluster plot in top right
  clusters_view <- callModule(
    scClusterPlot, 'cluster_plot',
    scseq = scForm$scseq,
    annot = scForm$annot,
    clusters = scForm$clusters,
    dataset_name = scForm$dataset_name,
    is_mobile = is_mobile,
    clusters_marker_view = clusters_marker_view,
    grid_abundance = grid_abundance,
    grid_expression_fun = scForm$grid_expression_fun,
    selected_gene = scForm$samples_gene,
    show_pbulk = scForm$show_pbulk,
    dataset_dir = scForm$dataset_dir)

  # cluster comparison plots ---


  clusters_marker_view <- callModule(
    scMarkerPlot, 'marker_plot_cluster',
    scseq = scForm$scseq,
    annot = scForm$annot,
    clusters = scForm$clusters,
    added_metrics = scForm$added_metrics,
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
             clusters = scForm$clusters,
             plots_dir =scForm$plots_dir,
             h5logs = scForm$h5logs,
             is_mobile = is_mobile)


  # sample comparison plots ---

  have_comparison <- reactive(length(scForm$groups()) == 2)

  # grid abundance layer data
  grid_abundance <- callModule(
    scGridAbundance, 'grid_abundance',
    scseq = scForm$scseq,
    meta = scForm$meta,
    groups = scForm$groups,
    dplots_dir = scForm$dplots_dir,
    sc_dir = sc_dir)

  test_markers_view <- callModule(
    scMarkerPlot, 'expr_test',
    scseq = scForm$scseq,
    annot = scForm$annot,
    meta = scForm$meta,
    groups = scForm$groups,
    clusters = scForm$clusters,
    added_metrics = scForm$added_metrics,
    selected_feature = scForm$samples_gene,
    h5logs = scForm$h5logs,
    group = 'test',
    is_mobile = is_mobile,
    show_plot = have_comparison,
    clusters_view = clusters_view,
    markers_view = ctrl_markers_view)

  ctrl_markers_view <- callModule(
    scMarkerPlot, 'expr_ctrl',
    scseq = scForm$scseq,
    annot = scForm$annot,
    meta = scForm$meta,
    groups = scForm$groups,
    clusters = scForm$clusters,
    added_metrics = scForm$added_metrics,
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
    toggleCssClass(id = "sample_comparison_row", 'invisible', condition = scForm$comparison_type() != 'samples')
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
scForm <- function(input, output, session, sc_dir, indices_dir, tx2gene_dir, gs_dir, is_mobile, add_sc, remove_sc, integrate_sc, export_sc) {

  set_readonly <- reactive({
    mobile <- is_mobile()
    !is.null(mobile) && mobile
  })

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
    progress$set(message = "Loading:", detail = 'logcounts', value = 1)
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

  resoln <- reactive(scResolution$resoln())

  # update scseq with cluster changes (from resolution)
  scseq_clusts <- reactive({
    scseq <- scDataset$scseq()
    if (is.null(scseq)) return(NULL)

    clusters <- scResolution$clusters()
    if (!is.null(clusters)) scseq$cluster <- clusters

    return(scseq)
  })


  # added metrics (custom and previously saved) needed for violin/marker plots
  added_metrics <- reactive({

    metrics <- list(
      scClusterGene$saved_metrics(),
      scClusterGene$custom_metrics()
    )

    metrics <- metrics[!sapply(metrics, is.null)]

    # to avoid sync issues when switch datasets
    nrows <- sapply(metrics, nrow)
    if (length(unique(nrows)) > 1) return(NULL)

    do.call(cbind, metrics)
  })


  # metrics to add to DT feature table (constant metrics and saved metrics)
  qc_metrics <- reactive({

    scseq <- scDataset$scseq()
    if(is.null(scseq)) return(NULL)

    cdata <- scseq@colData
    metrics <- cdata[, colnames(cdata) %in% c(const$features$qc, const$features$metrics)]

    saved_metrics <- scClusterGene$saved_metrics()
    if (!is.null(saved_metrics)) {
      if (!identical(row.names(metrics), row.names(saved_metrics))) return(NULL)
      metrics <- cbind.safe(metrics, saved_metrics)
    }

    qc <- colnames(metrics)
    names(qc) <- sapply(metrics, class)
    qc <- qc[names(qc) %in% c('numeric', 'logical')]

    samples <- unique(scseq$batch)
    nsamp <- length(samples)
    if (nsamp > 1) {
      names(samples) <- rep('logical', length(samples))
      qc <- c(qc, samples)
    }

    return(qc)
  })

  # show the toggle if dataset is integrated
  observe({
    toggle(id = "comparison_toggle_container",  condition = scDataset$is_integrated())
  })


  # show appropriate inputs based on comparison type
  observe({
    toggle(id = "cluster_comparison_inputs",  condition = comparisonType() == 'clusters')
    toggle(id = "sample_comparison_inputs",  condition = comparisonType() == 'samples')
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
                  'samples' = scSampleClusters$selected_cluster())

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
                          remove_sc = remove_sc,
                          export_sc = export_sc)


  # label transfer between datasets
  # show/hide label transfer forms
  observe({
    toggle(id = "label-resolution-form", anim = TRUE, condition = scDataset$show_label_resoln())
  })

  scLabelTransfer <- callModule(labelTransferForm, 'transfer',
                                sc_dir = sc_dir,
                                tx2gene_dir = tx2gene_dir,
                                set_readonly = set_readonly,
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
                             compare_groups = scSampleGroups$groups,
                             annot = annot)

  # dataset integration
  scIntegration <- callModule(integrationForm, 'integration',
                              sc_dir = sc_dir,
                              tx2gene_dir = tx2gene_dir,
                              datasets = scDataset$datasets,
                              selected_dataset = scDataset$dataset_name,
                              integrate_sc = integrate_sc)

  # dataset subset
  scSubset <- callModule(subsetForm, 'subset',
                         sc_dir = sc_dir,
                         set_readonly = set_readonly,
                         scseq = scseq_clusts,
                         saved_metrics = scClusterGene$saved_metrics,
                         annot = annot,
                         datasets = scDataset$datasets,
                         selected_dataset = scDataset$dataset_name,
                         show_subset = scDataset$show_subset,
                         is_integrated = scDataset$is_integrated,
                         tx2gene_dir = tx2gene_dir)


  # comparison type
  comparisonType <- callModule(comparisonType, 'comparison',
                               is_integrated = scDataset$is_integrated)



  # the selected cluster for cluster comparison
  scClusterComparison <- callModule(clusterComparison, 'cluster',
                                    sc_dir = sc_dir,
                                    set_readonly = set_readonly,
                                    dataset_dir = dataset_dir,
                                    dataset_name = scDataset$dataset_name,
                                    resoln_dir = resoln_dir,
                                    resoln = resoln,
                                    scseq = scseq_clusts,
                                    annot = annot,
                                    annot_path = annot_path,
                                    ref_preds = scLabelTransfer$pred_annot,
                                    clusters = scResolution$clusters,
                                    dgclogs = dgclogs)

  # the selected gene for cluster comparison
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



  # the selected groups for sample comparison
  scSampleGroups <- callModule(scSampleGroups, 'sample_groups',
                               dataset_dir = dataset_dir,
                               resoln_dir = resoln_dir,
                               input_scseq = scseq_clusts,
                               counts = counts,
                               dataset_name = scDataset$dataset_name,
                               show_pbulk = scSampleGene$show_pbulk)


  # the selected cluster for sample comparison
  scSampleClusters <- callModule(scSampleClusters, 'sample_clusters',
                                 input_scseq = scseq_clusts,
                                 set_readonly = set_readonly,
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
                                 tx2gene_dir = tx2gene_dir,
                                 is_integrated = scDataset$is_integrated,
                                 comparison_type = comparisonType,
                                 exclude_ambient = scSampleGene$exclude_ambient,
                                 applied = scResolution$applied,
                                 is_mobile = is_mobile)



  # the selected gene for sample comparison
  scSampleGene <- callModule(selectedGene, 'gene_samples',
                             scseq = scDataset$scseq,
                             h5logs = h5logs,
                             dataset_name = scDataset$dataset_name,
                             resoln_name = scResolution$resoln_name,
                             resoln_dir = resoln_dir,
                             tx2gene_dir = tx2gene_dir,
                             is_integrated = scDataset$is_integrated,
                             can_statistic = scSampleGroups$can_statistic,
                             selected_markers = scSampleClusters$top_table,
                             selected_cluster = scSampleClusters$selected_cluster,
                             type = 'samples',
                             ambient = scSampleClusters$ambient)



  return(list(
    scseq = scDataset$scseq,
    annot = annot,
    meta = scSampleGroups$meta,
    groups = scSampleGroups$groups,
    clusters = scResolution$clusters,
    samples_gene = scSampleGene$selected_gene,
    clusters_gene = scClusterGene$selected_gene,
    added_metrics = added_metrics,
    show_biogps = scClusterGene$show_biogps,
    show_pbulk = scSampleGene$show_pbulk,
    samples_violin_pfun = scSampleClusters$violin_pfun,
    grid_expression_fun = scSampleClusters$grid_expression_fun,
    clusters_cluster = scClusterComparison$selected_cluster,
    samples_cluster = scSampleClusters$selected_cluster,
    selected_cluster = selected_cluster,
    comparison_type = comparisonType,
    dataset_name = scDataset$dataset_name,
    species = scDataset$species,
    plots_dir = plots_dir,
    dplots_dir = dplots_dir,
    dataset_dir = dataset_dir,
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
  input_ids <- c('compare_groups', 'edit_groups')


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


  show_groups_table <- reactiveVal(FALSE)
  observeEvent(input$edit_groups, {

    showing <- show_groups_table()
    if (!showing) {
      show_groups_table(TRUE)
      return()
    }

    res <- rhandsontable::hot_to_r(input$groups_table)
    res <- data.frame(group = res$`Group name`, pair = NA, row.names = res$Sample)
    res[res == ''] <- NA

    no.group <- all(is.na(res$group))

    msg <- validate_up_meta(res, ref_meta())
    prev <- prev_meta()

    valid <- is.null(msg)
    if (no.group | valid) show_groups_table(FALSE)
    if (no.group) return(NULL)

    error_msg(msg)

    if (valid && !identical(res, prev)) {
      prev_meta(res)
      qs::qsave(res, meta_path())
    }
  })

  observe({
    toggle('groups_table_container', condition = show_groups_table())
    toggleCssClass('edit_groups', 'btn-primary',  condition = show_groups_table())
  })

  output$groups_table <- rhandsontable::renderRHandsontable({

    # force re-render on show to avoid disappearing issues
    req(show_groups_table())

    meta <- prev_meta()
    if (is.null(meta)) meta <- ref_meta()

    meta <- data.frame('Sample' = row.names(meta),
                       'Group name' = meta$group, check.names = FALSE)

    rhandsontable::rhandsontable(data = meta,
                                 width = '100%',
                                 height = '200px',
                                 stretchH = "all",
                                 colWidths = c('50%', '50%'),
                                 rowHeaders = FALSE,
                                 contextMenu = FALSE,
                                 manualColumnResize = TRUE) %>%
      rhandsontable::hot_col("Sample", readOnly = TRUE)
  })



  ref_meta <- reactive({
    scseq <- scseq()
    samples <- unique(scseq$batch)
    data.frame(
      group = rep(NA_character_, length(samples)),
      pair = NA_character_,
      check.names = FALSE,
      row.names = samples,
      stringsAsFactors = FALSE)
  })


  # previous annotation
  prev_meta <- reactiveVal()
  error_msg <- reactiveVal()

  # reset when change dataset
  observe({
    prev_meta(qread.safe(meta_path()))
  })

  observe({
    msg <- error_msg()
    html('error_msg', html = msg)
    toggleClass('validate-up', 'has-error', condition = isTruthy(msg))
  })


  group_choices <- reactive({
    meta <- prev_meta()
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
    if (!first) {
      updateSelectizeInput(session, 'compare_groups', selected = '')
    }


    first_set(FALSE)
  })

  observe({
    # group_choices may not change with dataset_name change
    dataset_name()
    choices <- group_choices()

    updateSelectizeInput(session,
                         'compare_groups',
                         choices = choices,
                         selected = prev_choices(),
                         server = TRUE,
                         options = group_options)
  })

  summed <- reactive(qs::qread(file.path(resoln_dir(), 'summed.qs')))
  species <- reactive(qread.safe(file.path(dataset_dir(), 'species.qs')))


  # save groups as previous
  observe({
    groups <- groups()
    if (!is.null(groups) && groups == 'reset') return(NULL)
    prev_path <- isolate(prev_path())
    if (is.null(prev_path)) return(NULL)
    qs::qsave(groups, prev_path)
  })


  summed_grid <- reactive({
    # grid sum is cluster (aka resolution) independent
    summed_path <- file.path(dataset_dir(), 'summed_grid.qs')

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
    meta <- prev_meta()
    if (!all(groups %in% meta$group)) return(NULL)
    meta <- meta[meta$group %in% groups, ]

    # add hash using uploaded metadata to detect changes
    # can re-use fit with change to contrast if groups the same
    meta_hash <- digest::digest(list(meta = meta), algo = 'murmur32')
    fit_file <- paste0('lm_fit_0svs_', meta_hash, '.qs')
    fit_path <- file.path(resoln_dir, fit_file)

    if (file.exists(fit_path)) {
      fit <- qs::qread(fit_path)

    } else {
      # make sure summed is current
      summed <- summed()
      if (!all(row.names(meta) %in% summed$batch)) return(NULL)
      summed <- summed[, summed$batch %in% row.names(meta)]

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

    groups <- input$compare_groups
    if (is.null(groups)) return(NULL)
    if (length(groups) != 2) return(NULL)

    # make sure meta is current
    meta <- prev_meta()
    if (!all(groups %in% meta$group)) return(NULL)
    meta <- meta[meta$group %in% groups, ]

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
      summed <- summed[, summed$batch %in% row.names(meta)]

      lm_fit <- run_limma_scseq(
        summed = summed,
        meta = meta,
        species = species(),
        trend = TRUE,
        method = 'RLE',
        with_fdata = FALSE,
        progress = progress,
        value = 1,
        min.total.count = 3,
        min.count = 1)

      progress$set(message = "Saving fits", detail = "", value = 5)
      qs::qsave(lm_fit, fit_path)

      enableAll(input_ids)
    }
    return(lm_fit)
  })


  can_statistic <- reactive({
    groups <- groups()
    meta <- prev_meta()

    sum(meta$group %in% groups) > 2
  })


  return(list(
    lm_fit = lm_fit,
    lm_fit_grid = lm_fit_grid,
    groups = groups,
    meta = prev_meta,
    can_statistic = can_statistic
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
scSampleClusters <- function(input, output, session, input_scseq, meta, lm_fit, groups, dataset_dir, resoln_dir, resoln, plots_dir, dataset_name, sc_dir, tx2gene_dir, gs_dir = NULL, set_readonly = function()TRUE, lm_fit_grid = function()NULL, input_annot = function()NULL, is_integrated = function()TRUE, is_sc = function()TRUE, exclude_ambient = function()FALSE, comparison_type = function()'samples', applied = function()TRUE, is_mobile = function()FALSE, h5logs = function()NULL, page = 'single-cell') {
  input_ids <- c('click_dl_anal', 'selected_cluster')

  cluster_options <- reactive({
    on_init <- NULL
    if (set_readonly()) on_init <- disableMobileKeyboard(session$ns('selected_cluster'))

    list(render = I('{option: contrastOptions, item: contrastItem}'),
         onInitialize = on_init)
  })


  contrast_dir <- reactiveVal()

  # update contrast_dir if groups or resoln changes
  observe({

    groups <- groups()
    dataset_dir <- isolate(dataset_dir())
    if (length(groups) != 2) {
      cdir <- NULL
    } else {

      # hash meta/groups to detect changed samples for contrast
      meta_hash <- hash_meta(meta(), groups)

      contrast <- paste0(groups, collapse = '_vs_')
      contrast <- paste0(contrast, '_', meta_hash)
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

    # hash meta/groups to get directory
    meta <- qread.safe(file.path(dataset_dir, 'meta.qs'))
    if (is.null(meta)) {
      contrast_dir(NULL)
      return()
    }
    meta_hash <- hash_meta(meta, groups)

    contrast <- paste0(groups, collapse = '_vs_')
    contrast <- paste0(contrast, '_', meta_hash)
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

    annot <- annot()
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
    choices <- cluster_choices()
    if (is.null(choices)) return(NULL)
    updateSelectizeInput(session, 'selected_cluster',
                         choices = rbind(NA, cluster_choices()),
                         options = cluster_options(), server = TRUE)
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
      annot <- annot()
      if(!isTruthyAll(sel, gene, scseq, annot)) return(NULL)

      levels(scseq$cluster) <- annot
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
      progress$set(message = "Detecting ambient genes", detail = paste('cluster', sel), value = 1)

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
    if (page == 'drugs') {
      tt <- top_tables_hs()[[sel]]
      return(list(tt))
    }

    ambient <- cluster_ambient()
    tt$ambient <- row.names(tt) %in% ambient

    # add ambient-excluded adjusted pvals
    tt$adj.P.Val.Amb[!tt$ambient] <- stats::p.adjust(tt$P.Value[!tt$ambient], method = 'BH')
    if (all(tt$ambient)) tt$adj.P.Val.Amb <- NA

    # need as list to check validity in markers table
    tt <- list(tt)
    names(tt) <- sel
    return(tt)
  })


  # need ambient for pathway and drugs
  ambient <- reactive({tt <- top_table()[[1]]; row.names(tt)[tt$ambient]})
  species <- reactive(qread.safe(file.path(dataset_dir(), 'species.qs')))

  # indicating if 'All Clusters' selected
  is.meta <- reactive(input$selected_cluster == tail(names(top_tables()), 1))

  top_tables_hs <- reactive({
    tts <- top_tables()
    species <- species()

    if (species != 'Homo sapiens') {
      # map from species symbols to hgnc
      species_tx2gene <- load_tx2gene(species, tx2gene_dir)
      hsapiens_tx2gene <- load_tx2gene('Homo sapiens', tx2gene_dir)

      symbols <- unique(unlist(lapply(tts, row.names)))

      map <- data.frame(
        row.names = symbols,
        hgnc = species_symbols_to_other(symbols, species_tx2gene, hsapiens_tx2gene)
      )

      # convert row names of top tables to hgnc
      tts <- lapply(tts, function(tt) {
        hgnc <- map[row.names(tt), 'hgnc']

        valid <- !is.na(hgnc) & !duplicated(hgnc)
        tt <- tt[valid, ]
        row.names(tt) <- hgnc[valid]
        return(tt)
      })
    }

    return(tts)
  })

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
      tts <- top_tables_hs()

      for (i in seq_along(tts)) {
        cluster <- names(tts)[i]
        tt <- tts[[cluster]]
        paths <- get_drug_paths(contrast_dir(), cluster)
        run_drug_queries(tt, paths, es)
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
        de <- top_table()[[1]]

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
          de <- top_table()[[1]]
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

    tt <- top_table()[[1]]
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
    tt <- top_table()[[1]]
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


disableMobileKeyboard <- function(id) {
  I(paste0('
        function(){
          $("#', id, '+ .selectize-control input").attr("readonly", "readonly");
        }
    '))

}

#' Logic for selected dataset part of scForm
#'
#' @keywords internal
#' @noRd
scSelectedDataset <- function(input, output, session, sc_dir, new_dataset, indices_dir, tx2gene_dir, add_sc, remove_sc, export_sc) {
  dataset_inputs <- c('selected_dataset', 'show_label_resoln', 'show_subset')

  options <- list(
    render = I('{option: scDatasetOptions, item: scDatasetItem}'),
    searchField = c('optgroup', 'label'))

  dataset_name <- reactiveVal()
  observe({
    sel_idx <- input$selected_dataset
    req(sel_idx)

    ds <- datasets()
    sel <- ds$name[ds$value == sel_idx]
    if (!length(sel)) sel <- NULL

    # catch deleted datasets
    prev <- isolate(dataset_name())
    if (!is.null(prev) && !prev %in% ds$name) sel <- NULL

    if (is.null(prev) || is.null(sel) || sel != prev) dataset_name(sel)
  })

  dataset_dir <- reactive(file.path(sc_dir, dataset_name()))
  snn_path <- reactive(file.path(dataset_dir(), 'snn_graph.qs'))

  dataset_exists <- reactive(isTruthy(dataset_name()))

  scseq <- reactive({
    dataset_dir <- dataset_dir()
    if (!isTruthy(dataset_dir)) return(NULL)
    disableAll(dataset_inputs)

    require(SingleCellExperiment)
    scseq <- load_scseq_qs(dataset_dir)
    gc()
    enableAll(dataset_inputs)

    return(scseq)
  })

  # load snn graph
  snn_graph <- reactive({
    snn_path <- snn_path()

    if (file.exists(snn_path)) {
      snn_graph <- qs::qread(snn_path)

    } else {
      scseq <- scseq()
      if (!isTruthy(scseq)) return(NULL)
      is.ref <- file.exists(file.path(dataset_dir(), 'ref_name.qs'))
      if (is.ref) return(NULL)

      types <- SingleCellExperiment::reducedDimNames(scseq)
      if (!any(c('corrected', 'PCA') %in% types)) return(NULL)

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
    new <- new[file.exists(new$datapath), ]

    up_all(rbind.data.frame(prev, new))
  })


  observeEvent(input$delete_row, {
    selected_row <- as.numeric(strsplit(input$delete_row, "_")[[1]][2])
    df <- up_all()
    samples <- up_samples()

    unlink(df$datapath[selected_row])
    df <- df[-selected_row, ]
    samples <- samples[-selected_row]
    if (!nrow(df)) df <- NULL

    up_all(df)
    up_samples(samples)
    removeClass('validate-up', 'has-error')
  })

  up_table <- reactive({
    df <- up_all()
    if (is.null(df)) return(NULL)

    df <- df[, c('name', 'size')]
    df$size <- sapply(df$size, utils:::format.object_size, units = 'auto')
    colnames(df) <- c('File', 'Size')

    df <- dplyr::mutate(df, ' ' = NA, Sample = NA, .before = 1)
    df$` ` <- getDeleteRowButtons(session, nrow(df))

    samples <- isolate(up_samples())
    if (!is.null(samples)) df$Sample <- samples
    return(df)
  })


  empty_table <- data.frame(' ' = character(0), Sample = character(0), File = character(0), Size = character(0), check.names = FALSE)
  output$up_table <- DT::renderDataTable({

    DT::datatable(empty_table,
                  class = 'cell-border',
                  rownames = FALSE,
                  escape = FALSE, # to allow HTML in table
                  selection = 'multiple',
                  options = list(
                    scrollX = TRUE,
                    ordering = FALSE,
                    dom = 't',
                    paging = FALSE
                  )) %>%
      DT::formatStyle('Size', `text-align` = 'right') %>%
      DT::formatStyle(c('File', 'Size'), color = 'gray')
  })

  proxy <- DT::dataTableProxy('up_table')

  observe({
    shinyjs::toggleCssClass('up_table_container', 'invisible-height', condition = is.null(up_table()))
  })

  observe({
    table <- up_table()
    samples <- up_samples()
    if (!is.null(samples)) table$Sample <- samples

    DT::replaceData(proxy, table, rownames = FALSE)
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
    new_dataset(paste0(remove_datasets, '_delete'))
  })

  observe({
    toggle('sample_name_container', condition = isTruthy(up_table()))
  })

  validate_add_sample <- function(sample, rows) {
    msg <- NULL

    if (is.null(rows)) {
      msg <- 'No rows selected.'
      return(msg)
    }

    if (!isTruthy(sample)) {
      msg <- 'No sample name provided'
      return(msg)
    }

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
    showModal(importSingleCellModal(session, isTruthy(up_table())))
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

    showModal(deleteModal(session, choices, type = 'Single Cell'))
  })

  # get auto sample names
  observeEvent(input$up_raw, {
    new <- input$up_raw
    new <- new[file.exists(new$datapath), ]

    prev <- up_samples()

    # initialize names using file prefixes
    pat <- paste0(
      '([_ -]+)?',
      c('barcodes[.]tsv(.+)?$',
        'features[.]tsv(.+)?$',
        'genes[.]tsv(.+)?$',
        'matrix[.]mtx(.+)?$',
        '[.]rds$',
        '[.]qs$',
        'filtered_feature_bc_matrix(.+)?[.]h(df)?5$',
        'filtered_gene_bc_matrices(.+)?[.]h(df)?5$',
        'raw_gene_bc_matrices(.+)?[.]h(df)?5$',
        '[.]h(df)?5$'
      ), collapse = '|')

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


  # import settings
  detected_species <- reactive({
    # otherwise detected
    up_df <- up_all()

    tryCatch(
      detect_import_species(up_df),
      error = function(e) NULL)
  })

  observe({
    updateSelectizeInput(session, 'import_species', selected = detected_species())
  })

  species_refs <- reactive({
    species <- input$import_species
    if (is.null(species)) return(NULL)

    # reference based not implemented for R object import
    up_df <- up_all()
    if (all(grepl('[.]qs$|[.]rds$', up_df$name))) return(NULL)

    get_refs_list(species)
  })

  robject_import <- reactive(any(grepl('[.]qs$|[.]rds$', up_all()$name)))

  observe({
    toggleClass('confirm_import_datasets', 'disabled', condition = !isTruthy(input$import_species))
  })

  observe({
    have_refs <- length(species_refs()) > 0
    toggle('ref_name_container', condition = have_refs)
  })

  observe({
    updateSelectizeInput(session, 'ref_name', choices = c('', species_refs()))
  })


  # ask for confirmation
  observeEvent(input$import_datasets, {

    up_df <- up_all()
    samples <- up_samples()
    req(up_df)

    msg <- validate_scseq_import(up_df, samples)

    html('error_msg', html = msg)
    toggleClass('validate-up', 'has-error', condition = isTruthy(msg))

    if (!is.null(msg)) return(NULL)

    showModal(confirmImportSingleCellModal(
      session,
      const$features$metrics,
      detected_species(),
      species_refs(),
      warn_robject = robject_import()))

  })

  observe({
    toggleState('import_datasets', condition = !is.null(up_all()))
  })


  # run single-cell quantification
  qargs <- reactiveValues()
  quants <- reactiveValues()
  pquants <- reactiveValues()
  deselect_dataset <- reactiveVal(0)

  observeEvent(input$confirm_import_datasets, {
    species <- input$import_species
    if (!isTruthy(species)) return(NULL)

    metrics <- input$qc_metrics
    # none, all, all and none: can't combine
    if (length(metrics) > 1 && !all(metrics %in% const$features$metrics)) return(NULL)

    if (!isTruthy(metrics)) metrics <- 'none'
    if (metrics[1] == 'all') metrics <- const$features$metrics

    removeModal()

    up <- up_all()
    samples <- up_samples()
    uniq_samples <- unique(na.omit(samples))

    ref_name <- input$ref_name
    if (!isTruthy(ref_name)) ref_name <- NULL

    for (dataset_name in uniq_samples) {
      upi <- up[samples %in% dataset_name,, drop = FALSE]

      uploaded_data_dir <- file.path(sc_dir, dataset_name)
      unlink(uploaded_data_dir, recursive = TRUE)
      dir.create(uploaded_data_dir)
      file.move(upi$datapath, file.path(uploaded_data_dir, upi$name))

      if (metrics[1] == 'all and none') {
        opts <- list(
          list(dataset_name = paste0(dataset_name, '_QC0'),
               metrics = NULL,
               founder = dataset_name),
          list(dataset_name = paste0(dataset_name, '_QC1'),
               metrics = const$features$metrics,
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
        uploaded_data_dir = uploaded_data_dir,
        sc_dir = sc_dir,
        indices_dir = indices_dir,
        tx2gene_dir = tx2gene_dir,
        ref_name = ref_name,
        species = species
      )
    }

    # clear uploaded
    up_all(NULL)
    up_samples(NULL)
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
        func = run_import_scseq,
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

  # handle dataset export
  observeEvent(export_sc(), {
    datasets <- datasets()
    datasets <- datasets_to_list(datasets)

    sel <- input$selected_dataset

    showModal(exportModal(session, choices = datasets, selected = sel, options = options))
  })


  scseq_export <- reactiveVal()

  export_name <- reactive({
    ds <- datasets()
    sel_idx <- input$export_dataset
    req(sel_idx)

    ds$name[ds$value == sel_idx]
  })

  observeEvent(input$confirm_export, {
    disable('confirm_export')

    dataset_dir <- file.path(sc_dir, export_name())

    progress <- Progress$new(session, min = 0, max = 3)
    progress$set(message = "Preparing export:", detail = export_name(), value = 1)
    on.exit(progress$close())

    # load scseq
    scseq <- load_scseq_qs(dataset_dir, with_counts = TRUE, with_logs = TRUE)
    scseq <- prep_scseq_export(scseq, dataset_dir)

    # save to disk
    progress$set(message = "Saving:", detail = paste0(export_name(), '.qs'), value = 2)
    fpath <- tempfile()
    qs::qsave(scseq, fpath)
    scseq_export(fpath)

    progress$set(3)
    shinyjs::click('download_dataset')
    enable('confirm_export')
  })

  output$download_dataset <- downloadHandler(
    filename = function() {
      paste0(export_name(), '.qs')
    },
    content = function(file) {
      file.copy(scseq_export(), file)
    }
  )


  # show/hide integration/label-transfer forms
  show_subset <- reactive(input$show_subset %% 2 == 1)
  show_label_resoln <- reactive(input$show_label_resoln %% 2 == 1)

  # hide integration/label-transfer buttons no dataset
  observe({
    toggle('show_label_resoln-parent', condition = dataset_exists())
  })

  # color label-resoln button when toggled
  observe({
    toggleClass('show_label_resoln', 'btn-primary', condition = show_label_resoln())
  })

  observe({
    toggleClass('show_subset', 'btn-primary', condition = show_subset())
  })


  return(list(
    dataset_name = dataset_name,
    scseq = scseq,
    snn_graph = snn_graph,
    datasets = datasets,
    show_subset = show_subset,
    show_label_resoln = show_label_resoln,
    is_integrated = is_integrated,
    dataset_exists = dataset_exists,
    species = species
  ))
}

get_refs_list <- function(species) {

  ref_names <- refs$name
  names(ref_names) <- refs$label
  split(ref_names, refs$type)
}

prep_scseq_export <- function(scseq, dataset_dir) {

  # store selected resolution
  resoln_path <- file.path(dataset_dir, 'resoln.qs')
  scseq@metadata$resoln <- qread.safe(resoln_path)

  # store reference name
  ref_path <- file.path(dataset_dir, 'ref_name.qs')
  scseq@metadata$ref_name <- qread.safe(ref_path)

  # store group metadata
  meta_path <- file.path(dataset_dir, 'meta.qs')
  scseq@metadata$meta <- qread.safe(meta_path)

  # store contrast groups
  group_path <- file.path(dataset_dir, 'prev_groups.qs')
  scseq@metadata$prev_groups <- qread.safe(group_path)

  # add cluster annotations
  resoln_dir <- load_resoln(dataset_dir)
  annot <- qs::qread(file.path(dataset_dir, resoln_dir, 'annot.qs'))
  levels(scseq$cluster) <- annot

  return(scseq)
}


detect_import_species <- function(up_df) {

  gene.file <- grep('features.tsv|genes.tsv', up_df$name)[1]
  h5.file <- grep('[.]h5$|[.]hdf5$', up_df$name)[1]

  # get user selection if can't detect
  if (is.na(gene.file) & is.na(h5.file)) return(NULL)

  if (!is.na(h5.file)) {
    infile <- hdf5r::H5File$new(up_df$datapath[h5.file], 'r')
    genomes <- names(infile)
    slot <- ifelse(hdf5r::existsGroup(infile, "matrix"), 'features/id', 'genes')

    genes <- infile[[file.path(genomes[1], slot)]][]
    genes <- data.frame(row.names = genes)
  } else {
    genes <- read.table(up_df$datapath[gene.file])
    genes <- genes[!is.na(genes$V1), ]
    row.names(genes) <- genes$V1
  }

  get_species(genes)
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
labelTransferForm <- function(input, output, session, sc_dir, tx2gene_dir, set_readonly, dataset_dir, resoln_dir, resoln_name, annot_path, datasets, dataset_name, scseq, species, clusters, show_label_resoln) {
  label_transfer_inputs <- c('overwrite_annot', 'ref_name', 'sc-form-resolution-resoln', 'sc-form-resolution-resoln_ref')
  asis <- c(FALSE, FALSE, TRUE, TRUE)

  disabled_demo <- getShinyOption('is_example', FALSE)
  observe(if (disabled_demo) addClass('overwrite_annot', 'disabled fa-disabled'))

  options <-  reactive({
    on_init <- NULL
    if (set_readonly()) on_init <- disableMobileKeyboard(session$ns('ref_name'))

    list(render = I('{option: transferLabelOption, item: scDatasetItemDF}'), onInitialize = on_init)
  })

  ref_preds <- reactiveVal()
  new_preds <- reactiveVal()
  new_annot <- reactiveVal()

  preds_dir <- reactive(file.path(resoln_dir(), 'preds'))


  # saved label transfer predictions
  preds <- reactive({
    new_preds()

    preds_dir <- preds_dir()
    pred_files <- list.files(preds_dir)

    preds <- list()
    for (pred_file in pred_files) {
      pred_name <- tools::file_path_sans_ext(pred_file)
      preds[[pred_name]] <- qs::qread(file.path(preds_dir, pred_file))
    }

    preds <- validate_preds(preds, sc_dir)
    return(preds)
  })

  observeEvent(resoln_name(), new_preds(NULL))

  # update annotation transfer choices
  observe({
    preds <- preds()

    datasets <- datasets()
    dataset_name <- dataset_name()
    req(preds, datasets)

    transfer_name <- new_preds()
    selected <- get_selected_from_transfer_name(transfer_name, dataset_name)

    choices <- get_label_transfer_choices(datasets, dataset_name, preds)
    updateSelectizeInput(session,
                         'ref_name',
                         choices = choices,
                         server = TRUE,
                         selected = selected,
                         options = options())
  })



  query <- reactive({
    query_path <- scseq_part_path(sc_dir, resoln_name(), 'scseq_sample')
    if (!file.exists(query_path)) return(NULL)

    qs::qread(query_path)
  })

  transfers <- reactiveValues()
  ptransfers <- reactiveValues()
  is_disabled <- reactiveVal(FALSE)

  # submit annotation transfer
  observeEvent(input$ref_name, {

    query_name <- dataset_name()
    ref_name <- input$ref_name
    resoln_name <- resoln_name()
    preds <- preds()

    req(ref_name != 'reset')
    req(query_name, ref_name, preds)
    req(!ref_name %in% names(preds))
    req(show_label_resoln())

    transfer_name <- paste0(ref_name, ' → ', query_name)

    disableAll(label_transfer_inputs, asis)
    is_disabled(TRUE)

    transfers[[transfer_name]] <- callr::r_bg(
      func = run_label_transfer,
      package = 'dseqr',
      args = list(
        sc_dir = sc_dir,
        tx2gene_dir = tx2gene_dir,
        resoln_name = resoln_name,
        query_name = query_name,
        ref_name = ref_name
      )
    )

    progress <- Progress$new(max=3)
    progress$set(message = paste0(transfer_name, ':'), value = 0)
    ptransfers[[transfer_name]] <- progress

  })

  # enable when complete (only one transfer at a time)
  observe({
    invalidateLater(1000, session)
    doing <- reactiveValuesToList(transfers)
    doing <- names(doing)[!sapply(doing, is.null)]

    if (!length(doing) && is_disabled()) {
      enableAll(label_transfer_inputs, asis)
      is_disabled(FALSE)
    }

  })

  observe({
    invalidateLater(5000, session)
    handle_sc_progress(transfers, ptransfers, new_preds)
  })


  # show transferred labels immediately upon selection if have
  observe({
    query_name <- resoln_name()
    ref_name <- input$ref_name
    req(query_name)
    preds <- preds()

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
    if (disabled_demo) return(NULL)
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

# TODO: save SingleR result as resolution independent and perform table on demand
run_label_transfer <- function(sc_dir, tx2gene_dir, resoln_name, query_name, ref_name, progress = NULL, value = 0) {

  if (is.null(progress)) {
    progress <- list(set = function(value, message = '', detail = '') {
      cat(value, message, detail, '...\n')
    })
  }

  query_path <- scseq_part_path(sc_dir, resoln_name, 'scseq_sample')
  preds_dir <- file.path(sc_dir, resoln_name, 'preds')
  dir.create(preds_dir, showWarnings = FALSE)

  preds_path <- file.path(preds_dir, paste0(ref_name, '.qs'))
  query <- qs::qread(query_path)

  # get arguments for SingleR
  tab <- NULL
  ref_date <- NULL
  senv <- loadNamespace('celldex')

  progress$set(value+1, detail = 'getting reference')

  if (ref_name %in% ls(senv)) {
    ref <- get(ref_name, envir = senv)()
    ref@metadata$species <- get_celldex_species(ref_name)
    labels <- ref$label.fine

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

      # aggregation removes metadata
      ref_species <- ref@metadata$species

      if (file.exists(ref_path)) {
        ref <- qs::qread(ref_path)

      } else {
        set.seed(100)
        ref <- SingleR::aggregateReference(ref, labels=ref$cluster)
        qs::qsave(ref, ref_path)
      }
      labels <- ref$label
      ref@metadata$species <- ref_species
    }
  }

  progress$set(value+2, detail = 'predicting')

  if (is.null(tab)) {

    # use homologous hgnc symbols if not the same species
    if (query@metadata$species != ref@metadata$species) {
      ref <- convert_species(ref, tx2gene_dir)
      query <- convert_species(query, tx2gene_dir)
    }

    # take best label for each cluster
    preds <- SingleR::SingleR(test = query, ref = ref, labels = labels)
    tab <- table(assigned = preds$pruned.labels, cluster = query$cluster)
  }

  preds <- row.names(tab)[apply(tab, 2, which.max)]

  # keep track of date that reference was used so that can invalidate if overwritten
  attr(preds, 'ref_date') <- as.numeric(ref_date)

  qs::qsave(preds, preds_path)
  return(TRUE)
}



#' Logic for leiden resolution slider
#'
#' @keywords internal
#' @noRd
resolutionForm <- function(input, output, session, sc_dir, resoln_dir, dataset_dir, dataset_name, scseq, counts, dgclogs, snn_graph, annot_path, show_label_resoln, compare_groups, annot) {
  resolution_inputs <- c('resoln', 'resoln_ref')

  disabled_demo <- getShinyOption('is_example', FALSE)
  observe(if (disabled_demo) {
    disable('resoln')
    disable('resoln_ref')
  })

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
    type <- ifelse(is_reference(), 'nclus_ref', 'nclus')
    shinyjs::html(type, nclus)
  })

  observeEvent(input[[rname()]], {
    set <- input[[rname()]]
    if (set == '') return(NULL)
    if (!is.numstring(set) || set >= 0.1 & set <= 5.1) {
      resoln(set)

      # prevent update to DE results after change resolution
      compare_groups('reset')
    }

  }, ignoreInit = TRUE)

  rname <- reactiveVal('fixed')
  is_reference <- reactiveVal(FALSE)

  observe({
    shinyjs::toggle('resoln_container', condition=!is_reference())
    shinyjs::toggle('resoln_ref_container', condition=is_reference())
  })

  observeEvent(resoln_dir(), {

    fpath <- file.path(resoln_dir(), 'provided_clusters.qs')
    shinyjs::toggle('provided_clusters_warning', condition = file.exists(fpath))
  })

  observeEvent(dataset_dir(),  {

    dataset_dir <- dataset_dir()
    req(dataset_dir)

    # restrict resolutions if reference
    ref_path <- file.path(dataset_dir(), 'ref_name.qs')
    is.ref <- file.exists(ref_path)
    is_reference(is.ref)

    rpath <- file.path(dataset_dir(), 'resoln.qs')
    resoln_path(rpath)
    init <- qread.safe(rpath, 1)

    resoln(init)
    resoln_fixed <- init == 'provided.clusters'

    if (resoln_fixed) {
      rname('fixed')

    } else if (is.ref) {
      rname('resoln_ref')
      cols <- colnames(scseq()@colData)
      choices <- get_ref_cols(cols, 'cluster')
      updateSelectizeInput(session, 'resoln_ref', choices = choices, selected = init)

    } else {
      rname('resoln')
      updateSliderInput(session, 'resoln', value = init)
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
      on.exit(progress$close())

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
      on.exit(progress$close())

      clusters <- get_clusters(g, resolution = resoln)
      qs::qsave(clusters, clusters_path)

      # transfer annotation from prev clusters to new
      qs::qsave(levels(clusters), annot_path())
      annot <- transfer_prev_annot(resoln, prev_resoln, dataset_name(), sc_dir)
      annot(annot)
    }


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
subsetForm <- function(input, output, session, sc_dir, set_readonly, scseq, saved_metrics, annot, datasets, show_subset, selected_dataset, cluster_choices, is_integrated, tx2gene_dir) {
  subset_inputs <- c('subset_name', 'submit_subset', 'subset_features', 'toggle_exclude', 'click_up')
  type <- name <- NULL

  disabled_demo <- getShinyOption('is_example', FALSE)
  observe(if (disabled_demo){
    addClass('submit_subset', 'disabled fa-disabled')
  })

  contrastOptions <- reactive({
    on_init <- NULL
    if (set_readonly()) on_init <- disableMobileKeyboard(session$ns('subset_features'))

    list(render = I('{option: contrastOptions, item: contrastItem}'),
         onInitialize = on_init)

  })

  subset_name <- reactive(input$subset_name)
  new_dataset <- reactiveVal()


  # show/hide integration forms
  observe({
    toggle(id = "subset-form", anim = TRUE, condition = show_subset())
  })

  is_include <- reactive({ input$toggle_exclude %% 2 != 0 })

  cluster_choices <- reactive({
    scseq <- scseq()
    annot <- annot()
    if (is.null(scseq) | is.null(annot)) return(NULL)
    get_cluster_choices(annot, scseq = scseq)
  })

  metric_choices <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)

    saved_metrics <- saved_metrics()
    if (!is.null(saved_metrics)) {
      if (!identical(row.names(saved_metrics), colnames(scseq))) return(NULL)
      scseq@colData <- cbind(scseq@colData, saved_metrics)
    }
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

  # set reference names
  species_refs <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)
    species <- scseq@metadata$species
    get_refs_list(species)
  })

  observe({
    species_refs <- species_refs()
    updateSelectizeInput(session, 'ref_name', choices = c('', species_refs))
  })

  observe({
    toggle('ref_name_container', condition = length(species_refs()) > 0)
  })

  # change UI of exclude toggle
  observe({
    toggleClass(id = 'toggle_icon', 'fa-plus text-success', condition = is_include())
    toggleClass(id = 'toggle_icon', 'fa-minus text-warning', condition = !is_include())
  })


  # run integration
  subsets <- reactiveValues()
  psubsets <- reactiveValues()

  new_dataset_name <- reactive({
    paste(founder(), input$subset_name, sep = '_')
  })

  founder <- reactive({
    from_dataset <- selected_dataset()
    get_founder(sc_dir, from_dataset)
  })

  subset_clusters <- reactive({
    intersect(cluster_choices()$value, input$subset_features)
  })

  subset_metrics <-  reactive({
    intersect(metric_choices()$value, input$subset_features)
  })
  observeEvent(input$submit_subset, {
    if (disabled_demo) return(NULL)


    error_msg <- validate_subset(selected_dataset(),
                                 input$subset_name,
                                 input$subset_features,
                                 is_include(),
                                 hvgs())

    if (is.null(error_msg)) {
      removeClass('name-container', 'has-error')
      showModal(
        confirmSubsetModal(session,
                           new_dataset_name(),
                           input$ref_name,
                           subset_metrics(),
                           subset_clusters(),
                           is_include())
      )

    } else {
      # show error message
      html('error_msg', html = error_msg)
      addClass('name-container', class = 'has-error')
    }
  })

  observeEvent(input$confirm_subset, {

    removeModal()

    is_include <- is_include()
    from_dataset <- selected_dataset()
    founder <- founder()
    dataset_name <- new_dataset_name()
    is_integrated <- is_integrated()
    hvgs <- hvgs()
    subset_metrics <- subset_metrics()

    ref_name <- input$ref_name
    if (!isTruthy(ref_name)) ref_name <- NULL


    exclude_clusters <- subset_clusters()
    if (is_include && length(exclude_clusters)) {
      exclude_clusters <- setdiff(cluster_choices()$value, exclude_clusters)
    }

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
        ref_name = ref_name,
        tx2gene_dir = tx2gene_dir
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

    # source of selectize.min.js javascript error seems to be in contrastOptions()
    updateSelectizeInput(session, 'subset_features',
                         choices = exclude_choices(),
                         selected = isolate(input$subset_features),
                         options = contrastOptions(),
                         server = TRUE)
  })

  return(new_dataset)
}


#' Logic for integration form toggled by showIntegration
#'
#' @keywords internal
#' @noRd
integrationForm <- function(input, output, session, sc_dir, tx2gene_dir, datasets, integrate_sc, selected_dataset) {
  type <- name <- NULL

  integration_inputs <- c('integration_datasets',
                          'integration_name',
                          'submit_integration',
                          'integration_types',
                          'ref_name')


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
  observeEvent(integrate_sc(), {
    showModal(integrationModal(session, choices = integration_choices()))
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
  observe(toggle(id = 'integration_options_container', condition = allow_integration()))

  # set references based on species
  species <- reactive({
    get_integration_species(input$integration_datasets, sc_dir)
  })

  species_refs <- reactive({
    species <- species()
    get_refs_list(species)
  })

  enable_ref <- reactive(length(species_refs()) > 0)

  observe(toggle('ref_name', condition = enable_ref()))
  observe({
    disabledChoices <- NULL
    if (!enable_ref()) disabledChoices <- 'reference'

    updateCheckboxGroupButtons(session,
                               'integration_types',
                               disabledChoices = disabledChoices)
  })

  observe({
    updateSelectizeInput(session, 'ref_name', choices = c('', species_refs()))
  })


  # run integration
  pintegs <- reactiveValues()
  integs <- reactiveValues()

  use_reference <- reactive({
    'reference' %in% input$integration_types &&
      enable_ref() &&
      allow_integration()
  })

  observe({
    toggle('ref_name_container', condition = use_reference())
  })

  observeEvent(input$submit_integration, {

    dataset_names <- input$integration_datasets
    types <- input$integration_types
    name <- input$integration_name

    ref_name <- input$ref_name
    if (!use_reference()) ref_name <- NULL

    error_msg <- validate_integration(types, name, ref_name, dataset_names, sc_dir)
    if (is.null(error_msg)) {
      removeClass('name-container', 'has-error')
      removeModal()

      integs[[name]] <- callr::r_bg(
        func = run_integrate_saved_scseqs,
        package = 'dseqr',
        args = list(
          sc_dir = sc_dir,
          tx2gene_dir = tx2gene_dir,
          dataset_names = dataset_names,
          integration_name = name,
          integration_types = types,
          ref_name = ref_name
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
comparisonType <- function(input, output, session, is_integrated) {

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
clusterComparison <- function(input, output, session, sc_dir, set_readonly, dataset_dir, dataset_name, resoln_dir, resoln, scseq, annot, annot_path, ref_preds, clusters, dgclogs) {
  cluster_inputs <- c('selected_cluster', 'rename_cluster', 'show_contrasts', 'show_rename')

  disabled_demo <- getShinyOption('is_example', FALSE)
  observe(if (disabled_demo) addClass('show_rename', 'disabled fa-disabled'))

  contrast_options <- reactive({
    on_init <- NULL
    if (set_readonly()) on_init <- disableMobileKeyboard(session$ns('selected_cluster'))

    list(render = I('{option: contrastOptions, item: contrastItem}'),
         onInitialize = on_init)
  })

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
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)

    if (show_contrasts()) {
      test <- isolate(test_cluster())
      choices <- get_contrast_choices(clusters, test)

    } else {
      choices <- get_cluster_choices(clusters, with_all = TRUE, scseq = scseq)
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


  observe({
    sel <- input$selected_cluster
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
    if (disabled_demo) return(NULL)
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


  prev_new_name <- reactiveVal('')
  observeEvent(input$new_cluster_name, {

    curr_new <- input$new_cluster_name
    prev_new <- prev_new_name()

    if (nchar(curr_new) && !nchar(prev_new)) {
      updateActionButton(session, "rename_cluster", icon = tags$i(class ='far fa-fw fa-check-square'))

    } else if (!nchar(curr_new) && nchar(prev_new)) {
      updateActionButton(session, "rename_cluster", icon = tags$i(class ='far fa-fw fa-window-close'))
    }

    prev_new_name(curr_new)
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
                         options = contrast_options(), server = TRUE)
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
      levs <- levels(scseq$cluster)
      if (length(levs) < 2) return(NULL)

      dataset_name <- dataset_name()

      # need dgclogs for presto
      logcounts(scseq) <- dgclogs()

      disableAll(cluster_inputs)
      progress <- Progress$new(session, min = 0, max = 3)
      progress$set(message = "Getting markers", value = 1)
      on.exit(progress$close())

      levels(scseq$cluster) <- seq_along(levs)

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
    sel <- selected_cluster()

    if (!is.null(sel)) {
      new <- markers()[sel]

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
selectedGene <- function(input, output, session, dataset_name, resoln_name, resoln_dir, tx2gene_dir, scseq, h5logs, is_integrated, selected_markers, selected_cluster, type, can_statistic = function()FALSE, cluster_markers = function()NULL, qc_metrics = function()NULL, ambient = function()NULL) {

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
  observe(toggleState('show_pbulk', condition = can_statistic()))
  observe(toggleState('show_biogps', condition = have_biogps()))

  saved_metrics <- reactiveVal()
  custom_metrics <- reactiveVal()
  exist_metric_names <- reactive(c(row.names(scseq()),
                                   colnames(scseq()@colData),
                                   colnames(saved_metrics())
  ))


  have_metric <- reactive(input$custom_metric %in% colnames(custom_metrics()))
  allow_save <- reactive(try(!input$custom_metric %in% row.names(scseq())))

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

    if (methods::is(res, 'DFrame')) {

      prev <- custom_metrics()
      if (!is.null(prev) && nrow(prev) != nrow(res)) {
        custom_metrics(NULL)
        return(NULL)
      }

      if (!is.null(prev)) {
        res <- cbind(prev, res)
        row.names(res) <- colnames(scseq)
      }
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

    if (metric %in% exist_metric_names()) {
      return(NULL)
    }

    # can remove custom metric by selecting as feature and saving empty
    if (!isTruthy(metric)) {
      feature <- selected_gene()
      req(feature %in% colnames(prev))
      res <- prev[, !colnames(prev) %in% feature, drop = FALSE]
      if (ncol(res) == 0) res <- NULL

    } else {
      res <- custom_metrics[, metric, drop = FALSE]
      if (!is.null(prev)) res <- cbind.safe(prev, res)
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

  # update marker genes based on cluster selection
  gene_table <- reactiveVal()

  observe({

    markers <- selected_markers()
    markers_cluster <- names(markers)
    selected_cluster <- selected_cluster()

    markers <- markers[[1]]
    qc_metrics <- qc_metrics()

    # will error if labels
    # also prevents intermediate redraws
    if (is.null(markers) & isTruthy(selected_cluster)) return(NULL)

    # prevent redraws on cluster comparison
    if (isTruthy(selected_cluster) && selected_cluster != markers_cluster) return(NULL)

    qc_first <- FALSE
    if (is.null(markers) | !isTruthy(selected_cluster)) {
      markers <- scseq_genes()
      qc_first <- TRUE
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
    if (qc_first) tables <- tables[c(2,1,3)]

    res <- data.table::rbindlist(tables, fill = TRUE)[, ..cols]

    prev <- isolate(gene_table())
    if (!identical(res, prev)) gene_table(res)
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

    qc_first <- all(colnames(gene_table) %in% c('Feature', 'feature'))
    if (qc_first) {
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
        )
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

  # fixes bug where changing dataset didn't update genes table until play with scrollbar
  observe({
    gene_table <- gene_table()
    if (is.null(gene_table)) return(NULL)
    DT::replaceData(DTproxy, gene_table, rownames = FALSE)
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

safe_set_annot <- function(scseq, annot) {
  if (is.null(scseq) | is.null(annot)) return(NULL)
  if (length(levels(scseq$cluster)) != length(annot)) return(NULL)
  levels(scseq$cluster) <- annot
  return(scseq)
}

safe_set_clusters <- function(scseq, clusters) {
  if (is.null(scseq)) return(NULL)
  if (is.null(clusters)) return(NULL)
  scseq$cluster <- clusters
  return(scseq)
}


safe_set_meta <- function(scseq, meta, groups) {

  if (!isTruthyAll(scseq, meta, groups)) return(NULL)
  if (!all(row.names(meta) %in% scseq$batch)) return(NULL)

  if (length(groups) != 2) return(NULL)
  # if (sum(meta$group %in% groups) < 3) return(NULL)
  if (length(intersect(groups, meta$group)) < 2) return(NULL)

  attach_meta(scseq, meta = meta, groups = groups)
}



#' Logic for cluster plots
#'
#' @keywords internal
#' @noRd
scClusterPlot <- function(input, output, session, scseq, annot, clusters, dataset_name, is_mobile, clusters_marker_view, grid_abundance, grid_expression_fun, selected_gene, show_pbulk, dataset_dir) {


  show_plot <- reactive(!is.null(scseq()))
  observe(toggleCssClass('cluster_plot_container', class = 'invisible', condition = !show_plot()))

  coords <- reactiveVal()
  observe({
    scseq <- scseq()
    if (is.null(scseq)) {
      coords(NULL)
      return(NULL)
    }

    reds <- SingleCellExperiment::reducedDimNames(scseq)
    reds <- reds[reds %in% c('UMAP', 'TSNE')]
    red <- ifelse('UMAP' %in% reds, 'UMAP', reds[1])

    new <- SingleCellExperiment::reducedDim(scseq, red)
    new <- as.data.frame(new)
    new <- new[keep(), ]
    prev <- isolate(coords())
    if (is.null(prev) || !identical(prev, new)) coords(new)
  })

  labels <- reactive({
    scseq <- scseq()
    annot <- annot()
    clusters <- clusters()

    scseq <- safe_set_clusters(scseq, clusters)
    if (is.null(scseq)) return(NULL)

    # fixes issue where no cluster plot after selecting newly imported dataset
    if (is.null(annot)) {
      dataset_dir <- dataset_dir()
      resoln <- load_resoln(dataset_dir)
      annot <- qread.safe(file.path(dataset_dir, resoln, 'annot.qs'))
    }

    scseq <- safe_set_annot(scseq, annot)
    if (is.null(scseq)) return(NULL)

    levels(scseq$cluster) <- format_violin_annot(annot)
    labels <- unname(scseq$cluster)[keep()]

    return(labels)
  })


  label_repels <- reactive({
    coords <- coords()
    labels <- labels()
    if (!isTruthyAll(coords, labels)) return(NULL)
    if (nrow(coords) != length(labels)) return(NULL)

    levels(labels) <- stringr::str_trunc(levels(labels), width = 25, side = 'center')
    labels <- as.character(labels)

    # show nums if too many labels/mobile
    label_coords <- get_label_coords(coords, labels)

    if (is_mobile() | nrow(label_coords) > 45)
      label_coords$label <- gsub('^(\\d+):.+?$', '\\1', label_coords$label)

    nlab <- nrow(label_coords)
    fontsize <- ifelse(nlab > 15, 14, 18)

    label_repels <- repel::repel_text(
      label_coords,
      xrange = range(coords[,1]),
      yrange = range(coords[,2]),
      mar = rep(0, 4),
      box.padding = 0.3,
      fontsize = fontsize,
      direction = 'y')


    return(label_repels)
  })

  label_coords <- reactive({
    label_repels <- label_repels()
    coords <- coords()
    show_grid <- show_grid()

    if (!isTruthyAll(label_repels, coords)) return(NULL)

    label_size <- ifelse(nrow(label_repels) > 15, 14, 18)

    label_repels$anchor <- 'middle'
    label_repels$baseline <- 'center'
    label_repels$size <- label_size

    # hide cluster labels if showing grid
    if (show_grid)
      label_repels$label <- ''

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


  colors <- reactive({
    labels <- labels()
    if (is.null(labels)) return(NULL)

    annot <- levels(labels)
    pal <- get_palette(annot, with_all = TRUE)
    colors <- pal[as.numeric(labels)]

    return(colors)
  })

  keep <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)

    set.seed(0)
    ncells <- ncol(scseq)
    if (ncells > const$max.cells) {
      keep <- sample(ncells, const$max.cells)
    } else {
      keep <- seq_len(ncells)
    }
    return(keep)
  })



  # don't understand magic but mostly stops intermediate color/label change
  # when dataset changes

  update_colors_proxy <- reactiveVal(FALSE)
  update_label_coords_proxy <- reactiveVal(FALSE)
  rendered_dataset <- reactiveVal()

  observeEvent(dataset_name(), {
    update_colors_proxy(FALSE)
    update_label_coords_proxy(FALSE)
  }, priority = 100)

  grid_legend_items = list(
    list(color = '#FF0000', label = '↑'),
    list(color = '#0000FF', label = '↓'),
    list(color = '#989898', label = 'p < .05'),
    list(color = '#EAEAEA', label = 'p ≥ .05')
  )

  output$cluster_plot <- picker::renderPicker({

    coords <- coords()
    if (!isTruthy(coords)) return(NULL)

    is_mobile <- isolate(is_mobile())
    ncells <- nrow(coords)

    scatter_props <- get_scatter_props(is_mobile, ncells)
    deck_props <- list()
    if (is_mobile) {
      deck_props <- list(
        '_typedArrayManagerProps' = list(overAlloc = 1, poolSize = 0)
      )
    }

    colors <- isolate(colors())
    if (is.null(colors)) return(NULL)

    labels <- isolate(labels())
    if (is.null(labels)) return(NULL)

    rendered_dataset(isolate(dataset_name()))

    picker::picker(coords,
                   colors = colors,
                   labels = labels,
                   label_coords = isolate(label_coords()),
                   polygons = isolate(polygons()),
                   point_color_polygons = "white",
                   show_controls = FALSE,
                   grid_legend_items = grid_legend_items,
                   deck_props = deck_props,
                   text_props = text_props,
                   scatter_props = scatter_props)
  })

  proxy <- picker::picker_proxy('cluster_plot')
  observe(picker::update_picker(proxy, clusters_marker_view()))
  observe(picker::update_picker(proxy, labels = labels()))

  observe({
    have <- rendered_dataset()
    sel <- dataset_name()
    if (!is.null(have) && !is.null(sel) && have != sel) return(NULL)
    picker::update_picker(proxy, polygons = polygons())
  })

  observe({
    title <- title()
    title <- ifelse(show_grid(), title, '')
    picker::update_picker(proxy, title = title)
  })

  observe({
    label_coords <- label_coords()
    if (!is.null(label_coords)) {
      have <- rendered_dataset()
      if (!is.null(have) && have != dataset_name()) return(NULL)

      allow <- isolate(update_label_coords_proxy())
      if (allow) picker::update_picker(proxy, label_coords = label_coords)
      update_label_coords_proxy(TRUE)
    }
  })

  observe({
    colors <- colors()
    if (!is.null(colors)) {
      allow <- isolate(update_colors_proxy())
      have <- rendered_dataset()
      if (!is.null(have) && have != dataset_name()) return(NULL)

      if (allow) picker::update_picker(proxy, colors = colors)
      update_colors_proxy(TRUE)
    }
  })

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

get_scatter_props <- function(is_mobile, ncells) {
  if (is_mobile) {
    pt.size <- 3
    if (ncells > 10000) pt.size <- 2
    if (ncells > 20000) pt.size <- 1
  } else {
    pt.size <- 5
    if (ncells > 10000) pt.size <- 3
  }

  scatter_props <- list(
    radiusMinPixels = pt.size,
    radiusMaxPixels = max(6, pt.size*2),
    stroked = pt.size >= 3
  )

  return(scatter_props)
}


hash_meta <- function(meta, groups) {
  meta <- meta[meta$group %in% groups, ]
  tohash <- list(meta = meta, groups = groups)
  digest::digest(tohash, algo = 'murmur32')
}

# get grid abundance data
scGridAbundance <- function(input, output, session, scseq, sc_dir, groups, dplots_dir, meta) {

  grid_abundance <- reactive({
    groups <- groups()
    meta <- meta()

    scseq <- safe_set_meta(scseq(), meta, groups)
    if (is.null(scseq)) return(NULL)
    if (sum(meta$group %in% groups) < 3) return(NULL)

    scseq <- subset_contrast(scseq)

    # add hash for if change groups
    meta_hash <- hash_meta(meta, groups)

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
scMarkerPlot <- function(input, output, session, scseq, annot, clusters, selected_feature, h5logs, clusters_view, is_mobile, meta = function()NULL, groups = function()NULL, show_plot = function()TRUE, markers_view = function()NULL, group = NULL, show_controls = FALSE, deck_props = NULL, added_metrics = function()NULL) {


  observe(toggleCssClass('marker_plot_container', 'invisible', condition = !(show_plot() && have_colors())))
  have_colors <- reactive(length(colors()))

  coords <- reactive({
    scseq <- scseq()
    if (!isTruthy(scseq)) return(NULL)

    reds <- SingleCellExperiment::reducedDimNames(scseq)
    reds <- reds[reds %in% c('UMAP', 'TSNE')]
    red <- ifelse('UMAP' %in% reds, 'UMAP', reds[1])

    coords <- SingleCellExperiment::reducedDim(scseq, red)
    data.frame(coords)
  })

  labels <- reactive({
    scseq <- scseq()
    annot <- annot()
    clusters <- clusters()
    cells <- cells()
    if (is.null(cells)) return(NULL)

    scseq <- safe_set_clusters(scseq, clusters)
    scseq <- safe_set_annot(scseq, annot)
    if (is.null(scseq)) return(NULL)

    cell.idx <- match(cells, colnames(scseq))
    scseq <- scseq[, cell.idx]

    labels <- unname(scseq$cluster)
    levels(labels) <- format_violin_annot(annot)
    return(labels)
  })


  output$marker_plot <- picker::renderPicker({
    if (!show_plot()) return(NULL)

    scseq <- scseq()
    coords <- coords()
    cells <- cells()

    if (!isTruthyAll(cells, scseq, coords)) return(NULL)
    if (!all(cells %in% colnames(scseq))) return(NULL)

    # use xrange and yrange for all data
    xrange <- range(coords[,1])
    yrange <- range(coords[,2])

    scatter_props <- get_scatter_props(is_mobile(), nrow(coords))

    # now subset
    cell.idx <- match(cells, colnames(scseq))
    coords <- coords[cell.idx, ]

    # show group name as plot label
    label_coords <- NULL

    deck_props <- list()
    if (is_mobile()) {
      deck_props <- list(
        '_pickable' = FALSE,
        '_typedArrayManagerProps' = list(overAlloc = 1, poolSize = 0)
      )
    }

    picker::picker(coords,
                   colors = isolate(colors()),
                   labels = isolate(labels()),
                   title = isolate(title()),
                   xrange = xrange,
                   yrange = yrange,
                   show_controls = FALSE,
                   label_coords = label_coords,
                   deck_props = deck_props,
                   scatter_props = scatter_props)
  })


  # column data (custom metric + stored)
  cdata <- reactive({
    scseq <- scseq()
    if (!is.null(group)) scseq <- safe_set_meta(scseq, meta(), groups())
    if (!isTruthy(scseq)) return(NULL)

    metrics <- added_metrics()
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

  samples <- reactive(unique(scseq()$batch))

  colors <- reactive({
    scseq <- scseq()
    if (!is.null(group)) scseq <- safe_set_meta(scseq, meta(), groups())

    cdata <- cdata()
    feature <- selected_feature()
    if (!isTruthyAll(feature, scseq, cdata)) return(NULL)

    is_gene <- feature %in% row.names(scseq)
    is_feature <- feature %in% colnames(cdata)
    is_sample <- feature %in% samples()
    if (!is_gene && !is_feature && !is_sample) return(NULL)

    if (is_sample) cdata[[feature]] <- scseq$batch == feature

    # get feature
    if (is_gene) {
      ft <- h5logs()[feature, ]

    } else {
      ft <- cdata[[feature]]
      names(ft) <- colnames(scseq)
    }

    # e.g. group or sample columns
    if (is.character(ft) | is.factor(ft)) return(NULL)

    # get group
    ids <- colnames(scseq)
    if (!is.null(group)) ids <- ids[which(scseq$orig.ident == group)]

    # order cell ids if logical
    bool.ft <- is.logical(ft)

    set.seed(0)
    if (bool.ft) ids <- ids[order(ft)]
    else ids <- sample(ids)

    # get title and colors
    ft.ids <- ft[ids]
    all.zero <- all(ft.ids == 0)
    ntot <- length(ft.ids)

    if (bool.ft) {
      # title is info
      ncells <- sum(ft.ids)
      pcells <- round(ncells/ntot*100, 1)
      title(sprintf("%s (%s cells :: %s%%)", feature, ncells, pcells))
    }

    if (bool.ft || all.zero) {
      cols <- const$colors$ft

      colors <- rep(cols[1], ntot)
      colors[ft.ids] <- cols[2]

    } else {

      # scale before subsetting to group
      ft.scaled <- scales::rescale(ft, c(0, 1))
      ft.scaled <- ft.scaled[ids]

      is_qc <- feature %in% const$features$qc

      cols <- if (is_qc) const$colors$qc else const$colors$ft
      colors <- scales::seq_gradient_pal(cols[1], cols[2])(ft.scaled)

      # title is group
      prev <- isolate(title())
      if (prev != group_title) title(group_title)
    }

    # down-sample after getting titles
    ncells <- ncol(scseq)
    if (ncells > const$max.cells) {
      set.seed(0)
      keep <- sample(ncells, const$max.cells)
      cells_keep <- colnames(scseq)[keep]
      is.keep <- ids %in% cells_keep
      ids <- ids[is.keep]
      colors <- colors[is.keep]
    }

    # update ids
    prev <- isolate(cells())
    changed.ids <- !identical(prev, ids)
    if (changed.ids) cells(ids)

    update_colors_proxy(!changed.ids)

    return(colors)
  })


  proxy <- picker::picker_proxy('marker_plot')
  observe(picker::update_picker(proxy, clusters_view()))
  observe(picker::update_picker(proxy, markers_view()))
  observe(picker::update_picker(proxy, labels = labels()))
  observe(picker::update_picker(proxy, title = title()))

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
scViolinPlot <- function(input, output, session, selected_gene, selected_cluster, scseq, annot, clusters, plots_dir, h5logs, is_mobile) {

  show_plot <- reactive(!is.null(plot()))
  observe(toggle('violin_plot', condition = show_plot()))


  violin_data <- reactive({
    gene <- selected_gene()
    cluster <- selected_cluster()
    if (!isTruthy(gene)) return(NULL)

    scseq <- scseq()
    annot <- annot()
    clusters <- clusters()
    scseq <- safe_set_clusters(scseq, clusters)
    scseq <- safe_set_annot(scseq, annot)
    if (is.null(scseq)) return(NULL)

    is.gene <- gene %in% row.names(scseq)
    is.num <- is.gene || is.numeric(scseq@colData[[gene]])
    if (!is.num) return(NULL)

    h5logs <- if (is.gene) h5logs() else NULL
    if (is.null(scseq)) return(NULL)
    vdat <- get_violin_data(gene, scseq, cluster, with_all = TRUE, h5logs=h5logs)

    return(vdat)
  })

  height <- reactiveVal()

  plot <- reactive({
    violin_data <- violin_data()
    if (is.null(violin_data)) return(NULL)
    if (all(violin_data$df$x == 0)) return(NULL)

    height(max(length(levels(violin_data$df$y))*38, 420))
    plot_violin(violin_data = violin_data, is_mobile = is_mobile())
  })


  output$violin_plot <- renderPlot(plot(), height=height)
}



# modal to upload single-cell dataset
importSingleCellModal <- function(session, show_init) {

  modalDialog(
    tags$div(
      class='alert alert-warning', role = 'alert',
      tags$div(tags$b("For each sample upload files:")),
      tags$br(),
      tags$div("- ", tags$code("filtered_feature_bc_matrix.h5"), "or"),
      tags$br(),
      tags$div("- ", tags$code("matrix.mtx"), ", ", tags$code("barcodes.tsv"), ", and", tags$code("features.tsv"), 'or'),
      tags$br(),
      tags$div("- ", tags$code("*.rds"), "or", tags$code("*.qs"), "with", tags$code("Seurat"), "or", tags$code("SingleCellExperiment"), "objects"),
      hr(),
      '🌱 Add prefixes e.g.', tags$i(tags$b('sample_matrix.mtx')), ' to auto-name samples:',
      tags$a(href = 'https://dseqr.s3.amazonaws.com/GSM3972011_involved.zip', target = '_blank', 'example files.')
    ),
    div(class='dashed-upload',
        fileInput(placeholder = 'drag files here',
                  session$ns('up_raw'),
                  label = '',
                  buttonLabel = 'Upload',
                  width = '100%',
                  accept = c('.h5', '.hdf5', '.tsv', '.fastq.gz', '.mtx', '.rds', '.qs'),
                  multiple = TRUE
        )
    ),
    tags$div(
      id = session$ns('sample_name_container'),
      style = ifelse(show_init, '', 'display: none;'),
      hr(),
      shinypanel::textInputWithButtons(
        session$ns('sample_name'),
        'Sample name for selected rows:',
        actionButton(session$ns('add_sample'), '', icon('plus', class='fa-fw')),
        container_id = session$ns('validate-up'),
        help_id = session$ns('error_msg')
      ),
      hr()
    ),
    div(
      id = session$ns('up_table_container'),
      class= ifelse(show_init, 'dt-container', 'invisible-height dt-container'),
      DT::dataTableOutput(session$ns('up_table'), width = '100%'),
    ),
    title = 'Upload Single Cell Datasets',
    size = 'l',
    footer = tagList(
      actionButton(
        inputId = session$ns("import_datasets"),
        label = "Import Datasets",
        class = ifelse(show_init, 'btn-warning', 'btn-warning disabled')
      ),
      tags$div(class='pull-left', modalButton('Cancel'))
    ),
    easyClose = FALSE,
  )
}



get_species_choices <- function(detected_species) {
  is_detected <- !is.null(detected_species)
  species <- unique(ensmap$species)

  # some common species first
  first_name <- c('Homo sapiens', 'Mus musculus', 'Rattus norvegicus', 'Canis lupus familiaris',
                  'Heterocephalus glaber', 'Macaca mulatta', 'Gorilla gorilla gorilla')

  first_idx <- sapply(first_name, function(x) which(species == x))
  other_idx <- setdiff(seq_along(species), first_idx)
  species <- species[c(first_idx, other_idx)]


  if (is_detected) {
    names(species) <- species
    names(species)[species == detected_species] <- paste(detected_species, '(detected)')
  }

  return(species)
}

confirmImportSingleCellModal <- function(session, metric_choices, detected_species, species_refs, warn_robject) {

  # indicate detected species
  is_detected <- !is.null(detected_species)
  have_refs <- !is.null(species_refs)
  species <- get_species_choices(detected_species)

  triangle <- tags$i(class = 'fas fa-exclamation-triangle', style='color: orange;')

  UI <- div(
    selectizeInput(
      session$ns('import_species'),
      'Select species:',
      width = '100%',
      choices = c('', species),
      selected = detected_species,
    ),
    selectizeInput(
      session$ns('qc_metrics'),
      HTML('Select <a href="https://docs.dseqr.com/docs/single-cell/quality-control/" target="_blank">QC</a> metrics:'),
      width = '100%',
      choices = c('all', 'all and none', 'none', metric_choices),
      selected = 'all',
      multiple = TRUE),
    div(
      id = session$ns('ref_name_container'),
      style = ifelse(have_refs, '', 'display: none;'),
      selectizeInput(
        session$ns('ref_name'),
        HTML('Select reference:'),
        width = '100%',
        choices = c('', species_refs),
        options = list(placeholder = 'optional'))
    ),
    div(
      style = ifelse(warn_robject, '', 'display: none;'),
      style='color: grey; font-style: italic;',
      tags$p(triangle, tags$b(' for R object import:')),
      tags$div(' - QC is skipped if multi-sample'),
      tags$div(' - Reference based analyses available after import'))

  )

  modalDialog(
    UI,
    title = 'Import settings',
    size = 'm',
    footer = tagList(
      actionButton(
        session$ns('confirm_import_datasets'),
        'Import',
        class = ifelse(is_detected, 'btn-warning', 'btn-warning disabled')),
      tags$div(class='pull-left', modalButton('Cancel'))
    )
  )
}

