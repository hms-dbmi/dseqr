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
             feature_plot_clusters = scForm$feature_plot_clusters,
             dgrlogs = scForm$dgrlogs)

  callModule(scBioGpsPlot, 'biogps_plot',
             selected_gene = scForm$clusters_gene,
             species = scForm$species)


  callModule(scRidgePlot, 'ridge_plot',
             selected_gene = scForm$clusters_gene,
             selected_cluster = scForm$clusters_cluster,
             scseq = scForm$scseq,
             annot = scForm$annot,
             plots_dir =scForm$plots_dir,
             dgrlogs = scForm$dgrlogs)


  # sample comparison plots ---

  # cluster plot in top right
  callModule(scAbundancePlot, 'abundance_plot',
             scseq = scForm$scseq,
             dataset_dir = scForm$dataset_dir,
             dplots_dir = scForm$dplots_dir,
             comparison_type = scForm$comparison_type,
             compare_groups = scForm$compare_groups,
             sc_dir = sc_dir)

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

  # dgRMatrix logcounts for fast row indexing
  dgrlogs <- reactive({
    dataset_dir <- dataset_dir()
    if (is.null(dataset_dir)) return(NULL)
    qs::qread(file.path(dataset_dir, 'dgrlogs.qs'))
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
  have_cluster <- reactive(isTruthy(scSampleClusters$selected_cluster()))
  observe(toggle(id = "sample_cluster_input",condition = have_contrast()))
  observe(toggle(id = "sample_gene_input",condition = have_cluster()))

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
                               dataset_name = scDataset$dataset_name,
                               show_dprimes = scSampleGene$show_dprimes)

  scSampleClusters <- callModule(scSampleClusters, 'sample_clusters',
                                 scseq = scSampleGroups$scseq,
                                 dgrlogs = dgrlogs,
                                 lm_fit = scSampleGroups$lm_fit,
                                 lm_fit_grid = scSampleGroups$lm_fit_grid,
                                 groups = scSampleGroups$groups,
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
                             selected_markers = scSampleClusters$top_table,
                             selected_cluster = scSampleClusters$selected_cluster,
                             type = 'samples',
                             ambient = scSampleClusters$ambient)



  scLabelsComparison <- callModule(scLabelsComparison, 'labels',
                                   cluster_choices = scSampleClusters$cluster_choices)



  return(list(
    scseq = scseq,
    samples_gene = scSampleGene$selected_gene,
    clusters_gene = scClusterGene$selected_gene,
    custom_metrics = scClusterGene$custom_metrics,
    show_ridge = scClusterGene$show_ridge,
    samples_pfun_left = scSampleClusters$pfun_left,
    samples_pfun_right = scSampleClusters$pfun_right,
    samples_pfun_right_bottom = scSampleClusters$pfun_right_bottom,
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
    feature_plot_clusters = scClusterPlots$feature_plot_clusters,
    cluster_plot = scClusterPlots$cluster_plot,
    annot = annot,
    compare_groups = scSampleGroups$groups,
    dgrlogs = dgrlogs
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
scSampleGroups <- function(input, output, session, dataset_dir, resoln_dir, dataset_name, input_scseq = function()NULL, show_dprimes = function()FALSE) {
  group_options <- list(render = I('{option: bulkContrastOptions, item: bulkContrastItem}'))
  input_ids <- c('click_dl_meta', 'click_up_meta', 'compare_groups')


  # need for drugs tab
  scseq_clust <- reactive({
    scseq <- input_scseq()
    if (!is.null(scseq)) return(scseq)

    dataset_dir <- dataset_dir()
    if (!isTruthy(dataset_dir)) return(NULL)

    scseq <- load_scseq_qs(dataset_dir)
    attach_clusters(scseq, resoln_dir())
  })

  scseq <- reactive({
    scseq <- scseq_clust()
    if (is.null(scseq)) return(NULL)

    meta <- up_meta()
    groups <- input$compare_groups
    if (is.null(meta)) return(scseq)
    if (length(groups) != 2) return(scseq)


    attach_meta(scseq, meta=meta, groups=groups)
  })


  prev_path <- reactive({
    if (is.null(dataset_dir())) return(NULL)
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
    groups <- input$compare_groups
    if (length(groups) != 2) return(NULL)
    if (is.null(prev_path())) return(NULL)

    qs::qsave(groups, prev_path())
  })


  summed_grid_path <- reactive({
    summed_path <- file.path(resoln_dir(), 'summed_grid.qs')

    if (!file.exists(summed_path)){
      scseq <- scseq()
      grid <- get_grid(scseq)
      scseq$cluster <- factor(grid$cluster)
      summed <- aggregate_across_cells(scseq)
      qs::qsave(summed, summed_path)
    }

    return(summed_path)
  })

  lm_fit <- reactive({
    fit_path <- file.path(resoln_dir(), 'lm_fit_0svs.qs')

    if (file.exists(fit_path)) {
      fit <- qs::qread(fit_path)

    } else {
      disableAll(input_ids)
      progress <- Progress$new(session, min = 0, max = 4)
      on.exit(progress$close())

      fit <- run_limma_scseq(meta = up_meta(),
                             fit_path = fit_path,
                             species = species(),
                             summed = summed(),
                             progress = progress)


      enableAll(input_ids)
    }
    return(fit)
  })


  lm_fit_grid <- reactiveVal()

  gfits <- reactiveValues()
  pgfits <- reactiveValues()

  observeEvent(show_dprimes(), {
    if (!show_dprimes()) return(NULL)
    if (!is.null(lm_fit_grid())) return(NULL)

    fit_path <- file.path(dataset_dir(), 'lm_fit_grid_0svs.qs')

    if (file.exists(fit_path)) {
      lm_fit_grid(qs::qread(fit_path))

    } else {
      dataset_name <- dataset_name()
      if (dataset_name %in% names(gfits)) return(NULL)

      disableAll(input_ids)
      progress <- Progress$new(session, min = 0, max = 5)
      progress$set(message = "Comparing groups:", detail = "pseudobulk", value = 1)
      summed_path <- summed_grid_path()
      enableAll(input_ids)

      browser()

      run_limma_scseq(
        meta = up_meta(),
        fit_path = fit_path,
        species = species(),
        summed_path = summed_path,
        value = 1,
        min.count = 1,
        min.total.count = 3
      )

      gfits[[dataset_name]] <- callr::r_bg(
        func = run_limma_scseq,
        package = 'dseqr',
        args = list(
          meta = up_meta(),
          fit_path = fit_path,
          species = species(),
          summed_path = summed_path,
          value = 1,
          min.count = 1,
          min.total.count = 3
        )
      )

      pgfits[[dataset_name]] <- progress
    }
  })

  new_gfit <- reactiveVal(NULL)

  observe({
    invalidateLater(5000, session)
    handle_sc_progress(gfits, pgfits, new_gfit)
  })

  observeEvent(new_gfit(), {
    gfit_name <- new_gfit()
  })

  gfit_running <- reactive({
    invalidateLater(5000, session)
    dataset_name() %in% names(gfits)
  })

  # reset lm_fit_grid when change dataset name
  observe({
    dataset_name()
    lm_fit_grid(NULL)
  })




  groups <- reactiveVal()
  observeEvent(dataset_name(), groups(NULL))
  observe(groups(input$compare_groups))

  return(list(
    lm_fit = lm_fit,
    gfit_running = gfit_running,
    lm_fit_grid = lm_fit_grid,
    groups = groups,
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
scSampleClusters <- function(input, output, session, scseq, lm_fit, lm_fit_grid, groups, dataset_dir, resoln_dir, plots_dir, feature_plot, dataset_name, sc_dir, input_annot = function()NULL, show_dprimes = function()TRUE, is_integrated = function()TRUE, is_sc = function()TRUE, exclude_ambient = function()FALSE, comparison_type = function()'samples', applied = function()TRUE, is_mobile = function()FALSE, dgrlogs = function()NULL) {
  cluster_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))
  input_ids <- c('click_dl_anal', 'selected_cluster')


  contrast_dir <- reactive({
    groups <- groups()
    if (length(groups) != 2) return(NULL)
    contrast <- paste0(groups, collapse = '_vs_')
    file.path(isolate(resoln_dir()), contrast)
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
  clusters_str <- reactive(collapse_sorted(input$selected_cluster))
  drug_paths <- reactive(get_drug_paths(contrast_dir(), clusters_str()))
  go_path <- reactive(file.path(contrast_dir(), paste0('go_', clusters_str(), '.qs')))
  kegg_path <- reactive(file.path(contrast_dir(), paste0('kegg_', clusters_str(), '.qs')))
  goana_path <- reactive(file.path(contrast_dir(), paste0('goana_', clusters_str(), '.qs')))
  kegga_path <- reactive(file.path(contrast_dir(), paste0('kegga_', clusters_str(), '.qs')))
  top_tables_paths <- reactive(file.path(contrast_dir(), 'top_tables.qs'))


  # plot functions for left
  sel <- reactive(input$selected_cluster)


  pfun_left <- reactive({
    req(is_integrated())

    function(gene) {
      if(!isTruthy(gene)) return(NULL)

      if (show_dprimes() & is_integrated()) {
        top_tables <- top_tables_grid()
        if(is.null(top_tables)) return(NULL)

        pname <- paste(gene, 'diff_grid.qs', sep='-')
        plots_dir <- plots_dir()
        if (is.null(plots_dir)) return(NULL)
        plot_path <- file.path(plots_dir, pname)

        scseq <- scseq()
        if (file.exists(plot_path)) {
          plot_data <- qs::qread(plot_path)

        } else {
          grid <- get_grid(scseq)
          plot_data <- get_gene_diff(gene, top_tables, grid)
          qs::qsave(plot_data, plot_path)
        }

        plot <- plot_scseq_diff(plot_data, gene)
        pfun <- list(plot = plot, height = 453)

      } else {
        scseq <- scseq()
        if(is.null(scseq)) return(NULL)
        gene_data <- fast_dgr_row(dgrlogs(), gene)

        # update base feature plot
        plot <- feature_plot()
        plot <- update_feature_plot(plot, gene_data, gene)
        plot <- plot_feature_sample(gene, scseq(), 'test', plot=plot)
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

      gene_data <- fast_dgr_row(dgrlogs(), gene)
      if (is.null(gene_data)) return(NULL)

      # update base feature plot
      plot <- feature_plot()
      plot <- update_feature_plot(plot, gene_data, gene)
      if (!show_dprimes() | !is_integrated()) {
        plot <- plot_feature_sample(gene, scseq(), 'ctrl', plot=plot) +
          ggplot2::theme(legend.position = 'none')
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

      try(ridge_data <- get_ridge_data(
        gene, scseq, sel, by.sample = TRUE, with_all = TRUE, dgrlogs=dgrlogs()))

      plot <- VlnPlot(ridge_data = ridge_data, with.height = TRUE, is_mobile = is_mobile())
      return(plot)
    }
  }) %>% debounce(20)


  summed <- reactive(qs::qread(file.path(resoln_dir(), 'summed.qs')))


  cluster_ambient <- reactive({
    compare <- input$compare_groups
    if (length(compare) != 2) return(NULL)

    sel <- sel()
    resoln_dir <- resoln_dir()
    ambience_path <- file.path(resoln_dir, paste0('ambience_', sel, '.qs'))

    if (file.exists(ambience_path)) {
      amb <- qs::qread(ambience_path)

    } else {
      # no cluster e.g. for 'All Clusters'
      summed <- summed()
      if (!sel %in% levels(summed$cluster)) return(NULL)

      # check that not disabled
      if (cluster_choices()[sel, 'disabled']) return(NULL)

      disableAll(input_ids)
      progress <- Progress$new(session, min = 0, max = 2)
      on.exit(progress$close())
      progress$set(message = paste("Calculating ambience: cluster", sel), value = 1)

      amb <- get_ambience(scseq())
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
      if (is.null(fit)) return(NULL)


      disableAll(input_ids)
      nmax <- length(fit)+1
      progress <- Progress$new(session, min = 0, max = nmax)
      on.exit(progress$close())
      progress$set(message = "Differential Expression:", value = 1)

      tts <- list()
      for (i in seq_along(fit)) {

        cluster <- names(fit)[i]
        progress$set(detail=paste('cluster', cluster), value = i)

        tts[[cluster]] <- crossmeta::get_top_table(
          fit[[cluster]],
          groups(),
          allow.no.resid = TRUE)
      }

      # TODO: RobustRankAggreg:  - 1v1 comparison
      # TODO: run_esmeta:        - metap to do pvalue meta-analyses

      # add 'All Clusters' result
      progress$set(detail='all clusters', value = nmax)
      annot <-  qs::qread(file.path(resoln_dir, 'annot.qs'))
      all <- as.character(length(annot)+1)
      es <- run_esmeta(tts)

      if (!is.null(es)) {
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
        progress$set(detail=paste('cluster', cluster), value = i)

        tt <- tryCatch(crossmeta::get_top_table(fit[[cluster]]),
                       error = function(e) return(NULL))
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

    ambient <- cluster_ambient()
    tt <- top_tables()[[sel]]
    if (is.null(tt)) return(NULL)

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

    saved_drugs <- any(grepl('^cmap_res_', list.files(contrast_dir)))

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
        paths <- get_drug_paths(contrast_dir(), cluster)
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

  observe(toggleState('click_dl_anal', condition = isTruthy(input$selected_cluster)))

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
    pfun_left = pfun_left,
    pfun_right = pfun_right,
    pfun_right_bottom = pfun_right_bottom
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
      plot <- plot_feature(scseq, row.names(scseq)[1]) +
        ggplot2::labs(x='', y='')
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
      plot <- plot_feature(scseq, row.names(scseq)[1]) +
        ggplot2::labs(x='', y='')
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
    load_scseq_qs(dataset_dir())
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
      is.azi <- file.exists(file.path(dataset_dir(), 'azimuth_ref.qs'))
      if (is.azi) return(NULL)

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
  observeEvent(new_dataset_dir(), showModal(confirmModal(session, 'quant', metric_choices)))


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
    if (!isTruthy(azimuth_ref)) azimuth_ref <- NULL

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


# modal to confirm adding single-cell dataset
confirmModal <- function(session, type = c('quant', 'subset'), metric_choices = NULL) {
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

  height <- reactive({
    h <- res()$height
    if (is.null(h)) return(1)
    else return(h)
  })



  filename <- function() {
    fname <- plot()$labels$title
    fname <- gsub(':', '', fname)
    paste0(fname, '.csv')
  }

  plot <- reactive({
    res <- res()
    if (is.null(res)) return(NULL)
    return(res$plot)
  })


  content <- function(file) {
    d <- plot()$data

    # clean up data for ridgeplots
    if (!'TSNE1' %in% colnames(d)) {
      d <- d[, c('x', 'y')]
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
    if (set >= 0.1 & set <= 5.1 & !first_set()) resoln(set)

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
    axist <- file.exists(apath)
    is_azimuth(axist)

    rname(ifelse(axist, 'resoln_azi', 'resoln'))

    rpath <- file.path(dataset_dir(), 'resoln.qs')
    resoln_path(rpath)
    init <- qread.safe(rpath, 1)
    resoln(init)
    updateNumericInput(session, rname(), value=init)
  }, priority = 1)

  resoln_name <- reactive(file.path(dataset_name(), paste0('snn', resoln())))

  # clusters after change resolution
  clusters_path <- reactive(file.path(resoln_dir(), 'clusters.qs'))

  clusters <- reactive({
    resoln <- resoln()
    clusters_path <- clusters_path()
    clusters <- qread.safe(clusters_path)

    if (!is.null(clusters)) {
      qs::qsave(resoln, resoln_path())

      resoln_name <-  paste0('snn', resoln)
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

  # change UI of exclude toggle
  observe({
    toggleClass(id = 'toggle_icon', 'fa-plus text-success', condition = is_include())
    toggleClass(id = 'toggle_icon', 'fa-minus text-warning', condition = !is_include())
  })


  # run integration
  subsets <- reactiveValues()
  psubsets <- reactiveValues()

  observeEvent(input$submit_subset, {
    showModal(confirmModal(session, 'subset'))
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

    choices <- ds %>%
      dplyr::group_by(type) %>%
      dplyr::summarise(names = list(name))

    names(choices$names) <- choices$type
    choices$names
  })

  new_dataset <- reactiveVal()


  is_include <- reactive({ input$toggle_exclude %% 2 != 0 })
  allow_integration <- reactive(length(input$integration_datasets > 1))

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
  observe(toggle(id = 'name-container', condition = allow_integration()))


  # run integration
  pintegs <- reactiveValues()
  integs <- reactiveValues()

  use_azimuth <- reactive('Azimuth' %in% input$integration_types)

  observe({
    toggle('azimuth_ref_container', condition = use_azimuth())
  })

  observeEvent(input$submit_integration, {

    dataset_names <- input$integration_datasets
    types <- input$integration_types
    name <- input$integration_name

    azimuth_ref <- input$azimuth_ref
    if (!use_azimuth()) azimuth_ref <- NULL

    error_msg <- validate_integration(types, name, azimuth_ref)
    if (is.null(error_msg)) {
      removeClass('name-container', 'has-error')

      integs[[integration_name]] <- callr::r_bg(
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

      disableAll(cluster_inputs)
      progress <- Progress$new(session, min = 0, max = 3)
      progress$set(message = "Getting markers", value = 1)
      on.exit(progress$close())

      levels(scseq$cluster) <- seq_along(levels(scseq$cluster))
      markers <- get_presto_markers(scseq)
      resoln_name <- paste0(dataset_name(), '/', 'snn', resoln())

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
      if (is.null(prev) || !identical(row.names(new), row.names(prev)))
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
selectedGene <- function(input, output, session, dataset_name, resoln_name, resoln_dir, scseq, is_integrated, selected_markers, selected_cluster, type, cluster_markers = function()NULL, qc_metrics = function()NULL, ambient = function()NULL) {
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
    gt <- gene_table()
    if (is.null(gt) | is.null(row)) return('')
    gt[row]$feature
  })


  # toggle for ridgeline
  gene_selected <- reactive({
    sel <- feature()
    scseq <- scseq()
    if (!isTruthyAll(sel, scseq)) return(FALSE)
    sel %in% row.names(scseq)
  })

  have_biogps <- reactive({
    feature() %in% biogps[, SYMBOL]
  })

  sel_ridge <- reactive(input$show_ridge %% 2 != 1)
  show_ridge <- reactive(sel_ridge() | !have_biogps())
  observe(toggleClass(id = "show_ridge", 'btn-primary', condition = !sel_ridge()))

  # toggle for showing custom metric
  show_custom_metric <- reactive(type != 'samples' && (input$show_custom_metric %%2 != 0))

  observe({
    toggle('custom_metric_panel', anim = TRUE, condition = show_custom_metric())
    if (show_custom_metric() & have_metric()) selected_gene(input$custom_metric)
  })

  observe(if (!show_custom_metric()) selected_gene(isolate(feature())))
  observe(toggleClass('show_custom_metric', class = 'btn-primary', condition = show_custom_metric()))

  # toggle for showing dprimes plot
  # TODO: show spinner and disable if calculating for dataset_name
  show_dprimes <- reactive(type == 'samples' && (input$show_dprimes %%2 != 0))
  observe(toggleClass(id = "show_dprimes", 'btn-primary', condition = show_dprimes()))


  # disable buttons when not valid
  observe({
    toggleState('show_dprimes', condition = gene_selected() & is_integrated())
    toggleState('show_ridge', condition = have_biogps())
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


  scseq_genes <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)

    markers <- cluster_markers()
    if (is.null(markers)) return(NULL)
    construct_top_markers(markers, scseq)
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
  tx2gene <- reactive(dseqr.data::load_tx2gene(species()))

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

    if (is.null(markers) || is.null(scseq())) return(NULL)

    tables <- list()
    tables[[1]] <- get_gene_table(markers, species = species(), tx2gene = tx2gene())
    tables[[2]] <- get_qc_table(qc_metrics)

    cols <- colnames(tables[[1]])
    if (qc_first()) tables <- tables[c(2,1)]

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

    vis_targ <- (length(cols)-1)
    search_targs <- 0:(vis_targ-1)

    # prevent sort when qc_first
    sort_targs <- 0
    if (qc_first()) sort_targs <- '_all'

    dt <- DT::datatable(
      gene_table,
      class = 'cell-border',
      rownames = FALSE,
      escape = FALSE, # to allow HTML in table
      selection = 'single',
      extensions = 'Scroller',
      options = list(
        deferRender = TRUE,
        scroller = TRUE,
        dom = '<"hidden"f>t',
        bInfo = 0,
        scrollY=250,
        search = list(regex = TRUE),
        language = list(search = 'Select feature to plot:'),
        columnDefs = list(
          list(visible = FALSE, targets = vis_targ),
          list(searchable = FALSE, targets = search_targs),
          list(sortable = FALSE, targets = sort_targs)
        ),
        headerCallback = DT::JS("function(thead, data, start, end, display) {$('th', thead).css('font-weight', 'normal');}")
      )
    )

    if (length(pct_targs)) dt <- DT::formatRound(dt, pct_targs, digits = 0)
    if (length(frac_targs)) dt <- DT::formatRound(dt, frac_targs, digits = 2)
    return(dt)

  }, server = TRUE)

  DTproxy <- DT::dataTableProxy("gene_table")

  observeEvent(input$gene_search, {
    DT::updateSearch(DTproxy, keywords = list('global' = input$gene_search))
  })

  observe({
    toggle('gene_search_input', condition = !is.null(gene_table()))
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

    omit_labels <- is_mobile() || length(annot) > 30
    if (omit_labels) annot <- as.character(seq_along(annot))


    if (isTruthy(cluster) && nclus >= as.numeric(cluster))
      hl <- as.numeric(cluster)

    pl <- update_cluster_plot(plot, annot, hl)

    if (omit_labels) {
      pl <- pl +
        ggplot2::ggtitle('Labels Omitted') +
        ggplot2::theme(
          plot.title.position = "plot",
          plot.title = ggplot2::element_text(
            color = 'lightgray', hjust = 1, size = 16, face = 'plain', margin = ggplot2::margin(b = 15)))
    }
    return(pl)
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


scAbundancePlot <- function(input, output, session, scseq, dataset_dir, sc_dir, comparison_type, compare_groups, dplots_dir) {

  show_plot <- reactive(length(compare_groups()) == 2 & comparison_type() == 'samples')
  observe(toggle('abundance_plot', condition = show_plot()))

  plot_data <- reactive({
    dataset_dir <- dataset_dir()
    scseq <- scseq()
    if (is.null(scseq)) return(NULL)

    groups <- compare_groups()
    if (length(groups) != 2) return(NULL)

    scseq_groups <- unique(scseq$orig.ident)
    if (length(scseq_groups) != 2) return(NULL)

    apath <- file.path(dplots_dir(), 'abundance_plot_data.qs')
    mpath <- file.path(dataset_dir, 'meta.qs')
    have <- file.exists(apath)
    valid <- file.info(apath)$ctime > file.info(mpath)$ctime

    if (have && valid) {
      plot_data <- qs::qread(apath)

    } else {
      plot_data <- get_abundance_diff(scseq)
      qs::qsave(plot_data, apath)
    }

    return(plot_data)
  })

  plot <- reactive({
    plot_data <- plot_data()
    if (is.null(plot_data)) return(NULL)
    plot_scseq_diff(plot_data, feature = 'abundance')
  })


  content <- function(file) {
    data <- plot()$data
    utils::write.csv(data, file)
  }


  output$abundance_plot <- shiny::renderPlot(plot())

  return(list(
    plot = plot
  ))
}


#' Logic for marker feature plots
#'
#' @keywords internal
#' @noRd
scMarkerPlot <- function(input, output, session, scseq, selected_feature, dataset_name, plots_dir, feature_plot_clusters, custom_metrics = function()NULL, dgrlogs = function()NULL) {


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

    is_bool <- is.logical(cdata[[feature]])
    is_num <- is_gene || is.numeric(cdata[[feature]])


    if (is_num) {
      plots_dir <- plots_dir()
      plot <- feature_plot_clusters()
      is.gene <- feature %in% row.names(scseq)
      dgrlogs <- dgrlogs()
      have.dgrlogs <- !is.null(dgrlogs)

      if (is.gene & have.dgrlogs) fdata <- fast_dgr_row(dgrlogs, feature)
      else if (is.gene) fdata <- SingleCellExperiment::logcounts(scseq)[feature, ]
      else fdata <- scseq[[feature]]

      pl <- update_feature_plot(plot, fdata, feature)

    } else if (is_bool) {

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
    if (!length(species)) return(NULL)
    if (species == 'Mus musculus') gene <- toupper(gene)
    if (!gene %in% biogps[, SYMBOL]) return(NULL)

    plot_biogps(gene)
  })
}


#' Logic for Ridge plot for clusters
#'
#' @keywords internal
#' @noRd
scRidgePlot <- function(input, output, session, selected_gene, selected_cluster, scseq, annot, plots_dir, dgrlogs) {

  height <- reactive({
    scseq <- scseq()
    if (is.null(scseq)) height <- 453
    else height <- max(length(levels(scseq$cluster))*38, 420)
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

    scseq <- scseq()
    if (is.null(scseq)) return(NULL)
    rdat <- get_ridge_data(gene, scseq, cluster, with_all = TRUE, dgrlogs=dgrlogs())

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

