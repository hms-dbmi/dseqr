#' Logic for Single Cell Exploration page
#' @export
#' @keywords internal
scPage <- function(input, output, session, sc_dir, new_dataset) {

  # the analysis and options
  scForm <- callModule(scForm, 'form',
                       new_dataset = new_dataset,
                       sc_dir = sc_dir)


  scCluster <- callModule(scClusterPlot, 'cluster_plot',
                          scseq = scForm$scseq,
                          plot_styles = scForm$plot_styles,
                          selected_group = NULL,
                          cached_plot = reactive(NULL))

  # showing cluster comparison ----
  scMarkerCluster <- callModule(scMarkerPlot, 'marker_plot_cluster',
                                scseq = scForm$scseq,
                                selected_gene = scForm$selected_gene_cluster,
                                plot_styles = scForm$plot_styles,
                                selected_group = 'all',
                                cached_plot = reactive(NULL))

  callModule(scBioGpsPlot, 'biogps_plot',
             selected_gene = scForm$selected_gene_cluster)

  # showing sample comparison -----
  scMarkerSample <- callModule(scMarkerPlot, 'marker_plot_test',
                               scseq = scForm$scseq,
                               selected_gene = scForm$selected_gene_sample,
                               plot_styles = scForm$plot_styles,
                               selected_group = 'test',
                               cached_plot = reactive(NULL))

  callModule(scMarkerPlot, 'marker_plot_ctrl',
             scseq = scForm$scseq,
             selected_gene = scForm$selected_gene_sample,
             plot_styles = scForm$plot_styles,
             selected_group = 'ctrl',
             cached_plot = scMarkerSample$plot)

  # showing labels comparison --------

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
             cached_plot = label_plot1,
             plot_styles = scForm$plot_styles)

  callModule(scClusterPlot, 'label_plot2',
             cached_plot = label_plot2,
             plot_styles = scForm$plot_styles)

  observe({
    toggle(id = "sample_comparison_row",  condition = scForm$comparison_type() == 'samples')
    toggle(id = "cluster_comparison_row", condition = scForm$comparison_type() == 'clusters')
    toggle(id = "label_comparison_row", condition = scForm$comparison_type() == 'labels')
  })

  return(NULL)
}

#' Logic for form on Single Cell Exploration page
#' @export
#' @keywords internal
scForm <- function(input, output, session, sc_dir, new_dataset) {

  # updates if new integrated dataset
  new_anal <- reactive({
    scIntegration()
  })


  # the analysis and options
  scAnal <- callModule(selectedAnal, 'anal',
                       sc_dir = sc_dir,
                       new_anal = new_anal,
                       new_dataset = new_dataset)



  # label transfer between datasets
  scLabelTransfer <- callModule(labelTransferForm, 'transfer',
                                sc_dir = sc_dir,
                                anal_options = scAnal$anal_options,
                                show_label_transfer = scAnal$show_label_transfer,
                                selected_anal = scAnal$selected_anal,
                                scseq = scAnal$scseq)

  # dataset integration
  scIntegration <- callModule(integrationForm, 'integration',
                              sc_dir = sc_dir,
                              anal_options = scAnal$anal_options,
                              show_integration = scAnal$show_integration)

  # show original labels for integrated datasets
  integrationAnnotAnals <- callModule(selectedAnnot, 'annot',
                                      scseq = scAnal$scseq,
                                      is.integrated = scAnal$is.integrated,
                                      sc_dir = sc_dir)

  # comparison type
  comparisonType <- callModule(comparisonType, 'comparison',
                               scseq = scAnal$scseq,
                               is.integrated = scAnal$is.integrated)


  # the selected cluster/gene for cluster comparison ----
  scClusterComparison <- callModule(clusterComparison, 'cluster',
                                    selected_anal = scAnal$selected_anal,
                                    scseq = scAnal$scseq,
                                    markers = scAnal$markers,
                                    annot_path = scAnal$annot_path,
                                    sc_dir = sc_dir,
                                    ref_preds = scLabelTransfer)

  scClusterGene <- callModule(selectedGene, 'gene_clusters',
                              selected_anal = scAnal$selected_anal,
                              scseq = scAnal$scseq,
                              selected_markers = scClusterComparison$selected_markers,
                              cluster_markers = NULL,
                              selected_cluster = scClusterComparison$selected_cluster,
                              annot_path = scAnal$annot_path,
                              comparison_type = comparisonType)

  # the selected clusters/gene for sample comparison ----

  scSampleComparison <- callModule(sampleComparison, 'sample',
                                   selected_anal = scAnal$selected_anal,
                                   scseq = scAnal$scseq,
                                   annot = scClusterComparison$annot,
                                   sc_dir = sc_dir)


  scSampleGene <- callModule(selectedGene, 'gene_samples',
                             selected_anal = scAnal$selected_anal,
                             scseq = scAnal$scseq,
                             selected_markers = scSampleComparison$selected_markers,
                             cluster_marker = scSampleComparison$cluster_markers,
                             selected_cluster = scSampleComparison$selected_cluster,
                             annot_path = scAnal$annot_path,
                             comparison_type = comparisonType)



  # update scseq with annotation changes and jitter ----
  scseq <- reactive({
    scseq <- scAnal$scseq()
    annot <- scClusterComparison$annot()
    jitter <- scAnal$plot_styles$jitter()
    shiny::req(scseq, annot, jitter)

    levels(scseq$seurat_clusters) <- annot
    Seurat::Idents(scseq) <- scseq$seurat_clusters

    if (jitter > 0)
      scseq <- jitter_umap(scseq, amount = jitter)


    return(scseq)
  })

  # show the toggle if dataset is integrated
  observe({
    toggle(id = "comparison_toggle_container",  condition = scAnal$is.integrated())
  })



  # show appropriate inputs based on comparison type
  observe({
    toggle(id = "label_comparison_inputs",  condition = comparisonType() == 'labels')
    toggle(id = "sample_comparison_inputs",  condition = comparisonType() == 'samples')
    toggle(id = "cluster_comparison_inputs", condition = comparisonType() == 'clusters')
  })

  # return values ----


  return(list(
    scseq = scseq,
    plot_styles = scAnal$plot_styles,
    selected_gene_cluster = scClusterGene$selected_gene,
    selected_gene_sample = scSampleGene$selected_gene,
    label_anals = integrationAnnotAnals,
    comparison_type = comparisonType

  ))
}

#' Logic for selected analysis part of scForm
#' @export
#' @keywords internal
selectedAnal <- function(input, output, session, sc_dir, new_anal, new_dataset) {

  selected_anal <- reactive({
    req(input$selected_anal)
    input$selected_anal
  })

  # get's used for saving annotation to disc
  annot_path <- reactive({
    scseq_part_path(sc_dir, selected_anal(), 'annot')
  })

  # load annotation for clusters
  annot <- reactive({
    readRDS(annot_path())
  })

  # load scseq
  scseq <- reactive({

    scseq_path <- scseq_part_path(sc_dir, selected_anal(), 'scseq')
    scseq <- readRDS(scseq_path)

    assay <- get_scseq_assay(scseq)

    if (Seurat::DefaultAssay(scseq) == 'integrated')
      Seurat::DefaultAssay(scseq) <- assay

    Seurat::Idents(scseq) <- scseq$seurat_clusters

    return(scseq)
  })

  is.integrated <- reactive({
    scseq <- scseq()
    req(scseq)

    return('integrated' %in% names(scseq@assays))
  })



  # load markers and name using annot
  markers <- reactive({
    req(annot())

    markers_path <- scseq_part_path(sc_dir, selected_anal(), 'markers')
    markers <- readRDS(markers_path)
    names(markers) <- annot()
    return(markers)
  })



  # available analyses
  anal_options <- reactive({
    # reactive to new anals and new sc datasets
    new_anal()
    new_dataset()

    # make sure integrated rds exists
    int_path <- file.path(sc_dir, 'integrated.rds')
    if (!file.exists(int_path)) saveRDS(NULL, int_path)

    # use saved anals as options
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
  })

  # update if options change
  observe({
    updateSelectizeInput(session, 'selected_anal', choices = anal_options())
  })

  # get styles and integration info
  plot_styles <- callModule(plotStyles, 'styles')
  show_integration <- callModule(showIntegration, 'integration')
  show_label_transfer <- callModule(showLabelTransfer, 'label-transfer')


  # return anal and options to app
  return(list(
    selected_anal = selected_anal,
    scseq = scseq,
    markers = markers,
    annot = annot,
    annot_path = annot_path,
    anal_options = anal_options,
    plot_styles = plot_styles,
    show_integration = show_integration,
    show_label_transfer = show_label_transfer,
    is.integrated = is.integrated
  ))
}

#' Logic for plot styles dropdown in selectedAnal
#' @export
#' @keywords internal
plotStyles <- function(input, output, session) {


  jitter <- reactive(input$point_jitter)
  size <- reactive(input$point_size)

  return(list(jitter = jitter, size = size))
}

#' Logic for show integration button in selectedAnal
#' @export
#' @keywords internal
showIntegration <- function(input, output, session) {

  show_integration <- reactive(input$show_integration %% 2 != 0)


  # show/hide integration form
  observe({
    toggleClass(id = "show_integration", 'btn-primary', condition = show_integration())
  })


  return(show_integration)
}

#' Logic for show label transfer button in selectedAnal
#' @export
#' @keywords internal
showLabelTransfer <- function(input, output, session) {

  show_label_transfer <- reactive(input$show_label_transfer %% 2 != 0)


  # show/hide label transfer form
  observe({
    toggleClass(id = "show_label_transfer", 'btn-primary', condition = show_label_transfer())
  })


  return(show_label_transfer)
}


#' Logic for label transfer between datasets
#' @export
#' @keywords internal
labelTransferForm <- function(input, output, session, sc_dir, anal_options, show_label_transfer, selected_anal, scseq) {
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
    query_name <- selected_anal()

    # load previously saved reference preds
    preds_path <- scseq_part_path(sc_dir, query_name, 'preds')
    if (file.exists(preds_path)) readRDS(preds_path) else list()
  })

  # update annotation transfer choices
  observe({
    preds <- preds()

    anal_options <- anal_options()
    selected_anal <- selected_anal()
    req(preds, anal_options)

    choices <- get_label_transfer_choices(anal_options, selected_anal, preds)
    updateSelectizeInput(session, 'ref_name', choices = choices, server = TRUE, selected = isolate(new_preds()), options = list(render = I('{option: transferLabelOption}')))
  })

  # submit annotation transfer
  observeEvent(input$submit_transfer, {

    query_name <- selected_anal()
    ref_name <- input$ref_name
    preds <- preds()
    req(query_name, ref_name, preds)

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

    # load anals
    query <- scseq()
    ref <- load_saved_scseq(ref_name, sc_dir)
    updateProgress(1/n)

    # transfer labels to query cells
    predictions <- transfer_labels(ref, query, updateProgress = updateProgress, n = n, n_init = 2)
    predictions$orig <- query$seurat_clusters

    # score as sum of percent in cluster with label and mean score
    pred_pcts <- predictions %>%
      group_by(orig, predicted.id) %>%
      summarise(mean.score = mean(prediction.score.max), n = n()) %>%
      arrange(desc(n), desc(mean.score)) %>%
      slice(1)


    preds_path <- scseq_part_path(sc_dir, query_name, 'preds')
    preds[[ref_name]] <- pred_pcts
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
    query_name <- selected_anal()
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
    anal_name <- selected_anal()

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
    anal_name <- selected_anal()

    req(anal_name)

    showModal(transferModal())
  })

  observeEvent(input$confirm_overwrite, {
    removeModal()
    ref_name <- input$ref_name
    ref_preds <- ref_preds()
    anal_name <- selected_anal()

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
integrationForm <- function(input, output, session, sc_dir, anal_options, show_integration) {

  integration_inputs <- c('ctrl_integration', 'integration_name', 'submit_integration', 'test_integration', 'exclude_clusters')


  integration_name <- reactive(input$integration_name)
  integration_options <- reactive(anal_options()$Individual)

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
    exclude_choices <- get_exclude_choices(anal_names, anal_colors, sc_dir)
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
    anal_options <- anal_options()

    error_msg <- validate_integration(test_anals, ctrl_anals, anal_name, anal_options)

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
    return(unique(scseq$project))
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
clusterComparison <- function(input, output, session, selected_anal, scseq, markers, annot_path, sc_dir, ref_preds) {


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
      choices <- get_cluster_choices(clusters, selected_anal(), sc_dir)
    }

    return(choices)
  })


  all_markers <- reactive({
    con_markers <- con_markers()
    markers <- markers()
    names(markers) <- seq(0, along.with = markers)
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
  observeEvent(selected_anal(), {
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
      con_markers[[sel]] <- get_scseq_markers(scseq(), ident.1 = con[1], ident.2 = con[2])
      con_markers(con_markers)
    }
  })


  observe({
    sel <- selected_cluster()
    req(sel)
    selected_markers <- all_markers()[[sel]]

    # remove redundant columns for saving csv
    selected_markers$gene <- NULL
    selected_markers$cluster <- NULL

    selected_markers(selected_markers)
  })


  return(list(
    annot = annot,
    selected_markers = selected_markers,
    selected_cluster = selected_cluster
  ))
}



#' Logic to for sample comparison input
#' @export
#' @keywords internal
sampleComparison <- function(input, output, session, selected_anal, scseq, annot, sc_dir) {
  contrast_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))


  # data.frame of markers for selected sample
  selected_markers <- reactiveVal()
  cluster_markers <- reactiveVal()

  cluster_choices <- reactive({
    req(annot())
    get_cluster_choices(annot(), selected_anal(), sc_dir, sample_comparison = TRUE)
  })

  # reset if switch analysis or annotation updates
  observe({
    selected_anal()
    annot()
    selected_markers(NULL)
    cluster_markers(NULL)
  }, priority = 1)



  # update UI for contrast/cluster choices
  observeEvent(cluster_choices(), {
    updateSelectizeInput(session, 'selected_clusters',
                         choices = cluster_choices(),
                         options = contrast_options, server = TRUE)
  })

  observeEvent(input$run_comparison, {

    selected_clusters <- input$selected_clusters
    req(selected_clusters)

    scseq <- scseq()
    anal_name <- selected_anal()

    res <- run_comparison(scseq, selected_clusters, sc_dir, anal_name)

    # set selected markers
    selected_markers(res$anal$top_table)
    cluster_markers(res$cluster_markers)
  })

  return(list(
    selected_markers = selected_markers,
    cluster_markers = cluster_markers,
    selected_cluster = reactive(input$selected_clusters)
  ))
}



#' Logic for selected gene to show plots for
#' @export
#' @keywords internal
selectedGene <- function(input, output, session, selected_anal, scseq, selected_markers, cluster_markers, selected_cluster, annot_path, comparison_type) {

  selected_gene <- reactiveVal(NULL)

  exclude_ambient <- reactive({
    if (is.null(input$exclude_ambient)) return(FALSE)
    input$exclude_ambient %% 2 != 0
  })

  # toggle for excluding ambient
  observe({
    toggleClass('exclude_ambient', class = 'btn-primary', condition = exclude_ambient())
  })


  filtered_markers <- reactive({
    markers <- selected_markers()
    if (is.null(markers)) return(NULL)

    if (exclude_ambient()) {
      cluster_markers <- cluster_markers()
      scseq <- scseq()
      markers <- ambient.omit(scseq,  markers = markers, cluster_markers = cluster_markers)
    }


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
    if (comparison_type == 'samples' & !'t' %in% colnames(markers)) return(NULL)
    if (comparison_type == 'clusters' & !'pct.1' %in% colnames(markers)) return(NULL)

    get_gene_choices(scseq,
                     markers = markers,
                     selected_cluster = selected_cluster,
                     comparison_type = comparison_type)

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
                     anal = selected_anal(),
                     comparison_type = comparison_type())

    },
    content = function(con) {
      write.csv(filtered_markers(), con)
    }
  )





  # reset selected gene if analysis changes
  observeEvent(selected_anal(), {
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
                         options = list(render = I('{option: geneOption, item: geneItem}')),
                         server = TRUE)
  })

  return(list(
    selected_gene = selected_gene
  ))

}


#' Logic for cluster plots
#' @export
#' @keywords internal
scClusterPlot <- function(input, output, session, scseq, plot_styles, selected_group, cached_plot) {

  plot <- reactive({
    cached_plot <- cached_plot()
    if (!is.null(cached_plot)) return(cached_plot)

    plot_umap_cluster(scseq(), pt.size = plot_styles$size())

  })

  output$cluster_plot <- renderPlot({
    plot()
  })

  return(list(
    plot = plot
  ))
}

#' Logic for marker gene plots
#' @export
#' @keywords internal
scMarkerPlot <- function(input, output, session, scseq, selected_gene, plot_styles, selected_group, cached_plot) {


  plot <- reactive({
    # cached plot if showing test samples
    cached <- cached_plot()
    if (!is.null(cached)) return(cached)

    req(selected_gene())
    plot_umap_gene(scseq(), selected_gene(), pt.size = plot_styles$size())
  })

  output$marker_plot <- renderPlot({
    pl <- plot()
    req(pl)

    if (selected_group != 'all')
      pl <- format_sample_gene_plot(pl, selected_group, selected_gene(), scseq())

    return(pl)
  })


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

