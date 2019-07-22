#' Logic for Single Cell Exploration page
#' @export
#' @keywords internal
scPage <- function(input, output, session, sc_dir) {

  # the analysis and options
  scForm <- callModule(scForm, 'form',
                       sc_dir = sc_dir)


  callModule(scClusterPlot, 'cluster_plot',
             scseq = scForm$scseq,
             plot_styles = scForm$plot_styles)

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

  show_samples <- reactive({
    scForm$comparison_type() == 'samples'
  })

  observe({
    toggle(id = "sample_comparison_row",  condition = show_samples())
    toggle(id = "cluster_comparison_row", condition = !show_samples())
  })



  return(NULL)
}

#' Logic for form on Single Cell Exploration page
#' @export
#' @keywords internal
scForm <- function(input, output, session, sc_dir) {

  # updates if new integrated dataset
  new_anal <- reactive({
    scIntegration()
  })


  # the analysis and options
  scAnal <- callModule(selectedAnal, 'anal',
                       sc_dir = sc_dir,
                       new_anal = new_anal)

  # integration stuff
  scIntegration <- callModule(integrationForm, 'integration',
                              sc_dir = sc_dir,
                              anal_options = scAnal$anal_options,
                              show_integration = scAnal$show_integration)

  # comparison type
  comparisonType <- callModule(comparisonType, 'comparison',
                               scseq = scAnal$scseq)

  # the selected cluster/gene for cluster comparison ----
  scClusterComparison <- callModule(clusterComparison, 'cluster',
                                    selected_anal = scAnal$selected_anal,
                                    scseq = scAnal$scseq,
                                    markers = scAnal$markers,
                                    annot_path = scAnal$annot_path)

  scClusterGene <- callModule(selectedGene, 'gene_clusters',
                              selected_anal = scAnal$selected_anal,
                              selected_cluster = scClusterComparison$selected_cluster,
                              scseq = scAnal$scseq,
                              selected_markers = scClusterComparison$selected_markers)

  # the selected clusters/gene for sample comparison ----



  scSampleComparison <- callModule(sampleComparison, 'sample',
                                   selected_anal = scAnal$selected_anal,
                                   scseq = scAnal$scseq,
                                   annot = scClusterComparison$annot)


  scSampleGene <- callModule(selectedGene, 'gene_samples',
                             selected_anal = scAnal$selected_anal,
                             selected_cluster = scSampleComparison$selected_cluster,
                             scseq = scAnal$scseq,
                             selected_markers = scSampleComparison$selected_markers)



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

  show_samples <- reactive({
    comparisonType() == 'samples'
  })


  observe({
    toggle(id = "sample_comparison_inputs",  condition = show_samples())
    toggle(id = "cluster_comparison_inputs", condition = !show_samples())
  })

  # return values ----


  return(list(
    scseq = scseq,
    plot_styles = scAnal$plot_styles,
    selected_gene_cluster = scClusterGene$selected_gene,
    selected_gene_sample = scSampleGene$selected_gene,
    comparison_type = comparisonType

  ))
}


#' Logic for selected analysis part of scForm
#' @export
#' @keywords internal
selectedAnal <- function(input, output, session, sc_dir, new_anal) {

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

  # load scseq and name using annot
  scseq <- reactive({

    scseq_path <- scseq_part_path(sc_dir, selected_anal(), 'scseq')
    scseq <- readRDS(scseq_path)

    if (Seurat::DefaultAssay(scseq) == 'integrated')
      Seurat::DefaultAssay(scseq) <- 'SCT'

    levels(scseq$seurat_clusters) <- annot()
    Seurat::Idents(scseq) <- scseq$seurat_clusters

    return(scseq)
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
    # reactive to new anals
    new_anal()

    # make sure integrated rds exists
    int_path <- file.path(sc_dir, 'integrated.rds')
    if (!file.exists(int_path)) saveRDS(NULL, int_path)

    # use saved anals as options
    integrated <- readRDS(file.path(sc_dir, 'integrated.rds'))
    individual <- setdiff(list.files(sc_dir), c(integrated, 'integrated.rds'))

    # must be a list if length one for option groups to work
    if (length(integrated) == 1) integrated <- list(integrated)
    if (length(individual) == 1) individual <- list(individual)

    list(Individual = individual, Integrated = integrated)
  })

  # update if options change
  observe({
    updateSelectizeInput(session, 'selected_anal', selected = 'sjia_lung', choices = anal_options())
  })

  # get styles and integration info
  plot_styles <- callModule(plotStyles, 'styles')
  show_integration <- callModule(showIntegration, 'integration')


  # return anal and options to app
  return(list(
    selected_anal = selected_anal,
    scseq = scseq,
    markers = markers,
    annot = annot,
    annot_path = annot_path,
    anal_options = anal_options,
    plot_styles = plot_styles,
    show_integration = show_integration
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


#' Logic for integration form toggled by showIntegration
#' @export
#' @keywords internal
integrationForm <- function(input, output, session, sc_dir, anal_options, show_integration) {


  integration_name <- reactive(input$integration_name)
  integration_options <- reactive(anal_options()$Individual)

  ctrl <- reactiveVal()
  test <- reactiveVal()
  new_anal <- reactiveVal()

  # show/hide integration form
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

  # run integration
  observeEvent(input$submit_integration, {

    test_anals <- test()
    ctrl_anals <- ctrl()
    anal_name <- input$integration_name
    anal_options <- anal_options()

    error_msg <- validate_integration(test_anals, ctrl_anals, anal_name, anal_options)

    if (is.null(error_msg)) {
      # clear error and disable button
      removeClass('validate', 'has-error')
      disable('submit_integration')

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
                             anal_name = anal_name,
                             updateProgress = updateProgress)


      # re-enable, clear inputs, and trigger update of available anals
      ctrl(NULL)
      test(NULL)
      new_anal(anal_name)
      updateTextInput(session, 'integration_name', value = '')
      enable('submit_integration')

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
comparisonType <- function(input, output, session, scseq) {

  # groups to show (e.g. ctrl and test)
  available_groups <- shiny::reactive({
    unique(as.character(scseq()$orig.ident))
  })

  show_groups <- reactive({
    length(available_groups()) > 1
  })

  # show UI component if more than one available groups
  observe({
    toggle(id = "comparison_type_container", condition = show_groups())
  })

  # always show clusters if not integrated
  observe({
    if( !show_groups())
      updateRadioGroupButtons(session, 'comparison_type', selected = 'clusters')
  })

  return(reactive(input$comparison_type))
}

#' @export
#' @keywords internal#' Logic for cluster comparison input

clusterComparison <- function(input, output, session, selected_anal, scseq, markers, annot_path) {


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
      choices <- get_cluster_choices(clusters, scseq())
    }

    return(choices)
  })


  all_markers <- reactive({
    con_markers <- con_markers()
    markers <- markers()
    names(markers) <- annot()
    return(c(markers, con_markers))
  })

  # update scseq with annotation
  annot_scseq <- reactive({
    scseq <- scseq()
    annot <- annot()
    req(annot, scseq)

    levels(scseq$seurat_clusters) <- annot
    Seurat::Idents(scseq) <- scseq$seurat_clusters
    return(scseq)
  })

  observe({ annot(readRDS(annot_path())) })

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
    mod_annot <- annot()
    sel_clust <- selected_cluster()
    sel_idx <- which(mod_annot == sel_clust)
    mod_annot[sel_idx] <- input$new_cluster_name
    mod_annot <- make.unique(mod_annot)

    # save on disc
    saveRDS(mod_annot, annot_path())

    # update annot and set selected cluster to new name
    annot(mod_annot)
    selected_cluster(input$new_cluster_name)

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
    if (!show_rename())
      updateTextInput(session, 'new_cluster_name', value = '', placeholder = paste('Type new name for', input$selected_cluster, '...'))
  })


  # get cluster if don't have (for comparing specific cluster)
  observeEvent(input$selected_cluster, {
    sel <- input$selected_cluster
    req(sel)

    if (!sel %in% names(all_markers())) {
      con <- strsplit(sel, ' vs ')[[1]]
      con_markers <- con_markers()
      con_markers[[sel]] <- get_scseq_markers(annot_scseq(), ident.1 = con[1], ident.2 = con[2])
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
    selected_markers = selected_markers
  ))
}

#' @export
#' @keywords internal#' Logic to for sample comparison input

sampleComparison <- function(input, output, session, selected_anal, scseq, annot) {
  contrast_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))


  # TODO: make it so that annotation changes don't require re-computing markers
  # TODO store markers that have been previously computed

  # data.frame of markers for selected sample
  selected_markers <- reactiveVal()

  cluster_choices <- reactive({
    req(annot())
    get_cluster_choices(annot(), scseq())
  })

  # reset if switch analysis or annotation updates
  observe({
    selected_anal()
    annot()
    selected_markers(NULL)
  }, priority = 1)

  # update scseq with annotation
  annot_scseq <- reactive({
    scseq <- scseq()
    annot <- annot()
    req(annot, scseq)

    levels(scseq$seurat_clusters) <- annot
    Seurat::Idents(scseq) <- scseq$seurat_clusters
    return(scseq)
  })


  # update UI for contrast/cluster choices
  observeEvent(cluster_choices(), {
    updateSelectizeInput(session, 'selected_clusters',
                         choices = cluster_choices(),
                         options = contrast_options, server = TRUE)
  })

  observeEvent(input$run_comparison, {
    req(input$selected_clusters)

    # set idents to ctrl and test
    scseq <- annot_scseq()
    Seurat::Idents(scseq) <- scseq$orig.ident

    # exclude non-selected clusters
    scseq <-  scseq[, scseq$seurat_clusters %in% input$selected_clusters]
    markers <- get_scseq_markers(scseq, ident.1 = 'test', ident.2 = 'ctrl', min.diff.pct = -Inf, test.use = 'limma')

    # set selected markers
    selected_markers(markers)
  })

  return(list(
    selected_markers = selected_markers
  ))
}


#' @export
#' @keywords internal#' Logic for selected gene to show plots for

selectedGene <- function(input, output, session, selected_anal, selected_cluster, scseq, selected_markers) {

  selected_gene <- reactiveVal(NULL)

  exclude_ambient <- reactive({
    if (is.null(input$exclude_ambient)) return(FALSE)
    input$exclude_ambient %% 2 != 0
  })

  # toggle for excluding ambient
  observe({
    toggleClass('exclude_ambient', class = 'btn-primary', condition = exclude_ambient())
  })



  # update marker genes based on cluster selection
  gene_choices <- reactive({
    selected_markers <- selected_markers()
    if (is.null(selected_markers)) return(NULL)

    choices <- row.names(selected_markers)

    if (exclude_ambient()) {
      scseq <- scseq()

      fts <- scseq[['SCT']]@meta.features
      test.ambient <- row.names(fts)[fts$test_ambient]
      ctrl.ambient <- row.names(fts)[fts$ctrl_ambient]

      # exclude test/ctrl ambient if positive/negative effect size
      # opposite would decrease extent of gene expression difference but not direction
      pos.test <- selected_markers[test.ambient, 't'] > 0
      neg.ctrl <- selected_markers[ctrl.ambient, 't'] < 0

      test.exclude <- test.ambient[pos.test]
      ctrl.exclude <- ctrl.ambient[neg.ctrl]

      choices <- setdiff(choices, c(test.exclude, ctrl.exclude))
    }

    # allow selecting non-marker genes (at bottom of list)
    choices <- c(choices, setdiff(row.names(scseq()), choices))

    return(choices)
  })


  observeEvent(input$genecards, {
    gene_link <- paste0('https://www.genecards.org/cgi-bin/carddisp.pl?gene=', input$selected_gene)
    runjs(paste0("window.open('", gene_link, "')"))
  })


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
    updateSelectizeInput(session, 'selected_gene', choices = gene_choices(), selected = NULL, server = TRUE)
  })

  return(list(
    selected_gene = selected_gene
  ))

}


#' @export
#' @keywords internal#' Logic for cluster plots

scClusterPlot <- function(input, output, session, scseq, plot_styles) {

  output$cluster_plot <- renderPlot({
    plot_umap_cluster(scseq(), pt.size = plot_styles$size())
  })
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


#' @export
#' @keywords internal#' Logic for BioGPS plot

scBioGpsPlot <- function(input, output, session, selected_gene) {
  # plot BioGPS data
  output$biogps_plot <- renderPlot({
    plot_biogps(selected_gene())
  })
}

