# page and input form logic ----

#' Logic for Single Cell Exploration page
#' @export
#' @keywords internal
scPage <- function(input, output, session, data_dir) {

  # the analysis and options
  scForm <- callModule(scForm, 'form',
                       data_dir = data_dir)

  callModule(scClusterPlot, 'cluster_plot',
             scseq = scForm$scseq,
             plot_styles = scForm$plot_styles)

  callModule(scMarkerPlot, 'marker_plot',
             scseq = scForm$scseq,
             selected_gene = scForm$selected_gene,
             selected_groups = scForm$selected_groups,
             plot_styles = scForm$plot_styles)


  callModule(scBioGpsPlot, 'biogps_plot',
             selected_gene = scForm$selected_gene)

  return(NULL)
}

scForm <- function(input, output, session, data_dir) {

  # updates if new integrated dataset
  new_anal <- reactive({
    scIntegration()
  })


  # the analysis and options
  scAnal <- callModule(selectedAnal, 'anal',
                       data_dir = data_dir,
                       new_anal = new_anal)

  # integration stuff
  scIntegration <- callModule(integrationForm, 'integration',
                              data_dir = data_dir,
                              anal_options = scAnal$anal_options,
                              show_integration = scAnal$show_integration)

  # the cluster
  scCluster <- callModule(selectedCluster, 'cluster',
                          selected_anal = scAnal$selected_anal,
                          scseq = scAnal$scseq,
                          markers = scAnal$markers,
                          annot_path = scAnal$annot_path)


  # the gene selection
  scGene <- callModule(selectedGene, 'gene',
                       selected_anal = scAnal$selected_anal,
                       selected_cluster = scCluster$selected_cluster,
                       scseq = scAnal$scseq,
                       selected_markers = scCluster$selected_markers)


  # the groups selection
  scGroups <- callModule(selectedGroups, 'groups', scseq = scAnal$scseq)


  # update scseq with annotation changes and jitter
  scseq <- reactive({
    scseq <- scAnal$scseq()
    annot <- scCluster$annot()
    jitter <- scAnal$plot_styles$jitter()
    shiny::req(scseq, annot, jitter)

    levels(scseq$seurat_clusters) <- annot
    Seurat::Idents(scseq) <- scseq$seurat_clusters

    if (jitter > 0)
      scseq <- jitter_umap(scseq, amount = jitter)


    return(scseq)
  })


  return(list(
    scseq = scseq,
    plot_styles = scAnal$plot_styles,
    selected_gene = scGene$selected_gene,
    selected_groups = scGroups$selected_groups
  ))
}

# selected analysis input logic ----
selectedAnal <- function(input, output, session, data_dir, new_anal) {

  selected_anal <- reactive({
    req(input$selected_anal)
    input$selected_anal
  })

  # get's used for saving annotation to disc
  annot_path <- reactive({
    scseq_part_path(data_dir, selected_anal(), 'annot')
  })

  # load annotation for clusters
  annot <- reactive({
    readRDS(annot_path())
  })

  # load scseq and name using annot
  scseq <- reactive({

    scseq_path <- scseq_part_path(data_dir, selected_anal(), 'scseq')
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

    markers_path <- scseq_part_path(data_dir, selected_anal(), 'markers')
    markers <- readRDS(markers_path)
    names(markers) <- annot()
    return(markers)
  })



  # available analyses
  anal_options <- reactive({
    # reactive to new anals
    new_anal()

    # make sure integrated rds exists
    int_path <- file.path(data_dir, 'integrated.rds')
    if (!file.exists(int_path)) saveRDS(NULL, int_path)

    # use saved anals as options
    integrated <- readRDS(file.path(data_dir, 'integrated.rds'))
    individual <- setdiff(list.files(data_dir), c(integrated, 'integrated.rds'))

    # must be a list if length one for option groups to work
    if (length(integrated) == 1) integrated <- list(integrated)
    if (length(individual) == 1) individual <- list(individual)

    list(Individual = individual, Integrated = integrated)
  })

  # update if options change
  observe({
    updateSelectizeInput(session, 'selected_anal', choices = anal_options())
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

plotStyles <- function(input, output, session) {


  jitter <- reactive(input$point_jitter)
  size <- reactive(input$point_size)

  return(list(jitter = jitter, size = size))
}

showIntegration <- function(input, output, session) {

  show_integration <- reactive(input$show_integration %% 2 != 0)


  # show/hide integration form
  observe({
    toggleClass(id = "show_integration", 'active', condition = show_integration())
  })


  return(show_integration)
}

# single cell dataset integration logic -----

#' Single Cell Integration form
#'
#' @return \code{reactiveVal} that is either NULL or contains the name of a new integrated analysis
integrationForm <- function(input, output, session, data_dir, anal_options, show_integration) {


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
      removeClass(id = 'validate', class = 'has-error')
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
      integrate_saved_scseqs(data_dir, ctrl_anals, ctrl_anals, anal_name, updateProgress = updateProgress)


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

# selected cluster/contrast/rename logic ----
selectedCluster <- function(input, output, session, selected_anal, scseq, markers, annot_path) {


  contrast_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))
  selected_cluster <- reactiveVal(NULL)
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
  contrast_choices <- reactive({
    clusters <- annot()
    req(clusters)

    if (show_contrasts()) {
      # group choices are as compared to other clusters
      test <- isolate(test_cluster())
      ctrls <- clusters[clusters != test]

      colours <- get_palette(clusters)
      names(colours) <- clusters

      contrast_choices <- data.frame(test = stringr::str_trunc(test, 17),
                                     ctrl = stringr::str_trunc(c('all', ctrls), 17),
                                     value = c(test, paste0(test, ' vs ', ctrls)),
                                     testColor = colours[test],
                                     ctrlColor = c('white', colours[ctrls]), row.names = NULL)


    } else {
      # show the cell numbers/percentages
      ncells <- tabulate(scseq()$seurat_clusters)
      pcells <- round(ncells / sum(ncells) * 100)
      pspace <- strrep('&nbsp;&nbsp;', 2 - nchar(pcells))

      # cluster choices are the clusters themselves
      testColor <- get_palette(clusters)
      contrast_choices <- data.frame(name = stringr::str_trunc(clusters, 27),
                                     value = clusters,
                                     label = clusters,
                                     testColor,
                                     ncells, pcells, pspace, row.names = NULL)
    }

    return(contrast_choices)
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
    toggleClass(id = "show_contrasts", 'active', condition = show_contrasts())
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
  observeEvent(contrast_choices(), {
    updateSelectizeInput(session, 'selected_cluster',
                         choices = contrast_choices(),
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

# selected marker gene logic ----
selectedGene <- function(input, output, session, selected_anal, selected_cluster, scseq, selected_markers) {

  selected_gene <- reactiveVal(NULL)

  # update marker genes based on cluster selection
  gene_choices <- reactive({
    selected_markers <- selected_markers()
    req(selected_markers)

    # allow selecting non-marker genes (at bottom of list)
    choices <- row.names(selected_markers)
    choices <- c(choices, setdiff(row.names(scseq()), choices))

    return(choices)
  })

  output$genecards <- renderUI({
    gene_link <- paste0('https://www.genecards.org/cgi-bin/carddisp.pl?gene=', input$selected_gene)
    ns <- session$ns
    withTags({
      a(class = 'btn btn-default',
        href = gene_link, target = '_blank',
        icon('external-link-alt', 'fa-fw')
      )

    })
  })

  # reset selected gene if analysis changes
  observeEvent(selected_anal(), {
    selected_gene(NULL)
  })

  # set return value based on input
  observeEvent(input$selected_gene, {
    selected_gene(input$selected_gene)
  })


  observe({
    updateSelectizeInput(session, 'selected_gene', choices = gene_choices(), selected = NULL, server = TRUE)
  })

  return(list(
    selected_gene = selected_gene
  ))

}

# selected groups (ctrl/test/all) logic -----
selectedGroups <- function(input, output, session, scseq) {

  # groups to show (e.g. ctrl and test)
  available_groups <- shiny::reactive({
    unique(as.character(scseq()$orig.ident))
  })

  show_groups <- reactive({
    length(available_groups()) > 1
  })

  selected_groups <- reactive({
    groups <- available_groups()
    # always show when just a single group
    if (length(groups) == 1 || input$selected_group == 'all') return(groups)

    return(input$selected_group)
  })

  # show UI component if more than one available groups
  observe({
    toggle(id = "selected_group_container", condition = show_groups())
  })


  return(list(
    selected_groups = selected_groups
  ))
}

# plot logic ----
scClusterPlot <- function(input, output, session, scseq, plot_styles) {

  output$cluster_plot <- renderPlot({
    plot_umap_cluster(scseq(), pt.size = plot_styles$size())
  })
}

scMarkerPlot <- function(input, output, session, scseq, selected_gene, selected_groups, plot_styles) {

  output$marker_plot <- renderPlot({
    req(selected_gene())
    req(selected_groups())
    plot_umap_gene(scseq(), selected_gene(), selected_idents = selected_groups(), pt.size = plot_styles$size())
  })
}

scBioGpsPlot <- function(input, output, session, selected_gene) {
  # plot BioGPS data
  output$biogps_plot <- renderPlot({
    plot_biogps(selected_gene())
  })
}
