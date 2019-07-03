server <- function(input, output, session) {
  shinyjs::useShinyjs(html = TRUE)

  # get arguments from calling function
  # defaults for server
  data_dir <- shiny::getShinyOption('data_dir', '/srv/shiny-server/drugseqr/scseq/sjia')
  pt.size  <- shiny::getShinyOption('pt.size', 2.5)


  # toggle to show dataset integration
  shinyjs::onclick("show_integration", {
  shinyjs::toggle(id = "integration", anim = TRUE)
  shinyjs::toggleClass(id = "show_integration", 'active')
  })

  # reactive values (can update and persist within session) -----
  new_anals_rv <- shiny::reactiveVal(NULL)
  annot_rv <- shiny::reactiveVal(NULL)
  annot_path_rv <- shiny::reactiveVal(NULL)
  con_markers_rv <- shiny::reactiveVal(list())
  selected_cluster_rv <- shiny::reactiveVal(NULL)
  selected_cluster_idx_rv <- shiny::reactiveVal(NULL)


  # reactive expressions (auto update) ----------

  # available analyses
  anal_options_r <- shiny::reactive({
    # reactive to new anals
    new_anals_rv()

    # use saved anals as options
    anal_files <- rev(list.files(data_dir))
    anal_files <- anal_files[!grepl('_annot.rds$', anal_files)]
    anal_options <- gsub('.rds$', '', anal_files)

    return(anal_options)
  })

  # analyses available for integration
  integration_options_r <- shiny::reactive({
    # TODO: integrated datasets not an option
    anal_options_r()
  })

  # groups to show (e.g. ctrl and test)
  available_groups_r <- shiny::reactive({
    scseq <- scseq_r()
    groups <- unique(as.character(scseq$orig.ident))
    return(groups)
  })

  # current analysis
  anal_r <- shiny::eventReactive(input$selected_anal, {

    selected_file <- paste0(input$selected_anal, '.rds')
    anal <- readRDS(file.path(data_dir, selected_file))

    if (Seurat::DefaultAssay(anal$scseq) == 'integrated')
      Seurat::DefaultAssay(anal$scseq) <- 'SCT'

    return(anal)
  }, ignoreInit = TRUE)

  # the markers
  markers_r <- shiny::reactive({
    anal <- anal_r()
    annot <- annot_rv()
    markers <- anal$markers
    names(markers) <- annot

    return(markers)
  })

  is_rename_r <- shiny::reactive({
    is_rename <- (input$rename_cluster + input$show_rename) %% 2 == 0
    return(is_rename)
  })


  # the seurat object
  scseq_r <- shiny::reactive({
    anal <- anal_r()
    markers <- markers_r()
    scseq <- anal$scseq

    levels(scseq$seurat_clusters) <- names(markers)
    Seurat::Idents(scseq) <- scseq$seurat_clusters

    return(scseq)
  })

  # the cluster names
  clusters_r <- shiny::reactive({
    markers <- markers_r()
    return(names(markers))
  })

  # currently selected test cluster (need to)
  test_cluster_r <- shiny::reactive({
    test_cluster <- input$selected_cluster
    test_cluster <- gsub(' vs .+?$', '', test_cluster)
    if (test_cluster == '') return(NULL)
    return(test_cluster)
  })

  # used to determine cell groups to show (e.g. ctrl and test)
  selected_groups_r <- shiny::reactive({

    groups <- available_groups_r()
    # always show when just a single group
    if (length(groups) == 1 || input$groups == 'all') return(groups)

    return(input$groups)
  })



  # observations (do stuff if something changes) -------
  # analysis options
  shiny::observe({
    shiny::updateSelectizeInput(session, 'selected_anal', choices = anal_options_r())
  })

  # integration analyses can be either control or test (not both)
  # update choices of opposite so that doesn't close current
  shiny::observe({

    ctrl <- input$ctrl_integration
    test <- shiny::isolate(input$test_integration)
    anal_options <- anal_options_r()

    shiny::updateSelectizeInput(session, 'test_integration', choices = anal_options[!anal_options %in% ctrl], selected = test)
  })
  shiny::observe({
    ctrl <- shiny::isolate(input$ctrl_integration)
    test <- input$test_integration
    anal_options <- anal_options_r()

    shiny::updateSelectizeInput(session, 'ctrl_integration', choices = anal_options[!anal_options %in% test], selected = ctrl)
  })


  # setup the initial cluster annotations/file path
  shiny::observeEvent(anal_r(), {
    anal <- anal_r()
    annot <- anal$annot

    if (is.null(annot))
      annot <- names(anal$markers)

    # make sure annot on disc so that can update quickly
    annot_path <- file.path(data_dir, paste0(input$selected_anal, '_annot.rds'))
    if (!file.exists(annot_path)) saveRDS(annot, annot_path)

    annot <- readRDS(annot_path)
    annot_rv(annot)
    annot_path_rv(annot_path)
    selected_cluster_rv(NULL)
  }, priority = 1)


  # update annot if rename a cluster
  shiny::observeEvent(input$rename_cluster, {

    if (is_rename_r()) {

      if (input$new_cluster_name != '') {

        # update reactive annotation
        annot <- annot_rv()
        sel.clust <- input$selected_cluster
        sel.idx   <- which(annot == sel.clust)
        annot[sel.idx] <- input$new_cluster_name
        annot <- make.unique(annot)
        annot_rv(annot)
        selected_cluster_idx_rv(sel.idx)

        # save on disc
        saveRDS(annot, annot_path_rv())
      }

    }
  })

  # update selected cluster if rename a cluster
  shiny::observeEvent(input$rename_cluster, {
    if (is_rename_r() & input$new_cluster_name != '') {
      annot <- annot_rv()
      sel.idx <- selected_cluster_idx_rv()
      selected_cluster_rv(annot[sel.idx])
    }
  })


  # update group choices/selected if they change or show_contrasts changes
  shiny::observe({

    scseq <- scseq_r()
    clusters <- clusters_r()

    if (input$show_contrasts %% 2 != 0) {
      # group choices are as compared to other clusters
      test <- shiny::isolate(test_cluster_r())
      ctrls <- clusters[clusters != test]

      colours <- get_palette(clusters)
      names(colours) <- clusters

      contrast_choices <- data.frame(test = stringr::str_trunc(test, 17),
                                     ctrl = stringr::str_trunc(c('all', ctrls), 17),
                                     value = c(test, paste0(test, ' vs ', ctrls)),
                                     testColor = colours[test],
                                     ctrlColor = c('white', colours[ctrls]))


    } else {
      # show the cell numbers/percentages
      ncells <- tabulate(scseq$seurat_clusters)
      pcells <- round(ncells/sum(ncells) * 100)
      pspace <- strrep('&nbsp;&nbsp;', 2-nchar(pcells))

      # cluster choices are the clusters themselves
      testColor <- get_palette(clusters)
      trunc <- stringr::str_trunc(clusters, 27)
      contrast_choices  <- data.frame(name = trunc,
                                      value = clusters,
                                      label = clusters,
                                      testColor,
                                      ncells,
                                      pcells, pspace)
    }

    contrast_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))

    shiny::updateSelectizeInput(session, 'selected_cluster',
                                choices = contrast_choices, selected = selected_cluster_rv(),
                                options = contrast_options, server = TRUE)
  })


  # update group buttons if show contrasts is toggled
  shiny::observeEvent(input$show_contrasts, {
    # update icon on toggle
    icon <- ifelse(input$show_contrasts %% 2 != 0, 'chevron-down', 'chevron-right')

    shiny::updateActionButton(session, 'show_contrasts', icon = shiny::icon(icon, 'fa-fw'))
    shinyjs::toggleState('show_rename')
  })

  # update selected cluster if show contrasts is toggled
  shiny::observeEvent(input$show_contrasts, {
    if (input$show_contrasts %% 2 != 0) {
      selected_cluster_rv(NULL)
    } else {
      selected_cluster_rv(test_cluster_r())
    }

  })

  # update marker genes based on cluster selection
  shiny::observeEvent(input$selected_cluster, {
    scseq <- scseq_r()
    con_markers <- con_markers_rv()
    markers <- c(markers_r(), con_markers)

    sel <- input$selected_cluster
    if (sel == '') return(NULL)

    # get cluster if don't have (for comparing specific cluster)
    if (!sel %in% names(markers)) {
      con <- strsplit(sel, ' vs ')[[1]]
      con_markers[[sel]] <- markers[[sel]] <- get_scseq_markers(scseq, ident.1 = con[1], ident.2 = con[2])

      # update con markers reactive value
      con_markers_rv(con_markers)
    }

    # allow selecting non-marker genes (at bottom of list)
    cluster_markers <- markers[[sel]]
    choices <- row.names(cluster_markers)
    choices <- c(choices, setdiff(row.names(scseq), choices))

    shiny::updateSelectizeInput(session, 'gene', choices = choices, selected = NULL, server = TRUE)
  })


  # dynamic UI elements (change based on state of app) ----

  # ui to show e.g. ctrl or test cells
  output$groups_toggle <- shiny::renderUI({
    # if more than one group allow showing cells based on groups
    groups <- available_groups_r()
    groups_toggle <- NULL
    if (length(groups) > 1) {
      groups_toggle <- shiny::tags$div(
        shiny::radioButtons("groups", "Show cells for group:",
                            choices = c('all', groups), inline = FALSE),
        shiny::br()
      )
    }
    return(groups_toggle)
  })

  # ui for renaming a cluster
  shiny::observeEvent(input$show_rename, {
    if (!is_rename_r())
      shiny::updateTextInput(session, 'new_cluster_name', value = '', placeholder = paste('Type new name for', input$selected_cluster, '...'))
  })


  # plots ------
  # show tSNE plot coloured by expression values
  output$marker_plot <- shiny::renderPlot({

    scseq <- scseq_r()

    if (input$gene == '' || !input$gene %in% row.names(scseq)) return(NULL)
    plot_umap_gene(scseq, input$gene, selected_idents = selected_groups_r(), pt.size = pt.size)
  })


  # show plot of predicted cell clusters
  output$cluster_plot <- shiny::renderPlot({
    scseq <- scseq_r()


    legend_title <-'Cluster'
    plot_umap_cluster(scseq, pt.size = pt.size)
  })


  # plot BioGPS data
  output$biogps <- shiny::renderPlot({
    plot_biogps(input$gene)
  })


  # link to GeneCards page for gene ----
  gene_link_r <- shiny::reactive({
    return(paste0('https://www.genecards.org/cgi-bin/carddisp.pl?gene=', input$gene))
  })

  shiny::observeEvent(input$genecards, {
    utils::browseURL(gene_link_r())
  })

  shiny::observeEvent(input$done, {
    shiny::stopApp()
  })

}
