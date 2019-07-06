 # toggle to show dataset integration
  shinyjs::onclick("show_integration", {
    shinyjs::toggle(id = "integration-form", anim = TRUE)
    shinyjs::toggleClass(id = "show_integration", 'active')
  })

  # reactive values (can update and persist within session) -----
  new_anal_rv <- shiny::reactiveVal(NULL)
  annot_rv <- shiny::reactiveVal(NULL)
  annot_path_rv <- shiny::reactiveVal(NULL)
  con_markers_rv <- shiny::reactiveVal(list())
  selected_cluster_rv <- shiny::reactiveVal(NULL)
  selected_cluster_idx_rv <- shiny::reactiveVal(NULL)


  # reactive expressions (auto update) ----------

  # available analyses
  anal_options_r <- shiny::reactive({
    # reactive to new anals
    new_anal_rv()

    # make sure integrated rds exists
    int_path <- file.path(data_dir, 'integrated.rds')
    if (!file.exists(int_path)) saveRDS(NULL, int_path)

    # use saved anals as options
    int_options <- readRDS(file.path(data_dir, 'integrated.rds'))
    ind_options <- setdiff(list.files(data_dir), c(int_options, 'integrated.rds'))

    # must be a list if length one for option groups to work
    if (length(int_options) == 1) int_options <- list(int_options)
    if (length(ind_options) == 1) ind_options <- list(ind_options)

    anal_options <- list(Individual = ind_options, Integrated = int_options)

    return(anal_options)
  })

  # analyses available for integration
  integration_options_r <- shiny::reactive({
    # TODO: integrated datasets not an option
    anal_options_r()$Individual
  })

  # groups to show (e.g. ctrl and test)
  available_groups_r <- shiny::reactive({
    scseq <- scseq_r()
    groups <- unique(as.character(scseq$orig.ident))
    return(groups)
  })

  # current analysis
  anal_r <- shiny::eventReactive(input$selected_anal, {

    # load analysis parts
    anal_name <- input$selected_anal
    scseq_path <- scseq_part_path(data_dir, anal_name, 'scseq')
    annot_path <- scseq_part_path(data_dir, anal_name, 'annot')
    markers_path <- scseq_part_path(data_dir, anal_name, 'markers')

    anal <- list(scseq = readRDS(scseq_path),
                 markers = readRDS(markers_path),
                 annot = readRDS(annot_path))

    if (Seurat::DefaultAssay(anal$scseq) == 'integrated')
      Seurat::DefaultAssay(anal$scseq) <- 'SCT'

    # reset session-persistent things
    annot_rv(anal$annot)
    annot_path_rv(annot_path)
    selected_cluster_rv(NULL)

    # whether analysis is integrated or not

    return(anal)
  }, ignoreInit = TRUE)


  # used to determine if showing integrated specific options
  output$is_integrated <- shiny::reactive({
    scseq <- anal_r()$scseq
    length(levels(scseq$orig.ident)) == 2
  })
  shiny::outputOptions(output, "is_integrated", suspendWhenHidden = FALSE)


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

    jitter <- input$point_jitter
    if (jitter > 0)
      scseq <- jitter_umap(scseq, factor = jitter)

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
    if (length(groups) == 1 || input$selected_group == 'all') return(groups)

    return(input$selected_group)
  })



  # observations (do stuff if something changes) -------


  # analysis options
  shiny::observe({
    shiny::updateSelectizeInput(session, 'selected_anal', choices = anal_options_r(), selected = 'lung_sjia')
  })

  # integration analyses can be either control or test (not both)
  # update choices of opposite so that doesn't close current
  shiny::observe({

    ctrl <- input$ctrl_integration
    test <- shiny::isolate(input$test_integration)
    anal_options <- integration_options_r()

    shiny::updateSelectizeInput(session, 'test_integration', choices = anal_options[!anal_options %in% ctrl], selected = test)
  })
  shiny::observe({
    ctrl <- shiny::isolate(input$ctrl_integration)
    test <- input$test_integration
    anal_options <- integration_options_r()

    shiny::updateSelectizeInput(session, 'ctrl_integration', choices = anal_options[!anal_options %in% test], selected = ctrl)
  })

  shiny::observeEvent(input$submit_integration, {
    test <- input$test_integration
    ctrl <- input$ctrl_integration
    anal_name <- input$integration_name
    anal_options <- anal_options_r()

    error_msg <- validate_integration(test, ctrl, anal_name, anal_options)

    if (is.null(error_msg)) {
      # clear error and disable button
      shinyjs::removeClass(selector = '#integration-form .validate-wrapper', class = 'has-error')
      shinyjs::disable('submit_integration')

      # run integration
      integrate_saved_scseqs(data_dir, test, ctrl, anal_name)

      # re-enable, clear inputs, close, and trigger update of available/selected anal
      shinyjs::enable('submit_integration')
      shiny::updateSelectizeInput('test_integration', selected = NULL)
      shiny::updateSelectizeInput('ctrl_integration', selected = NULL)
      shiny::updateTextInput('integration_name', value = NULL)
      new_anal_rv(anal_name)

    } else {
      # show error message
      shinyjs::html(selector = '#integration-form .validate-wrapper .help-block', html = error_msg)
      shinyjs::addClass(selector = '#integration-form .validate-wrapper', class = 'has-error')
    }

  })


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


  

