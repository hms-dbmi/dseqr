  # current analysis
  anal_r_contrasts <- shiny::eventReactive(input$selected_anal_contrasts, {

    # load analysis parts
    anal_name <- input$selected_anal_contrasts
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


   # analysis options
  shiny::observe({
    shiny::updateSelectizeInput(session, 'selected_anal_contrasts', choices = anal_options_r(), selected = 'lung_sjia')
  })


    # show plot of predicted cell clusters
  output$cluster_plot_contrasts <- shiny::renderPlot({
    scseq <- scseq_r()
    plot_umap_cluster(scseq, pt.size = input$point_size_contrasts)
  })