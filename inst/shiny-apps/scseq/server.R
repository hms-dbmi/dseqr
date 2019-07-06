scPage <- function(input, output, session, data_dir, new_anal, new_annot) {


  # the analysis and options
  scForm <- callModule(scForm, 'form',
                       data_dir = data_dir,
                       new_anal = new_anal,
                       new_annot = new_annot)

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

scForm <- function(input, output, session, data_dir, new_anal, new_annot) {


  # the analysis and options
  scAnal <- callModule(selectedAnal, 'anal',
                       data_dir = data_dir,
                       new_anal = new_anal,
                       new_annot = new_annot)

  # integration stuff
  scIntegration <- callModule(integrationForm, 'integration',
                              anal_options = scAnal$anal_options,
                              show_integration = scAnal$show_integration)

  # the cluster
  scCluster <- callModule(selectedCluster, 'cluster',
                          scseq = scAnal$scseq,
                          markers = scAnal$markers,
                          annot = scAnal$annot)

  # the gene selection
  scGene <- callModule(selectedGene, 'gene',
                       selected_cluster = scCluster$selected_cluster,
                       scseq = scAnal$scseq,
                       markers = scAnal$markers,
                       con_markers = scCluster$con_markers)


  # the groups selection
  scGroups <- callModule(selectedGroups, 'groups', scseq = scAnal$scseq)


  return(list(
    scseq = scAnal$scseq,
    plot_styles = scAnal$plot_styles,
    selected_gene = scGene$selected_gene,
    selected_groups = scGroups$selected_groups
  ))
}

scBioGpsPlot <- function(input, output, session, selected_gene) {
  # plot BioGPS data
  output$biogps_plot <- renderPlot({
    plot_biogps(selected_gene())
  })
}


selectedGroups <- function(input, output, session, scseq) {

  # groups to show (e.g. ctrl and test)
  available_groups <- shiny::reactive({
    unique(as.character(scseq()$orig.ident))
  })

  selected_groups <- reactive({
    groups <- available_groups()
    # always show when just a single group
    if (length(groups) == 1 || input$selected_group == 'all') return(groups)

    return(input$selected_group)
  })


  return(list(
    selected_groups = selected_groups
  ))
}

selectedGene <- function(input, output, session, selected_cluster, scseq, markers, con_markers) {



  # update marker genes based on cluster selection
  gene_choices <- reactive({
    sel <- selected_cluster()
    req(sel)

    # allow selecting non-marker genes (at bottom of list)
    cluster_markers <- c(markers(), con_markers())[[sel]]
    choices <- row.names(cluster_markers)
    choices <- c(choices, setdiff(row.names( scseq()), choices))

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


  observe({
    updateSelectizeInput(session, 'selected_gene', choices = gene_choices(), selected = NULL, server = TRUE)
  })

  return(list(
      selected_gene = reactive(input$selected_gene)
  ))

}

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

integrationForm <- function(input, output, session, anal_options, show_integration) {

  # show/hide integration form
  observe({
    toggle(id = "integration-form", anim = TRUE, condition = show_integration())
    toggleClass(id = "show_integration", 'active', condition = show_integration())
  })

  integration_name <- reactive(input$integration_name)
  integration_options <- reactive(anal_options()$Individual)

  observe({
    ctrl <- input$ctrl_integration
    test <- isolate(input$test_integration)
    options <- integration_options()

    updateSelectizeInput(session, 'test_integration', choices = options[!options %in% ctrl], selected = test)
  })

  observe({
    ctrl <- isolate(input$ctrl_integration)
    test <- input$test_integration
    options <- integration_options()

    updateSelectizeInput(session, 'ctrl_integration', choices = options[!options %in% test], selected = ctrl)
  })

}

selectedCluster <- function(input, output, session, scseq, markers, annot) {


  con_markers <- reactiveVal(list())

  show_contrasts <- reactive({ input$show_contrasts %% 2 != 0 })

  selected_cluster <- reactive(input$selected_cluster)

  show_rename <- reactive({
    (input$rename_cluster + input$show_rename) %% 2 != 0
  })

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

  test_cluster <- reactive({
    test_cluster <- input$selected_cluster
    gsub(' vs .+?$', '', test_cluster)
  })

  # update cluster/contrast choices
  observe({
    scseq <- scseq()
    clusters <- annot()

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
                                     ctrlColor = c('white', colours[ctrls]))


    } else {
      # show the cell numbers/percentages
      ncells <- tabulate(scseq$seurat_clusters)
      pcells <- round(ncells / sum(ncells) * 100)
      pspace <- strrep('&nbsp;&nbsp;', 2 - nchar(pcells))

      # cluster choices are the clusters themselves
      testColor <- get_palette(clusters)
      contrast_choices <- data.frame(name = stringr::str_trunc(clusters, 27),
                                     value = clusters,
                                     label = clusters,
                                     testColor,
                                     ncells,
                                     pcells, pspace)
    }

    contrast_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))

    updateSelectizeInput(session, 'selected_cluster',
                         choices = contrast_choices,
                         options = contrast_options, server = TRUE)


  })

  # update ui for renaming a cluster
  observe({
    if (!show_rename())
      updateTextInput(session, 'new_cluster_name', value = '', placeholder = paste('Type new name for', input$selected_cluster, '...'))
  })

  # get cluster if don't have (for comparing specific cluster)
  observe({
    sel <- selected_cluster()
    con_markers <- con_markers()
    req(sel)


    if (!sel %in% names(c(con_markers, markers()))) {
      con <- strsplit(sel, ' vs ')[[1]]
      con_markers[[sel]] <- get_scseq_markers(scseq(), ident.1 = con[1], ident.2 = con[2])
      con_markers(con_markers)
    }
  })


  return(list(
    selected_cluster = selected_cluster,
    con_markers = con_markers

  ))
}


selectedAnal <- function(input, output, session, data_dir, new_anal, new_annot) {


  # load scseq
  scseq <- reactive({
    req(input$selected_anal)

    scseq_path <- scseq_part_path(data_dir, input$selected_anal, 'scseq')
    scseq <- readRDS(scseq_path)

    if (Seurat::DefaultAssay(scseq) == 'integrated')
      Seurat::DefaultAssay(scseq) <- 'SCT'

    return(scseq)
  })

  # load markers
  markers <- reactive({
    req(input$selected_anal)

    markers_path <- scseq_part_path(data_dir, input$selected_anal, 'markers')
    readRDS(markers_path)
  })

  # load annotation for clusters
  annot <- reactive({
    req(input$selected_anal)

    annot_path <- scseq_part_path(data_dir, input$selected_anal, 'annot')
    readRDS(annot_path)
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
    updateSelectizeInput(session, 'selected_anal', choices = anal_options(), selected = 'lung_sjia')
  })

  # get styles and integration info
  plot_styles <- callModule(plotStyles, 'styles')
  show_integration <- callModule(showIntegration, 'integration')


  # return anal and options to app
  return(list(
    scseq = scseq,
    markers = markers,
    annot = annot,
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
  return(reactive({
    input$show_integration %% 2 != 0
  }))
}


server <- function(input, output, session) {

  # get arguments from calling function
  # defaults for server
  data_dir <- getShinyOption('data_dir', '/srv/shiny-server/drugseqr/scseq/sjia')

  # things that the global app needs to know about
  rv <- reactiveValues(new_anal = FALSE, new_annot = FALSE)

  # scAnal <- callModule(selectedAnal, 'sc', data_dir, reactive(new_anal))

  print(data_dir)

  # single cell analysis and options
  scPage <- callModule(scPage, 'sc',
                       data_dir = data_dir,
                       new_anal = reactive(rv$new_anal),
                       new_annot = reactive(rv$new_annot))



  # reactive values (can update and persist within session) -----
  # new_anal_rv <- reactiveVal(NULL)

  # # reactive expressions (auto update) ----------

  # # available analyses
  # anal_options_r <- reactive({
  #   browser()
  #   # reactive to new anals
  #   new_anal_rv()

  #   # make sure integrated rds exists
  #   int_path <- file.path(data_dir, 'integrated.rds')
  #   if (!file.exists(int_path)) saveRDS(NULL, int_path)

  #   # use saved anals as options
  #   int_options <- readRDS(file.path(data_dir, 'integrated.rds'))
  #   ind_options <- setdiff(list.files(data_dir), c(int_options, 'integrated.rds'))

  #   # must be a list if length one for option groups to work
  #   if (length(int_options) == 1) int_options <- list(int_options)
  #   if (length(ind_options) == 1) ind_options <- list(ind_options)

  #   anal_options <- list(Individual = ind_options, Integrated = int_options)

  #   return(anal_options)
  # })



  # source('server/single-cell.R', local = TRUE)
  # source('server/contrasts.R', local = TRUE)


}
