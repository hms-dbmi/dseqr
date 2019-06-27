#' Explore Single Cell Clusters
#'
#' @param scseq \code{SingleCellExperiment} or \code{Seurat} object
#' @param markers Named list of \code{data.frame}s where \code{row.names} are the marker genes. One list per cluster in \code{scseq}
#'  with the same name as the cluster.
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' # import kallisto quants
#' data_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'
#' scseq <- load_scseq(data_dir)
#'
#' # subset by whitelist norm/stabilize using good cell only
#' scseq <- scseq[, scseq$whitelist]
#' scseq <- preprocess_scseq(scseq)
#'
#' # get clusters and run tSNE
#' scseq <- add_scseq_clusters(scseq, resolution = 1.6)
#' scseq <- run_tsne(scseq)
#'
#' markers <- get_scseq_markers(scseq)
#'
#' explore_scseq_clusters(scseq, markers)
#'

explore_scseq_clusters <- function(data_dir, pt.size = 3) {

  anal_files <- list.files(data_dir)
  anal_options <- gsub('.rds$', '', anal_files)

  #  user interface ----

  ui <- miniUI::miniPage(
    shiny::tags$head(
      shiny::tags$style("#genecards, #show_contrasts, #show_rename, #rename_cluster {border-top-left-radius: 0; border-bottom-left-radius: 0; position: relative; z-index: 2; margin-left: -5px; margin-top: -26px;}"),
      shiny::tags$style("div.textinput-buttons #rename_cluster {margin-top: -3px;}"),
      shiny::tags$style("div.selectize-buttons .selectize-input, #show_contrasts, #new_cluster_name {border-top-right-radius: 0; border-bottom-right-radius: 0;}")
    ),
    miniUI::gadgetTitleBar("Explore Single-Cell Clusters", left = NULL, right = NULL),
    miniUI::miniContentPanel(
      shiny::fluidRow(
        shiny::column(6,
                      shiny::selectInput('selected_anal', 'Selected analysis:', choice = anal_options, width = '392.5px'),
                      shiny::br(),
                      shiny::conditionalPanel('input.show_rename == false',
                                              shiny::tags$div(class = 'selectize-buttons', style='height:79px;',
                                                              shiny::tags$div(style = "display:inline-block;",
                                                                              shiny::selectizeInput("selected_cluster", 'Show marker genes for:', choices = NULL, width = '306.5px')),
                                                              shinyBS::bsButton('show_contrasts', label = '', type = 'toggle', value = FALSE, icon = shiny::icon('chevron-right', 'fa-fw'), style='default', title = 'Toggle single group comparisons'),
                                                              shinyBS::bsButton('show_rename', label = '', type = 'toggle', value = FALSE, icon = shiny::icon('tag', 'fa-fw'), style='default', title = 'Toggle rename cluster'))
                      ),
                      shiny::conditionalPanel('input.show_rename == true', shiny::uiOutput('rename_ui')),
                      shiny::br(),
                      shiny::tags$div(class = 'selectize-buttons',
                                      shiny::tags$div(style = "display:inline-block;",
                                                      shiny::selectizeInput('gene', 'Show expression for:', choices = NULL, width = '350px')),
                                      shinyBS::bsButton('genecards', label = '', icon = shiny::icon('external-link-alt', 'fa-fw'), style='default', title = 'Go to GeneCards')
                      ),
                      shiny::br(),
                      shiny::uiOutput('groups_toggle')
        ),
        shiny::column(6,
                      shiny::plotOutput("cluster_plot")
        )
      ),
      shiny::hr(),
      shiny::tags$div(class = 'col-sm-6 col-lg-6 col-lg-push-6',
                      shiny::plotOutput("marker_plot")
      ),
      shiny::tags$div(class = 'col-sm-6 col-lg-6 col-lg-pull-6',
                      shiny::plotOutput('biogps')
      )
    )
  )

  # server ----

  server <- function(input, output, session) {

    # reactive values (can update and persist within session) -----
    annot_rv <- shiny::reactiveVal(NULL)
    con_markers_rv <- shiny::reactiveVal(list())


    # reactive expressions (auto update) ----------

    # groups to show (e.g. ctrl and test)
    available_groups_r <- shiny::reactive({
      scseq <- scseq_r()
      groups <- unique(as.character(scseq$orig.ident))
      return(groups)
    })

    # current analysis
    anal_r <- shiny::reactive({
      selected_file <- paste0(input$selected_anal, '.rds')

      anal <- readRDS(file.path(data_dir, selected_file))

      if (!is.null(anal$annot)) {
        levels(anal$scseq$seurat_clusters) <- anal$annot
        Idents(anal$scseq) <- anal$scseq$seurat_clusters
      }

      if (Seurat::DefaultAssay(anal$scseq) == 'integrated')
        Seurat::DefaultAssay(anal$scseq) <- 'SCT'

      #  annot
      annot_rv(anal$annot)
      return(anal)
    })


    # the markers
    markers_r <- shiny::reactive({
      anal <- anal_r()
      annot <- annot_rv()
      markers <- anal$markers

      if (!is.null(annot)) {
        names(markers) <- annot
      }

      return(markers)
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
      return(test_cluster)
    })

    # used to determine cell groups to show (e.g. ctrl and test)
    selected_groups_r <- shiny::reactive({

      groups <- available_groups_r()
      # always show when just a single group
      if (length(groups) == 1 || input$groups == 'all') return(groups)

      return(input$groups)
    })



    # update cluster/contrast choices based on show_contrasts toggle and renaming
    group_choices_r <- shiny::reactive({

      clusters <- clusters_r()

      if (isTRUE(input$show_contrasts)) {
        # group choices are as compared to other clusters
        test <- shiny::isolate(test_cluster_r())
        ctrls <- clusters[clusters != test]

        contrast_choices <- c(test, paste0(test, ' vs ', ctrls))
        # make sure not too long
        names(contrast_choices) <- stringr::str_trunc(paste0(test, ' vs ', c('all', ctrls)), 40)

      } else {
        # cluster choices are the clusters themselves
        contrast_choices  <- clusters
      }

      return(contrast_choices)
    })


    # observations (do stuff if something changes) -------
    # update annot if rename a cluster
    shiny::observeEvent(input$rename_cluster, {

      if (isTRUE(input$rename_cluster)) {

        if (input$new_cluster_name != '') {

          # update reactive annotation
          annot <- annot_rv()
          sel.clust <- input$selected_cluster
          sel.idx   <- which(annot == sel.clust)
          annot[sel.idx] <- input$new_cluster_name
          annot_rv(annot)
        }

        # reset toggles to initial state
        shinyBS::updateButton(session, 'show_rename', value = FALSE)
        shinyBS::updateButton(session, 'rename_cluster', value = FALSE)
      }
    })


    # update group choices if they change
    shiny::observe({
      shiny::updateSelectizeInput(session, 'selected_cluster', choices = group_choices_r())
    })

    # update group choices/buttons if show contrasts is toggled
    shiny::observeEvent(input$show_contrasts, {
      # update icon on toggle
      if (input$show_contrasts) {
        disable_rename <- TRUE
        icon <- 'chevron-down'
        selected <- NULL

      } else {
        disable_rename <- FALSE
        icon <- 'chevron-right'
        selected <- test_cluster_r()
      }
      shiny::updateSelectizeInput(session, 'selected_cluster', choices = group_choices_r(), selected = selected)
      shinyBS::updateButton(session, 'show_contrasts', icon = shiny::icon(icon, 'fa-fw'))
      shinyBS::updateButton(session, 'show_rename', disabled = disable_rename)
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

      shiny::updateSelectizeInput(session, 'gene', choices = choices, server = TRUE)
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
    output$rename_ui <- shiny::renderUI({
      current_name <- input$selected_cluster
      shiny::tags$div(class = 'textinput-buttons', style='height:79px;',
                      shiny::tags$div(style = "display:inline-block;",
                                      shiny::textInput("new_cluster_name", 'New cluster name:', width = '349px', placeholder = paste('Type new name for', current_name, '...'))),
                      shinyBS::bsButton('rename_cluster', label = '', type = 'toggle', value = FALSE, icon = shiny::icon('plus', 'fa-fw'), style='default', title = 'Rename cluster'))
    })



    # plots ------
    # show tSNE plot coloured by expression values
    output$marker_plot <- shiny::renderPlot({

      if (input$gene == '') return(NULL)
      plot_umap_gene(scseq_r(), input$gene, selected_idents = selected_groups_r(), pt.size = pt.size)
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
  shiny::runGadget(shiny::shinyApp(ui, server), viewer = shiny::browserViewer())
}
