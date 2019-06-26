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

  test_anal <- readRDS('data-raw/single-cell/example-data/test_lung_anal.rds')
  levels(test_anal$scseq$seurat_clusters) <- test_anal$annot
  Idents(test_anal$scseq) <- test_anal$scseq$seurat_clusters
  names(test_anal$markers) <- test_anal$annot

  scseq <- test_anal$scseq
  markers <- test_anal$markers

  if (Seurat::DefaultAssay(scseq) == 'integrated') Seurat::DefaultAssay(scseq) <- 'SCT'

  # setup ----
  if (is.null(markers))
    markers <- get_scseq_markers(scseq)

  # name cluster choices for drop down
  cluster_choices <- names(markers)

  test_cluster <- NULL


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
                      shiny::selectInput('selected_anal', 'Selected sample:', choice = c('sjia_lung', 'sjia_bm', 'sjia_pbmc'), width = '392.5px'),
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

    # reactive objects -----
    rvs <- shiny::reactiveValues(renaming = 0)

    # used to determine available groups to show (e.g. ctrl and test)
    available_groups_r <- shiny::reactive({
      groups <- unique(as.character(scseq$orig.ident))
      return(groups)
    })

    # used to determine cell groups to show (e.g. ctrl and test)
    selected_groups_r <- shiny::reactive({

      groups <- available_groups_r()
      # always show when just a single group
      if (length(groups) == 1 || input$groups == 'all') return(groups)

      return(input$groups)
    })

    # update cluster choices based on show_contrasts toggle and renaming
    cluster_choices_r <- shiny::reactive({


      if (isTRUE(input$rename_cluster)) {

        if (input$new_cluster_name != '') {
          # for ease remove any newly calculated markers
          markers  <<- markers[cluster_choices]

          # rename markers
          sel.clust <- shiny::isolate(input$selected_cluster)
          sel.idx   <- which(names(markers) == sel.clust)
          names(markers)[sel.idx] <<- input$new_cluster_name
          cluster_choices <<- names(markers)

          levels(scseq$seurat_clusters) <<- cluster_choices
          Seurat::Idents(scseq) <<- scseq$seurat_clusters

          shiny::isolate(rvs$renaming <- rvs$renaming + 1)
          test_cluster <<- input$new_cluster_name

        }

        # reset toggles to initial state
        shinyBS::updateButton(session, 'show_rename', value = FALSE)
        shinyBS::updateButton(session, 'rename_cluster', value = FALSE)
      }


      if (isTRUE(input$show_contrasts)) {
        # cluster choices are as compared to other clusters
        test <- shiny::isolate(input$selected_cluster)
        ctrls <- cluster_choices[cluster_choices != test]

        contrast_choices <- c(test, paste0(test, ' vs ', ctrls))
        # make sure not too long
        names(contrast_choices) <- stringr::str_trunc(paste0(test, ' vs ', c('all', ctrls)), 40)

        # update global so that can return to same cluster when toggle back
        test_cluster <<- test
        disable_rename <- TRUE

      } else {
        # cluster choices are the clusters themselves
        contrast_choices  <- cluster_choices
        disable_rename <- FALSE
      }

      # update icon on toggle
      icon <- ifelse(isTRUE(input$show_contrasts), 'chevron-down', 'chevron-right')
      shinyBS::updateButton(session, 'show_contrasts', icon = shiny::icon(icon, 'fa-fw'))
      shinyBS::updateButton(session, 'show_rename', disabled = disable_rename)

      return(contrast_choices)
    })

    # update cluster choices
    shiny::observe({
      # when first hit toggle, selected contrast doesn't change
      shiny::updateSelectizeInput(session, 'selected_cluster', choices = cluster_choices_r(), selected = test_cluster)
    })


    # update marker genes based on cluster selection -----
    shiny::observeEvent(input$selected_cluster, {
      sel <- input$selected_cluster
      if (sel == '') return(NULL)

      # get cluster if don't have (for comparing specific cluster)
      if (!sel %in% names(markers)) {
        con <- strsplit(sel, ' vs ')[[1]]
        markers[[sel]] <<- get_scseq_markers(scseq, ident.1 = con[1], ident.2 = con[2])
      }

      cluster_markers <- markers[[sel]]
      choices <- row.names(cluster_markers)
      choices <- c(choices, setdiff(row.names(scseq), choices))

      shiny::updateSelectizeInput(session, 'gene', choices = choices, server = TRUE)
    })



    # toggle for e.g. showing ctrl and/or test cells ----
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


    # conditional UI elements to either selected cluster or rename it ----
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
      plot_umap_gene(scseq, input$gene, selected_idents = selected_groups_r(), pt.size = pt.size)
    })


    # show plot of predicted cell clusters
    output$cluster_plot <- shiny::renderPlot({
      # update if renaming cluster
      rvs$renaming

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
