#' Explore Single Cell Clusters
#'
#' @param scseq \code{SingleCellExperiment} or \code{Seurat} object
#' @param markers Named list of \code{data.frame}s where \code{row.names} are the marker genes. One list per cluster in \code{scseq}
#'  with the same name as the cluster.
#' @param assay.type The assay data that gets used by \code{\link[scater]{plotTSNE}}. If \code{scseq} is a \code{Seurat}
#' object then the default ('logcounts') will be from the \code{data} slot of the \code{Seurat::DefaultAssay}.
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

explore_scseq_clusters <- function(scseq, markers = NULL, assay.type = 'logcounts') {

  # setup ----
  if (is.null(markers))
    markers <- get_scseq_markers(scseq, assay.type)

  # name cluster choices for drop down
  cluster_choices <- names(markers)
  names(cluster_choices) <- paste('Cluster', names(markers))

  # placeholder for default cluster to compare to
  default_contrast <- 'all'
  names(default_contrast) <- 'All other clusters'

  #  user interface ----

  ui <- miniUI::miniPage(
    shiny::tags$head(
      shiny::tags$style("#genecards {border-top-left-radius: 0; border-bottom-left-radius: 0; position: relative; z-index: 2; margin-left: -8px; margin-top: -26px}"),
      shiny::tags$style("button.btn.btn-default.dropdown-toggle {background-color: transparent !important;border-top-right-radius: 0; border-bottom-right-radius: 0}")
    ),
    miniUI::gadgetTitleBar("Explore Single-Cell Clusters"),
    miniUI::miniContentPanel(
      shiny::fluidRow(
        shiny::column(6,
                      shiny::uiOutput('groups_toggle'),
                      shiny::selectizeInput("selected_cluster", 'Show marker genes for:', choices = cluster_choices, width = '291.5px'),
                      shiny::selectizeInput("contrast_cluster", 'In comparison to:', choices = NULL, width = '291.5px'),
                      shiny::br(),
                      shiny::tags$div(
                        shiny::tags$div(style = "display:inline-block;", id='gene-container',
                                        shiny::selectizeInput('gene', 'Show expression for:', choices = NULL, width = '250px')
                                        ),
                        shinyBS::bsButton('genecards', label = '', icon = shiny::icon('external-link-alt', 'fa-fw'), style='default', title = 'Go to GeneCards')
                      )
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

    reseting_con <- reactiveVal(FALSE)

    # used to determine available groups to show (e.g. ctrl and test)
    available_groups_r <- shiny::reactive({
      groups <- unique(as.character(scseq$orig.ident))
      return(groups)
    })

    # used to determine cell groups to show (e.g. ctrl and test)
    selected_groups_r <- shiny::reactive({

      groups <- available_groups_r()
      # always show when just a single group
      if (length(groups) == 1) return(groups)

      return(input$groups)
    })

    # make it so that selected cluster can't contrast to itself
    contrast_choices_r <- shiny::reactive({
      contrast_choices <- c(default_contrast, cluster_choices)
      contrast_choices <- contrast_choices[contrast_choices != input$selected_cluster]
      return(contrast_choices)
    })


    # change cluster to compare against
    shiny::observeEvent(input$contrast_cluster, {

      # default contrast is all (named as selected contrast)
      id2 <- input$contrast_cluster
      id1 <- contrast <- input$selected_cluster
      if (id2 == '') return(NULL)

      if (input$contrast_cluster != 'all') {

        contrast <- paste0(id1, '-', id2)
        if (!contrast %in% names(markers))
          markers[[contrast]] <<- get_scseq_markers(scseq, ident.1 = id1, ident.2 = id2)
      }

      cluster_markers <- markers[[contrast]]
      choices <- row.names(cluster_markers)
      shiny::updateSelectizeInput(session, 'gene', choices = choices, server = TRUE)
    })


    # update marker genes based on cluster selection -----
    shiny::observe({
      cluster_markers <- markers[[input$selected_cluster]]
      choices <- row.names(cluster_markers)
      shiny::updateSelectizeInput(session, 'gene', choices = choices, server = TRUE)

      # reset contrast selections
      shiny::updateSelectizeInput(session, 'contrast_cluster', choices = contrast_choices_r())

    })



    # toggle for e.g. showing ctrl and/or test cells ----
    output$groups_toggle <- shiny::renderUI({
      # if more than one group allow showing cells based on groups
      groups <- available_groups_r()
      groups_toggle <- NULL
      if (length(groups) > 1) {
        groups_toggle <- shiny::tags$div(
          shinyWidgets::checkboxGroupButtons("groups", "Show cells for group:", choices = groups, selected = groups),
          shiny::br()
        )
      }
      return(groups_toggle)
    })


    # plots ------
    # show tSNE plot coloured by expression values
    output$marker_plot <- shiny::renderPlot({
      if (input$gene == '') return(NULL)
      plot_umap_gene(scseq, input$gene, selected_groups = selected_groups_r())
    })


    # show plot of predicted cell clusters
    output$cluster_plot <- shiny::renderPlot({
      legend_title <-'Cluster'
      plot_umap_cluster(scseq, selected_groups_r())
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

  }
  shiny::runGadget(shiny::shinyApp(ui, server), viewer = shiny::browserViewer())
}
