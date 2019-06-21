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
#' # import alevin quants
#' data_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'
#' scseq <- load_scseq(data_dir, type = 'Seurat')
#'
#' # subset by alevin whitelist norm/stabilize using good cell only
#' scseq <- scseq[, scseq$whitelist]
#' scseq <- preprocess_scseq(scseq)
#'
#' # get clusters and run tSNE
#' scseq <- add_scseq_clusters(scseq)
#' scseq <- run_tsne(scseq)
#'
#' explore_scseq_clusters(scseq)
#'

explore_scseq_clusters <- function(scseq, markers = NULL, assay.type = 'logcounts') {

  # setup ----
  if (is.null(markers))
    markers <- get_scseq_markers(scseq, assay.type)

  # plots all based on SingleCellExperiment
  # use non-integrated but corrected assay for marker gene plots
  sce <- srt_to_sce(scseq, "SCT")

  prev_gene <- NULL
  cluster_choices <- names(markers)
  names(cluster_choices) <- paste('Cluster', names(markers))

  #  user interface ----

  ui <- miniUI::miniPage(
    shiny::tags$head(
      shiny::tags$style("#genecards, #subcluster {border-top-left-radius: 0; border-bottom-left-radius: 0; position: relative; z-index: 2; margin-left: -6px;}"),
      shiny::tags$style("button.btn.btn-default.dropdown-toggle {background-color: transparent !important;border-top-right-radius: 0; border-bottom-right-radius: 0}"),
      shiny::tags$style(".sub-chev {font-size: 10px; padding: 0 10px; color: darkgray;}")
    ),
    miniUI::gadgetTitleBar("Explore Single-Cell Clusters"),
    miniUI::miniContentPanel(
      shiny::fluidRow(
        shiny::column(6,
                      shiny::uiOutput('groups_toggle'),
                      shiny::selectizeInput("selected_cluster", 'Show marker genes for:', choices = cluster_choices, width = '291.5px'),
                      shiny::br(),
                      shiny::tags$div(
                        shiny::tags$div(style = "display:inline-block;", id='gene-container', shinyWidgets::pickerInput('gene', 'Show expression for:', choices = NULL, width = '250px', options = shinyWidgets::pickerOptions(liveSearch = TRUE))),
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

    # used to determine available groups to show (e.g. ctrl and test)
    available_groups_r <- shiny::reactive({
      groups <- unique(as.character(sce$orig.ident))
      return(groups)
    })

    # used to determine cell groups to show (e.g. ctrl and test)
    selected_groups_r <- shiny::reactive({

      groups <- available_groups_r()
      # always show when just a single group
      if (length(groups) == 1) return(groups)

      return(input$groups)
    })


    # update marker genes based on cluster/subcluster selection -----
    shiny::observe({
      cluster_markers <- markers[[input$selected_cluster]]
      choices <- row.names(cluster_markers)
      shinyWidgets::updatePickerInput(session, 'gene', choices = choices)
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
      if (is.null(input$gene)) return(NULL)
      plot_tsne_gene(sce, input$gene,selected_groups = selected_groups_r(), assay.type = assay.type)
    })


    # show plot of predicted cell clusters
    output$cluster_plot <- shiny::renderPlot({
      legend_title <-'Cluster'
      plot_tsne_cluster(sce, legend_title, selected_groups_r())
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
