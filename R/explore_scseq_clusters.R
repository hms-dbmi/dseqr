

#' Explore Single Cell Clusters
#'
#' @param sce \code{SingleCellExperiment}
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' # import alevin quants
#' data_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'
#' sce <- load_scseq(data_dir)
#'
#' # subset by alevin whitelist norm/stabilize using good cell only
#' sce <- sce[, sce$whitelist]
#' sce <- norm_scseq(sce)
#' sce <- stabilize_scseq(sce)
#'
#' # get clusters and run tSNE
#' sce <- add_scseq_clusters(sce)
#'
#' set.seed(1000)
#' sce <- scater::runTSNE(sce, use_dimred="PCA")
#'
#' explore_scseq_clusters(sce)
#'

explore_scseq_clusters <- function(sce) {

  # setup ----

  # only upregulated as more useful for positive id of cell type
  markers <- scran::findMarkers(sce, clusters=sce$cluster, direction="up")

  cluster_choices <- names(markers)
  prev_gene <- NULL

  #  user interface ----

  ui <- miniUI::miniPage(
    shinyjs::useShinyjs(),
    # title bar
    miniUI::gadgetTitleBar("Explore Single-Cell Clusters"),
    miniUI::miniContentPanel(
      shiny::fillCol(flex = c(NA, NA, 1),
                     shiny::tags$div(
                       shinyWidgets::radioGroupButtons(
                         inputId = "cluster",
                         label = "Sort genes based on cluster:",
                         choices = cluster_choices
                       ),
                       shiny::selectizeInput('gene', 'Select gene:', choices = NULL)
                     ),
                     shiny::hr(),
                     shiny::fluidRow(
                       shiny::splitLayout(cellWidths = c("50%", "50%"),
                                          shiny::plotOutput("marker_plot"),
                                          shiny::plotOutput("cluster_plot"))
                     )
      )
    )
  )

  # server ----

  server <- function(input, output, session) {

    output$marker_plot <- shiny::renderPlot({
      # incase remove choice
      gene <- input$gene
      if (gene != '') {
        prev_gene <<- gene
      } else {
        gene <- prev_gene
      }

      if (is.null(gene)) return(NULL)

      suppressMessages(scater::plotTSNE(sce, colour_by = gene, point_size = 3, point_alpha = 1, theme_size = 14) +
                         ggplot2::scale_fill_distiller(palette = 'YlOrRd', name = '', direction = 1))



    })
    output$cluster_plot <- shiny::renderPlot({
      scater::plotTSNE(sce, colour_by = "cluster",  point_size = 3, point_alpha = 1, theme_size = 14)
    })

    shiny::observeEvent(input$cluster, {
      shiny::updateSelectizeInput(session, 'gene', choices = row.names(markers[[input$cluster]]), server = TRUE)
    })

  }


  shiny::runGadget(shiny::shinyApp(ui, server), viewer = shiny::browserViewer())


}
