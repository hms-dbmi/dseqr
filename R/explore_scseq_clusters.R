

#' Explore Single Cell Clusters
#'
#' @param sce \code{SingleCellExperiment}
#'
#' @return
#' @export
#'
#' @examples
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
  biogps <- readRDS(system.file('extdata', 'biogps.rds', package = 'drugseqr'))

  # only upregulated as more useful for positive id of cell type
  markers <- scran::findMarkers(sce, clusters=sce$cluster, direction="up")

  cluster_choices <- names(markers)
  prev_gene <- NULL

  #  user interface ----

  ui <- miniUI::miniPage(
    shiny::tags$head(
      shiny::tags$style("#wiki {margin-top: -26px; border-top-left-radius: 0; border-bottom-left-radius: 0; position: relative; z-index: 2; margin-left: -8px;}")
    ),
    miniUI::gadgetTitleBar("Explore Single-Cell Clusters"),
    miniUI::miniContentPanel(
      shiny::splitLayout(cellWidths = c("45%", "50%"),
                         shiny::tags$div(
                           style = 'height: 400px',
                           shinyWidgets::radioGroupButtons("cluster", "Sort genes based on cluster:", choices = cluster_choices),
                           shiny::tags$div(
                             shiny::tags$div(style = "display:inline-block", shiny::selectizeInput('gene', 'Select gene:', choices = NULL)),
                             shinyBS::bsButton('wiki', label = '', icon = shiny::icon('external-link-alt'), style='default', title = 'Go to Wikipedia')
                         )),
                         shiny::plotOutput("cluster_plot", width = '99%')
      ),
      shiny::hr(),
      shiny::fluidRow(
        shiny::splitLayout(cellWidths = c("45%", "50%"),
                           shiny::plotOutput('biogps'),
                           shiny::plotOutput("marker_plot")
        )
      )
    )
  )

  # server ----

  server <- function(input, output, session) {

    # selected gene: use previously selected if delete selection ----
    gene <- shiny::reactive({
      gene <- input$gene
      if (gene == '') {
        gene <- prev_gene
      } else {
        prev_gene <<- gene
      }

      return(gene)
    })

    # show tSNE plot coloured by expression values -----

    output$marker_plot <- shiny::renderPlot({

      gene <- gene()
      if (is.null(gene)) return(NULL)

      suppressMessages(scater::plotTSNE(sce, colour_by = gene, point_size = 3, point_alpha = 1, theme_size = 14) +
                         ggplot2::scale_fill_distiller(palette = 'YlOrRd', name = gene, direction = 1))



    })

    # show plot of predicted cell clusters -----
    output$cluster_plot <- shiny::renderPlot({

      # make selected cluster stand out
      point_alpha <- rep(0.1, ncol(sce))
      point_alpha[sce$cluster == input$cluster] <- 1

      scater::plotTSNE(sce, colour_by = "cluster",  point_size = 3, point_alpha = point_alpha, theme_size = 14) +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank())
    })

    # link to Wikipedia page for gene ----
    gene_link <- shiny::reactive({
      enid <- biogps[gene(), ENTREZID]
      return(paste0('http://genewiki.sulab.org/map/wiki/', enid, '/'))
    })

    # Click link out to Wikipedia ----
    shiny::observeEvent(input$wiki, {
      utils::browseURL(gene_link())
    })

    # plot BioGPS data -----
    output$biogps <- shiny::renderPlot({
      gene <- gene()
      if (is.null(gene) || !gene %in% biogps[, SYMBOL]) return(NULL)

      gene_dat <- unlist(biogps[gene, -c('ENTREZID', 'SYMBOL')])
      gene_dat <- sort(gene_dat, decreasing = TRUE)[1:20]
      gene_dat <- tibble::tibble(mean = gene_dat,
                                 source = factor(names(gene_dat), levels = rev(names(gene_dat))))

      ggplot2::ggplot(gene_dat, ggplot2::aes(x = source, y = mean, fill = mean)) +
        ggplot2::theme_minimal() +
        ggplot2::geom_bar(stat = "identity", color = 'black', size = 0.1, width = 0.7) +
        ggplot2::scale_fill_distiller(palette = 'YlOrRd', name = '', direction = 1) +
        ggplot2::xlab('') +
        ggplot2::ylab('') +
        ggplot2::ggtitle('BioGPS Human Gene Atlas Expression') +
        ggplot2::scale_y_continuous(expand = c(0, 0)) + # Set the axes to cross at 0
        ggplot2::coord_flip() +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
              panel.border = ggplot2::element_blank(),
              plot.title = ggplot2::element_text(),
              axis.text = ggplot2::element_text(size = 12),
              axis.text.x = ggplot2::element_blank(),
              legend.position = "none")


    })

    # Change cluster to sort genes -----

    shiny::observeEvent(input$cluster, {
      shiny::updateSelectizeInput(session, 'gene', choices = row.names(markers[[input$cluster]]), server = TRUE)
    })

  }


  shiny::runGadget(shiny::shinyApp(ui, server), viewer = shiny::browserViewer())


}
