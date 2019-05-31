

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
                         ggplot2::scale_fill_distiller(palette = 'YlOrRd', name = gene, direction = 1))



    })
    output$cluster_plot <- shiny::renderPlot({

      # make selected cluster stand out
      point_alpha <- rep(0.1, ncol(sce))
      point_alpha[sce$cluster == input$cluster] <- 1

      scater::plotTSNE(sce, colour_by = "cluster",  point_size = 3, point_alpha = point_alpha, theme_size = 14) +
        ggplot2::theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    })


    gene_link <- shiny::reactive({
      # incase remove choice
      gene <- input$gene
      if (gene != '') {
        prev_gene <<- gene
      } else {
        gene <- prev_gene
      }

      enid <- biogps[gene, ENTREZID]
      paste0('http://genewiki.sulab.org/map/wiki/', enid, '/')
    })

    shiny::observeEvent(input$wiki, {
      utils::browseURL(gene_link())
    })

    # plot BioGps plot -----
    output$biogps <- shiny::renderPlot({

      # incase remove choice
      gene <- input$gene
      if (gene != '') {
        prev_gene <<- gene
      } else {
        gene <- prev_gene
      }

      if (is.null(gene) || !gene %in% biogps[, SYMBOL]) return(NULL)

      gene_dat <- unlist(biogps[gene, -c('ENTREZID', 'SYMBOL')])
      gene_dat <- sort(gene_dat, decreasing = TRUE)[1:20]
      gene_dat <- tibble::tibble(mean = gene_dat,
                                 source = factor(names(gene_dat), levels = rev(names(gene_dat))))

      ggplot(gene_dat, aes(x = source, y = mean, fill = mean)) +
        theme_minimal() +
        geom_bar(stat = "identity", color = 'black', size = 0.1, width = 0.7) +
        scale_fill_distiller(palette = 'YlOrRd', name = '', direction = 1) +
        xlab('') +
        ylab('') +
        ggtitle('BioGPS Human Gene Atlas Expression') +
        scale_y_continuous(expand = c(0, 0)) + # Set the axes to cross at 0
        coord_flip() +
        theme(panel.grid = element_blank(),
              panel.border = element_blank(),
              # axis.ticks.y = element_line(size = 0.1),
              plot.title = element_text(),
              axis.text = element_text(size = 12),
              axis.text.x = element_blank(),
              legend.position = "none")


    })

    shiny::observeEvent(input$cluster, {
      shiny::updateSelectizeInput(session, 'gene', choices = row.names(markers[[input$cluster]]), server = TRUE)
    })

  }


  shiny::runGadget(shiny::shinyApp(ui, server), viewer = shiny::browserViewer())


}
