

#' Explore Single Cell Clusters
#'
#' @param scseq \code{SingleCellExperiment} or \code{Seurat} object
#'
#' @return
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

explore_scseq_clusters <- function(scseq, markers = NULL, assay.type = 'logcounts', use_dimred = 'TSNE') {

  # setup ----
  biogps <- readRDS(system.file('extdata', 'biogps.rds', package = 'drugseqr'))

  if (is.null(markers))
    markers <- get_scseq_markers(scseq, assay.type)

  # plots all based on SingleCellExperiment
  if (class(scseq) == 'Seurat') {
    sce <- srt_to_sce(scseq)

  } else if (class(scseq) == 'SingleCellExperiment') {
    sce <- scseq

  } else {
    stop('scseq must be either a Seurat or SingleCellExperiment object.')
  }

  # if more than one sample allow showing cells based on sample
  samples <- unique(sce$orig.ident)
  samples_toggle <- NULL
  if (length(samples) > 1) {
    samples_toggle <- shinyWidgets::checkboxGroupButtons("samples", "Highlight samples:", choices = samples, selected = samples)
  }

  cluster_choices <- names(markers)
  names(cluster_choices) <- paste('Cluster', cluster_choices)
  prev_gene <- NULL

  #  user interface ----

  ui <- miniUI::miniPage(
    shiny::tags$head(
      shiny::tags$style("#wiki {margin-top: -26px; border-top-left-radius: 0; border-bottom-left-radius: 0; position: relative; z-index: 2; margin-left: -16px;}"),
      shiny::tags$style("#gene-container .selectize-input {border-top-left-radius: 0; border-bottom-left-radius: 0; margin-left: -8px;}"),
      shiny::tags$style("#gene-container .selectize-dropdown {margin-left: -8px;}")
    ),
    miniUI::gadgetTitleBar("Explore Single-Cell Clusters"),
    miniUI::miniContentPanel(
      shiny::fluidRow(
        shiny::column(6,
                      shiny::tags$div(
                        samples_toggle,
                        shiny::tags$div(
                          shiny::tags$div(style = "display:inline-block; text-overflow:", shiny::selectInput("cluster", 'Cluster and marker gene:', choices = cluster_choices, width = '200px')),
                          shiny::tags$div(style = "display:inline-block;", id='gene-container', shiny::selectizeInput('gene', '', choices = NULL, width = '200px')),
                          shinyBS::bsButton('wiki', label = '', icon = shiny::icon('external-link-alt'), style='default', title = 'Go to Wikipedia')
                        ))

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

      suppressMessages(scater::plotTSNE(sce, by_exprs_values = assay.type, colour_by = gene, point_size = 3, point_alpha = 1, theme_size = 14) +
                         ggplot2::scale_fill_distiller(palette = 'YlOrRd', name = gene, direction = 1) +
                         ggplot2::theme(legend.position = 'none'))


    })

    # show plot of predicted cell clusters -----
    output$cluster_plot <- shiny::renderPlot({

      # make selected cluster stand out
      point_alpha <- rep(0.1, ncol(sce))
      point_alpha[sce$cluster == input$cluster] <- 1
      point_alpha[!sce$orig.ident %in% input$samples] <- 0.1


      scater::plotTSNE(sce, by_exprs_values = assay.type, colour_by = "cluster",  point_size = 3, point_alpha = point_alpha, theme_size = 14) +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank())

    })

    # link to Wikipedia page for gene ----
    gene_link <- shiny::reactive({
      enid <- biogps[gene(), ENTREZID]
      if (is.na(enid))
        return(paste0('https://en.wikipedia.org/wiki/', gene()))
      else
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
