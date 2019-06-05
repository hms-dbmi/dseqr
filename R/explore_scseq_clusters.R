

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

  # if more than one group allow showing cells based on groups
  groups <- unique(sce$orig.ident)
  groups_toggle <- NULL
  if (length(groups) > 1) {
    groups_toggle <- shiny::tags$div(
      shinyWidgets::checkboxGroupButtons("groups", "Show cells for:", choices = groups, selected = groups),
      shiny::br()
    )
  }

  current_markers <- markers
  subcluster_markers <- list()
  cluster_choices <- names(markers)
  names(cluster_choices) <- paste('Cluster', cluster_choices)
  prev_gene <- NULL
  prev_cluster <- NULL

  #  user interface ----

  ui <- miniUI::miniPage(
    shiny::tags$head(
      shiny::tags$style("#wiki, #subcluster {border-top-left-radius: 0; border-bottom-left-radius: 0; position: relative; z-index: 2; margin-left: -6px;}"),
      shiny::tags$style("button.btn.btn-default.dropdown-toggle {background-color: transparent !important;border-top-right-radius: 0; border-bottom-right-radius: 0}"),
      shiny::tags$style(".sub-chev {font-size: 10px; padding: 0 10px; color: darkgray;}")
    ),
    miniUI::gadgetTitleBar("Explore Single-Cell Clusters"),
    miniUI::miniContentPanel(
      shiny::fluidRow(
        shiny::column(6,
                      groups_toggle,
                      shiny::uiOutput('subcluster_ui'),
                      shiny::br(),
                      shiny::tags$div(
                        shiny::tags$div(style = "display:inline-block;", id='gene-container', shinyWidgets::pickerInput('gene', 'Show expression for:', choices = NULL, width = '250px')),
                        shinyBS::bsButton('wiki', label = '', icon = shiny::icon('external-link-alt', 'fa-fw'), style='default', title = 'Go to Wikipedia')
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

    # selected gene: use previously selected if delete selection ----
    gene <- shiny::reactive({
      gene <- input$gene
      if (is.null(gene) || gene == '') {
        gene <- prev_gene
      } else {
        prev_gene <<- gene
      }
      return(gene)
    })

    groups_r <- shiny::reactive({
      # always show when just a single group
      if (length(groups) == 1) return(groups)

      return(input$groups)
    })

    output$subcluster_ui <- shiny::renderUI({

      if (subcluster()) {
        title <- 'show all clusters'
        icon <- shiny::icon('chevron-left', 'fa-fw')
        selected <- NULL
      } else {
        title <- 'explore this cluster'
        icon <- shiny::icon('chevron-right', 'fa-fw')
        selected <- cluster()
      }

      shiny::tags$div(
        shiny::tags$div(style = "display:inline-block; text-overflow:",
                        shinyWidgets::pickerInput("cluster", 'Show marker genes for:',
                                                  choices = cluster_choices, width = '250px', selected = selected)),
        shinyBS::bsButton('subcluster', label = '', value = subcluster(), type = 'toggle', icon = icon, title = title)
      )
    })

    # show tSNE plot coloured by expression values -----

    output$marker_plot <- shiny::renderPlot({

      # show selected cluster only
      if (subcluster())
        sce <- sce[, sce$cluster == cluster()]

      gene <- gene()
      if (is.null(gene)) return(NULL)

      # make selected groups stand out
      point_alpha <- rep(1, ncol(sce))
      point_alpha[!sce$orig.ident %in% groups_r()] <- 0.1

      suppressMessages(scater::plotTSNE(sce, by_exprs_values = assay.type, colour_by = gene, point_size = 3, point_alpha = point_alpha, theme_size = 14) +
                         ggplot2::scale_fill_distiller(palette = 'YlOrRd', name = gene, direction = 1))

    })

    subcluster <- shiny::reactive({
      if (!length(input$subcluster)) return(FALSE)
      return(input$subcluster)
    })

    cluster <- shiny::reactive({
      cluster <- input$cluster
      if (is.null(cluster) || (!subcluster() && !cluster %in% sce$orig.ident))
        prev_cluster <<- cluster

      return(prev_cluster)
    })

    # show plot of predicted cell clusters -----
    output$cluster_plot <- shiny::renderPlot({

      # default
      colour_by <- 'cluster'

      # show selected cluster only and color by group
      if (subcluster()) {
        sce <- sce[, sce$cluster == cluster()]
        colour_by <- 'orig.ident'
      }

      # make selected cluster stand out and groups
      point_alpha <- rep(0.1, ncol(sce))
      point_alpha[sce$cluster == cluster()] <- 1
      point_alpha[!sce$orig.ident %in% groups_r()] <- 0.1

      legend_title <- ifelse(subcluster(), 'Group', 'Cluster')

      scater::plotTSNE(sce, by_exprs_values = assay.type, colour_by = colour_by,  point_size = 3, point_alpha = point_alpha, theme_size = 14) +
        ggplot2::guides(fill = ggplot2::guide_legend(legend_title)) +
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
      if (!length(gene) || !gene %in% biogps[, SYMBOL]) return(NULL)

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

    # Change cluster/subcluster for genes -----
    shiny::observeEvent(subcluster(), {

      # if subcluster get subcluster markers
      if (subcluster()) {
        if (!cluster() %in% names(subcluster_markers)) {
          scseq <- scseq[, scseq$seurat_clusters == cluster()]
          Seurat::Idents(scseq) <- scseq$orig.ident
          subcluster_markers[[cluster()]] <<- get_scseq_markers(scseq)
        }
        current_markers <<- subcluster_markers[[cluster()]]
        cluster_choices <- names(current_markers)
        names(cluster_choices) <- cluster_choices
        choicesOpt <- list(content = paste0('<span><span style="color: darkgray;">', 'Cluster ', cluster(), ' </span>', shiny::icon('chevron-right', 'sub-chev'), cluster_choices, '</span>'))
        selected <- NULL

      } else {
        choicesOpt <- NULL
        current_markers <<- markers
        selected <- cluster()
      }

      shinyWidgets::updatePickerInput(session, 'cluster', choices = cluster_choices, selected = selected, choicesOpt = choicesOpt)
    })

    shiny::observeEvent(input$cluster, {
      shinyWidgets::updatePickerInput(session, 'gene', choices = row.names(current_markers[[input$cluster]]))
    })

  }
  shiny::runGadget(shiny::shinyApp(ui, server), viewer = shiny::browserViewer())
}
