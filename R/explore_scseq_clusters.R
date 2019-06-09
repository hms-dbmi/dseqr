#' Explore Single Cell Clusters
#'
#' @param scseq \code{SingleCellExperiment} or \code{Seurat} object
#' @param markers Named list of \code{data.frame}s where \code{row.names} are the marker genes. One list per cluster in \code{scseq}
#'  with the same name as the cluster.
#' @param assay.type The assay data that gets used by \code{\link[scater]{plotTSNE}}. If \code{scseq} is a \code{Seurat}
#' object then the default ('logcounts') will be from the \code{data} slot of the \code{Seurat::DefaultAssay}.
#' @param colour_by The slot to colour the top right plot by. The default is 'cluster'. Another reasonable value is to add cell labels
#'  to \code{scseq} that were determined from two seperate samples and then see how Seurat data integration warps those groups.
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

explore_scseq_clusters <- function(scseq, markers = NULL, assay.type = 'logcounts', colour_by = 'cluster') {

  # setup ----
  if (is.null(markers))
    markers <- get_scseq_markers(scseq, assay.type)

  # plots all based on SingleCellExperiment
  # use non-integrated but corrected assay for marker gene plots
  sce <- srt_to_sce(scseq, "SCT")

  current_markers <- c()
  subclusters <- list()
  subcluster_markers <- list()
  lack_subclusters <- c()
  prev_gene <- NULL
  prev_cluster <- NULL

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
                      shiny::uiOutput('subcluster_ui'),
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
    # the current subset of the original SingleCellExperiment
    # if exploring subset, will look for new clusters/markers
    sce_r <- shiny::reactive({
      sce_orig <- sce

      if (in_subcluster_r() & !is.null(cluster_r())) {
        # show selected subcluster only
        sce <- sce[, sce$cluster == cluster_r()]

        if (isFALSE(cluster_r() %in% names(subcluster_markers))) {
          # subset original data to selected cluster
          scseq <- scseq[, scseq$seurat_clusters == cluster_r()]

          groups <- unique(scseq$orig.ident)
          if (length(groups) > 1) {
            # if there are user-labeled groups, compare them
            Seurat::Idents(scseq) <- scseq$orig.ident

          } else {
            # otherwise calculate and compare subclusters
            scseq <- tryCatch(add_scseq_clusters(scseq), error = function(e) return(scseq))
            subcl <- letters[scseq$seurat_clusters]

            # error or one subcluster
            if(length(unique(subcl)) == 1)
              lack_subclusters <<- c(lack_subclusters, cluster_r())

            names(subcl) <- colnames(scseq)
            # not sure if all of these are necessary
            scseq$seurat_clusters <- scseq$orig.ident <- Seurat::Idents(scseq) <- factor(subcl)
          }
          subclusters[[cluster_r()]] <<- scseq$orig.ident
          subcluster_markers[[cluster_r()]] <<- get_scseq_markers(scseq)
        }

        if (cluster_r() %in% lack_subclusters) {
          sce <- sce_orig
          sce@metadata$colour_by <- colour_by
          return(sce)
        }

        sce$orig.ident <- subclusters[[cluster_r()]]
        sce@metadata$colour_by <- 'orig.ident'

      } else {
        sce@metadata$colour_by <- colour_by
      }
      return(sce)
    })

    # selected gene: use previously selected if delete selection
    gene_r <- shiny::reactive({
      gene <- input$gene
      if (is.null(gene) || gene == '') {
        gene <- prev_gene
      } else {
        prev_gene <<- gene
      }
      return(gene)
    })

    # used to determine available groups to show (e.g. ctrl and test)
    available_groups_r <- shiny::reactive({
      sce <- sce_r()
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

    # boolean indicating if currently exploring subcluster
    in_subcluster_r <- shiny::reactive({
      if (!length(input$subcluster)) return(FALSE)
      return(input$subcluster)
    })

    # used to get the high level cluster (not the subcluster)
    cluster_r <- shiny::reactive({
      cluster <- input$cluster
      if (is.null(cluster) ||
          # after switching back to all clusters subcluster is FALSE but input$cluster is a subcluster
          # this test prevents saving the subcluster as prev_cluster and prevents returning to Cluster 0
          (!in_subcluster_r() && !cluster %in% sce$orig.ident && !cluster %in% letters))
        prev_cluster <<- cluster

      return(prev_cluster)
    })

    # update marker genes based on cluster/subcluster selection -----
    shiny::observeEvent(input$cluster, {
      cluster_markers <- current_markers[[input$cluster]]
      choices <- row.names(cluster_markers)
      shinyWidgets::updatePickerInput(session, 'gene', choices = choices)
    })


    # dynamic UI inputs -------
    # toggle for e.g. showing ctrl and/or test cells
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

    # UI elements change based on if exploring all or subcluster
    output$subcluster_ui <- shiny::renderUI({

      if (in_subcluster_r() && isFALSE(cluster_r() %in% lack_subclusters)) {
        # in subcluster and exists

        # show subclusters markers
        current_markers <<- subcluster_markers[[cluster_r()]]
        cluster_choices <- names(current_markers)

        # sometimes a race condition happens where updating UI before have subcluster markers
        # this prevents failure
        if (is.null(cluster_choices)) cluster_choices <- ''

        # use either selected subcluster or first if none
        selected <- ifelse(isTRUE(input$cluster %in% cluster_choices),
                           input$cluster, cluster_choices[1])

        # style options to indicate top level cluster
        choicesOpt <- list(content = paste0('<span><span style="color: darkgray;">',
                                            'Cluster ', cluster_r(), ' </span>',
                                            shiny::icon('chevron-right', 'sub-chev'),
                                            cluster_choices, '</span>'))

        # cues to return to all clusters
        title <- 'show all clusters'
        icon <- shiny::icon('chevron-left', 'fa-fw')
        value <- TRUE
        disabled <- FALSE

      } else {
        # either not in subcluster or don't exist

        # use top level cluster markers and previous cluster
        current_markers <<- markers
        cluster_choices <- names(current_markers)
        names(cluster_choices) <- paste('Cluster', cluster_choices)
        selected <- ifelse(is.null(cluster_r()),
                           cluster_choices[1], cluster_r())

        # no styling for options
        choicesOpt <- NULL

        # cues to enter subcluster
        title <- 'explore this cluster'
        icon <- shiny::icon('chevron-right', 'fa-fw')
        value <- FALSE
        disabled <- FALSE

        if (isTRUE(cluster_r() %in% lack_subclusters)) {
          # don't have subcluster

          # cues to prevent entering subclusters
          title <- 'no subclusters'
          icon <- shiny::icon('ban', 'fa-fw')
          disabled <- TRUE
        }
      }

      shiny::tags$div(
        shiny::tags$div(style = "display:inline-block; text-overflow:",
                        shinyWidgets::pickerInput("cluster", 'Show marker genes for:', width = '250px',
                                                  choices = cluster_choices,
                                                  choicesOpt = choicesOpt,
                                                  selected = selected)),

        shinyBS::bsButton('subcluster', label = '', type = 'toggle',
                          title = title, icon = icon, value = value, disabled = disabled)
      )
    })

    # plots ------
    # show tSNE plot coloured by expression values
    output$marker_plot <- shiny::renderPlot({
      sce <- sce_r()
      gene <- gene_r()
      if (is.null(gene)) return(NULL)

      # make selected groups stand out
      point_alpha <- rep(1, ncol(sce))
      point_alpha[!sce$orig.ident %in% selected_groups_r()] <- 0


      suppressMessages(scater::plotTSNE(sce,
                                        by_exprs_values = assay.type,
                                        colour_by = gene,
                                        point_size = 3,
                                        point_alpha = point_alpha,
                                        theme_size = 14) +
                         ggplot2::scale_fill_distiller(palette = 'Reds', name = gene, direction = 1))

    })


    # show plot of predicted cell clusters
    output$cluster_plot <- shiny::renderPlot({
      sce <- sce_r()

      # make selected cluster and groups stand out
      point_alpha <- rep(1, ncol(sce))
      point_alpha[sce$cluster == cluster_r()] <- 1
      point_alpha[!sce$orig.ident %in% selected_groups_r()] <- 0

      legend_title <- ifelse(in_subcluster_r(), 'Group', 'Cluster')

      scater::plotTSNE(sce,
                       by_exprs_values = assay.type,
                       colour_by = sce@metadata$colour_by,
                       point_size = 3,
                       point_alpha = point_alpha,
                       theme_size = 14) +
        ggplot2::guides(fill = ggplot2::guide_legend(legend_title)) +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank())

    })


    # plot BioGPS data
    output$biogps <- shiny::renderPlot({
      plot_biogps(gene_r())
    })


    # link to GeneCards page for gene ----
    gene_link_r <- shiny::reactive({
      return(paste0('https://www.genecards.org/cgi-bin/carddisp.pl?gene=', gene_r()))
    })

    shiny::observeEvent(input$genecards, {
      utils::browseURL(gene_link_r())
    })

  }
  shiny::runGadget(shiny::shinyApp(ui, server), viewer = shiny::browserViewer())
}



#' Plot BioGPS data for a HGNC symbol
#'
#' @param gene Character vector of gene name.
#' @keywords internal
#'
#' @return ggplot
#' @export
#'
#' @examples
plot_biogps <- function(gene) {
  if (!length(gene) || !gene %in% biogps[, SYMBOL]) return(NULL)

  gene_dat <- unlist(biogps[gene, -c('ENTREZID', 'SYMBOL')])
  gene_dat <- sort(gene_dat, decreasing = TRUE)[1:20]
  gene_dat <- tibble::tibble(mean = gene_dat,
                             source = factor(names(gene_dat), levels = rev(names(gene_dat))))

  ggplot2::ggplot(gene_dat, ggplot2::aes(x = source, y = mean, fill = mean)) +
    ggplot2::theme_minimal() +
    ggplot2::geom_bar(stat = "identity", color = 'black', size = 0.1, width = 0.7) +
    ggplot2::scale_fill_distiller(palette = 'Reds', name = '', direction = 1) +
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
}
