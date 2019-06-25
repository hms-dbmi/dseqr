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

explore_scseq_clusters <- function(scseq, markers = NULL, pt.size = 3) {

  # setup ----
  if (is.null(markers))
    markers <- get_scseq_markers(scseq)

  # name cluster choices for drop down
  cluster_choices <- names(markers)
  prefix <- ifelse(is.na(as.integer(cluster_choices[1])), '', 'Cluster ')
  names(cluster_choices) <- paste0(prefix, names(markers))

  test_cluster <- NULL

  #  user interface ----

  ui <- miniUI::miniPage(
    shiny::tags$head(
      shiny::tags$style("#genecards, #show_contrasts {border-top-left-radius: 0; border-bottom-left-radius: 0; position: relative; z-index: 2; margin-left: -8px; margin-top: -26px;}")
    ),
    miniUI::gadgetTitleBar("Explore Single-Cell Clusters"),
    miniUI::miniContentPanel(
      shiny::fluidRow(
        shiny::column(6,
                      shiny::tags$div(
                        shiny::tags$div(style = "display:inline-block;",
                                        shiny::selectizeInput("selected_cluster", 'Show marker genes for:', choices = cluster_choices, width = '250px')),
                        shinyBS::bsButton('show_contrasts', label = '', type = 'toggle', icon = shiny::icon('chevron-right', 'fa-fw'), style='default', title = 'Toggle single group comparisons')
                      ),
                      shiny::br(),
                      shiny::tags$div(
                        shiny::tags$div(style = "display:inline-block;",
                                        shiny::selectizeInput('gene', 'Show expression for:', choices = NULL, width = '250px')),
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
      if (length(groups) == 1 || input$groups == 'all') return(groups)

      return(input$groups)
    })

    # update cluster choices based on show_contrasts toggle
    cluster_choices_r <- shiny::eventReactive(input$show_contrasts, {

      if (input$show_contrasts) {
        # cluster choices are as compared to other clusters
        test <- input$selected_cluster
        ctrls <- cluster_choices[cluster_choices != test]

        contrast_choices <- c(test, paste0(test, ' - ', ctrls))
        names(contrast_choices) <- paste0(prefix, test, ' vs ', c('all', ctrls))

        # update global so that can return to same cluster when toggle back
        test_cluster <<- test

      } else {
        # cluster choices are the clusters themselves
        contrast_choices  <- cluster_choices
      }

      # update icon on toggle
      icon <- ifelse(input$show_contrasts, 'chevron-down', 'chevron-right')
      shinyBS::updateButton(session, 'show_contrasts', icon = shiny::icon(icon, 'fa-fw'))

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
      # get cluster if don't have (for comparing specific cluster)
      if (!sel %in% names(markers)) {
        con <- strsplit(sel, ' - ')[[1]]
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
                                          choices = c('all', groups), inline = TRUE),
          shiny::br()
        )
      }
      return(groups_toggle)
    })


    # plots ------
    # show tSNE plot coloured by expression values
    output$marker_plot <- shiny::renderPlot({
      if (input$gene == '') return(NULL)
      plot_umap_gene(scseq, input$gene, selected_idents = selected_groups_r(), pt.size = pt.size)
    })


    # show plot of predicted cell clusters
    output$cluster_plot <- shiny::renderPlot({
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

  }
  shiny::runGadget(shiny::shinyApp(ui, server), viewer = shiny::browserViewer())
}
