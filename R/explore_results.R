
#' Explore query results
#'
#' @param query_res \code{data.frame} returned from \code{\link{append_pdata}}
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' # load CMAP02 data
#' cmap_path <- system.file('extdata', 'cmap_es_ind.rds', package = 'drugseqr')
#' cmap_es <- readRDS(cmap_path)
#'
#' # load previous differential expression analysis
#' data_dir <- 'data-raw/example-data'
#' anal <- readRDS(file.path(data_dir, 'diff_expr_symbol.rds'))
#'
#' # get dprime effect size values for analysis
#' dprimes <- get_dprimes(anal)
#'
#' # get correlations between query and drug signatures
#' res <- query_drugs(dprimes, cmap_es)
#'
#' # append drug pdata
#' res <- append_pdata(res, 'CMAP02')
#'
#' explore_results(res)
#'
explore_results <- function(query_res) {

  # setup ----

  # order by increasing correlation
  query_res <- query_res[order(query_res$Correlation), ]

  # add linkout to pubchem
  cids <- query_res$`Pubchem CID`
  query_res$`Pubchem CID` <- paste0('<a href="https://pubchem.ncbi.nlm.nih.gov/compound/',  cids, '" target="_blank">', cids, '</a>')

  # replace correlation with svg element
  cors <- query_res$Correlation
  query_res$Correlation <- paste0('<svg class="simplot" width="180" height="38">
                            <line x1="90" x2="90" y1="0" y2="38" style="stroke: rgb(221, 221, 221); shape-rendering: crispEdges; stroke-width: 1px; stroke-dasharray: 3, 3;"></line>
                            <g><text x="', calcx(cors, range(cors)), '" y="38" class="x text" dy="-2">', signif(cors, 3), '</text></g>
                            <g><circle cx="', calcx(cors, range(cors)), '" cy="19" r="5" class="cor"></circle></g>
                            </svg>')



  #  user interface ----

  ui <- miniUI::miniPage(
    shinyjs::useShinyjs(),
    shiny::tags$head(
      shiny::tags$style("#query_res {white-space: nowrap;}"), # table text on 1 line
      shiny::tags$style("td.sim-cell {min-width: 180px !important; padding: 0px !important; !important; position: relative;}"),
      shiny::tags$style(".simplot {position: absolute; z-index: 1, top: 0; left: 0; right: 0; bottom: 0;}"),
      shiny::tags$style("text.x {fill: #ddd; font: 10px Arial, sans-serif; text-anchor: middle;}"),
      shiny::tags$style(".simplot circle {fill: transparent; stroke-width: 1.1px; stroke: rgba(0, 0, 0, 0.75);}"),
      shiny::tags$style(".simplot:hover text.x {fill: #443;}"),
      shiny::tags$style(".simplot:hover circle {fill: red; stroke: #443;}")
    ),
    # title bar
    miniUI::gadgetTitleBar("Explore Results"),
    miniUI::miniContentPanel(
      shiny::fillCol(flex = c(NA, NA, 1),
                     shiny::selectizeInput('study',
                                           'Select study:',
                                           choices = c('CMAP02', 'L1000')),
                     shiny::hr(),
                     DT::dataTableOutput("query_res")
      )
    )
  )

  # server ----

  server <- function(input, output, session) {

    # show query data
    output$query_res <- DT::renderDataTable({

      DT::datatable(
        query_res,
        class = 'cell-border',
        rownames = FALSE,
        selection = 'none',
        escape = FALSE, # to allow HTML in table
        options = list(
          ordering=FALSE,
          columnDefs = list(list(className = 'dt-nopad sim-cell', height=38, width=120, targets = 0)),
          scrollY = TRUE,
          pageLength = 20,
          paging = TRUE,
          bInfo = 0,
          dom = 'ftp'
        )
      )
    }, server = TRUE)


    # click 'Done' ----

    shiny::observeEvent(input$done, {
      shiny::stopApp()
    })

  }

  shiny::runGadget(shiny::shinyApp(ui, server), viewer = shiny::browserViewer())
}

#' Calculate x position for correlation plot
#'
#' @param cor Numeric vector of correlation values.
#' @param range Numeric vector of length 2 specifying maximum and minimum values of \code{cor}.
#' @param width Plot width to scale correlation values to.
#' @param pad Numeric value that is respectively, subtracted and added to values in \code{range}. Make it so that circles and
#' correlation text values don't get cut off.
#'
#' @return Numeric vector giving x position for correlation plot in \code{\link{explore_search}}
#' @export
#'
#' @examples
calcx <- function(cor, range = c(-1, 1), width = 180, pad = 0.1) {
  range[1] <- range[1] - pad
  range[2] <- range[2] + pad
  (cor - range[1])/diff(range) * width
}
