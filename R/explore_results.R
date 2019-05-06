
#' Explore query results
#'
#' @param cmap_res Named numeric vector returned from \code{\link{query_drugs}} using CMAP02 data.
#' @param l1000_res Named numeric vector returned from \code{\link{query_drugs}} using L1000 data.
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
#' cmap_res <- query_drugs(dprimes, cmap_es)
#'
#' explore_results(cmap_res)
#'
explore_results <- function(cmap_res = NULL, l1000_res = NULL) {


  # setup ----

  null_cmap <- is.null(cmap_res)
  null_l1000 <- is.null(l1000_res)

  if (null_cmap & null_l1000)
    stop('Must provide one of cmap_res or l1000_res')

  # append pdata to results and add html for correlation plots
  study_choices <- c('CMAP02', 'L1000')[!is.null(c(cmap_res, l1000_res))]
  study_tables <- list()

  if (!null_cmap)
    study_tables[['CMAP02']] <- study_table(cmap_res, 'CMAP02')

  # if cmap, don't setup l1000 until selection (for initial load speed)
  if (!null_l1000 & null_cmap)
    study_tables[['L1000']] <- study_table(l1000_res, 'L1000')

  # initial query_res
  query_res <- study_tables[[1]]

  #  user interface ----

  ui <- miniUI::miniPage(
    shinyjs::useShinyjs(),
    shiny::tags$head(
      shiny::tags$style("#clinical {margin-top: -26px;}"), # to align clinical button and study selection
      shiny::tags$style("#query_res {white-space: nowrap;}"), # table text on 1 line

      # css for correlation plots
      shiny::tags$style("td.sim-cell {min-width: 180px !important; padding: 0px !important; !important; position: relative;}"),
      shiny::tags$style(".simplot {position: absolute; z-index: 1, top: 0; left: 0; right: 0; bottom: 0;}"),
      shiny::tags$style(".simplot text.x {fill: transparent; font: 10px Arial, sans-serif; text-anchor: middle;}"),
      shiny::tags$style(".simplot circle {fill: transparent; stroke-width: 1.1px; stroke: rgba(0, 0, 0, 0.35);}"),
      shiny::tags$style(".cor-point:hover text.x {fill: #443;}"),
      shiny::tags$style(".cor-point:hover circle {fill: red; stroke: #443;}"),

      shiny::includeScript(system.file("js/toggleClinicalTitle.js", package = "drugseqr"))
    ),
    # title bar
    miniUI::gadgetTitleBar("Explore Results"),
    miniUI::miniContentPanel(
      shiny::fillCol(flex = c(NA, NA, 1),
                     shiny::tags$div(
                       shiny::tags$div(style = "display:inline-block", shiny::selectizeInput('study', 'Select study:', choices = study_choices)),
                       shinyBS::bsButton('clinical', label = '', icon = shiny::icon('pills'), style='default', onclick='toggleClinicalTitle(this)', title = 'only show compounds with a clinical phase')
                     ),
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
        query_res(),
        class = 'cell-border',
        rownames = FALSE,
        selection = 'none',
        escape = FALSE, # to allow HTML in table
        options = list(
          columnDefs = list(list(className = 'dt-nopad sim-cell', height=38, width=120, targets = 0),
                            list(targets = 8:14, render = DT::JS(
                              "function(data, type, row, meta) {",
                              "return type === 'display' && data !== null && data.length > 12 ?",
                              "'<span title=\"' + data + '\">' + data.substr(0, 12) + '...</span>' : data;",
                              "}"))),
          scrollY = TRUE,
          scrollX = TRUE,
          pageLength = 50,
          paging = TRUE,
          bInfo = 0,
          dom = 'ftp'
        )
      )
    }, server = TRUE)

    isClinical <- shiny::reactiveVal(FALSE)

    # query_res reactive to study choice
    query_res <- shiny::reactive({
      # initial setup for l1000 only after selection
      if (input$study == 'L1000' && !'L1000' %in% names(study_tables))
        study_tables[['L1000']] <<- study_table(l1000_res, 'L1000')

      query_res <- study_tables[[input$study]]

      # for removing entries without a clinical phase
      if (isClinical()) {
        query_res <- tibble::as_tibble(query_res)
        query_res <- dplyr::filter(query_res, !is.na(.data$`Clinical Phase`))
      }

      return(query_res)

    })

    # click 'Clinical' ----
    observeEvent(input$clinical,{
      toggle <- (input$clinical %% 2) + 1

      # update clinical button styling
      shinyBS::updateButton(session, 'clinical', style = c('default', 'primary')[toggle])

      # update boolean reactiveVal
      isClinical(toggle - 1)
    })

    # click 'Done' ----

    shiny::observeEvent(input$done, {
      shiny::stopApp()
    })

  }

  shiny::runGadget(shiny::shinyApp(ui, server), viewer = shiny::browserViewer())
}

#' Title
#'
#' @param query_res
#' @param study
#'
#' @return
#' @export
#'
#' @examples
study_table <- function(query_res, study) {
  query_res <- append_pdata(query_res, study)
  query_res <- summarize_compound(query_res)
  query_res <- add_table_html(query_res)
  return(query_res)
}


#' Title
#'
#' @param query_res
#'
#' @return
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
summarize_compound <- function(query_res) {

  # group by compound
  query_res <- query_res %>%
    dplyr::group_by(Compound)

  # put all correlations together in list
  # keep minimum correlation for sorting
  query_cors <- query_res %>%
    dplyr::select(Correlation, Compound) %>%
    dplyr::summarise(min_cor = min(Correlation), Correlation = I(list(Correlation))) %>%
    dplyr::select(Correlation, min_cor)

  # keep furthest clinical phase
  query_phase <- query_res %>%
    dplyr::select(clinical_phase, Compound) %>%
    dplyr::summarise(clinical_phase = max(clinical_phase)) %>%
    dplyr::pull(clinical_phase)


  # summarize cell lines etc
  query_rest <- query_res %>%
    dplyr::select(-Correlation, -clinical_phase) %>%
    dplyr::summarize_all(function(x) {
      unqx <- na.omit(x)
      unqx <- as.character(unqx)
      unqx <- unlist(strsplit(unqx, '|', fixed = TRUE))
      unqx <- unique(unqx)

      # keep as NA if they all are
      if (!length(unqx)) return(NA_character_)

      # collapse distinct non-NA entries
      return(paste(unqx, collapse = ' | '))
    }) %>%
    tibble::add_column('Clinical Phase' = query_phase, .after = 'Compound')

  query_res <- dplyr::bind_cols(query_cors, query_rest) %>%
    dplyr::arrange(min_cor) %>%
    dplyr::select(-min_cor)

  return(query_res)
}

#' Add HTML to query results table
#'
#' @param query_res \code{data.frame} returned by \code{\link{append_pdata}}
#'
#' @return \code{query_res} with pubchem cid links and correlation plot html.
#' @export
#'
#' @examples
add_table_html <- function(query_res) {

  # add linkout to pubchem
  cids <- query_res$`Pubchem CID`
  have_cid <- !is.na(cids)
  query_res$`Pubchem CID`[have_cid] <- paste0('<a href="https://pubchem.ncbi.nlm.nih.gov/compound/',
                                              cids[have_cid],
                                              '" target="_blank">',
                                              cids[have_cid], '</a>')

  # replace correlation with svg element
  cors <- query_res$Correlation
  cors_range <- range(unlist(cors))
  query_res$Correlation <- paste0('<svg class="simplot" width="180" height="38">
                            <line x1="90" x2="90" y1="0" y2="38" style="stroke: rgb(221, 221, 221); shape-rendering: crispEdges; stroke-width: 1px; stroke-dasharray: 3, 3;"></line>',
                            add_cors_html(cors, cors_range),
                            '</svg>')

  return(query_res)
}

add_cors_html <- function(cors, cors_range) {

  cors_html <- sapply(cors, function(x) {
    paste0('<g class="cor-point">
              <g><text x="', calcx(x, cors_range), '" y="38" class="x text" dy="-2">', signif(x, 3), '</text></g>
              <g><circle cx="', calcx(x, cors_range), '" cy="19" r="5" class="cor"></circle></g>
            </g>', collapse = '\n')
  })

  return(cors_html)
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
