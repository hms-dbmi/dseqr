
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
#' l1000_path <- system.file('extdata', 'l1000_es.rds', package = 'drugseqr')
#' l1000_es <- readRDS(l1000_path)
#'
#' # run differential expression analysis
#' data_dir <- 'data-raw/example-data'
#' eset <- readRDS(system.file('extdata', 'IBD', 'eset.rds', package = 'drugseqr'))
#' anal <- diff_expr(eset, data_dir)
#'
#' # alternatively load previous analysis
#' anal <- readRDS(file.path(data_dir, 'diff_expr_symbol.rds'))
#'
#' # get dprime effect size values for analysis
#' dprimes <- get_dprimes(anal)
#'
#' # get correlations between query and drug signatures
#' cmap_res <- query_drugs(dprimes, cmap_es)
#' l1000_res <- query_drugs(dprimes, l1000_es)
#'
#' explore_results(cmap_res, l1000_res)
#'
explore_results <- function(cmap_res = NULL, l1000_res = NULL) {


  # setup ----

  null_cmap <- is.null(cmap_res)
  null_l1000 <- is.null(l1000_res)

  if (null_cmap & null_l1000)
    stop('Must provide one of cmap_res or l1000_res')

  # append pdata to results and add html for correlation plots
  study_choices <- c('CMAP02', 'L1000')[!c(is.null(cmap_res), is.null(l1000_res))]
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
      shiny::tags$style(".cor-point:hover {cursor: crosshair;}"),
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
                       shinyBS::bsButton('clinical', label = '', icon = shiny::icon('pills'), style='default', onclick='toggleClinicalTitle(this)', title = 'only show compounds with a clinical phase'),
                       shiny::plotOutput(outputId = "histPlot", width = 800)
                     ),
                     shiny::hr(),
                     DT::dataTableOutput("query_res")
      )
    )
  )

  # server ----

  server <- function(input, output, session) {

    output$histPlot <- shiny::renderPlot({

      if (input$study == 'CMAP02') x <- cmap_res else x <- l1000_res
      bins <- seq(min(x), max(x), length.out = 51)

      hist(x, breaks = bins, col = "#75AADB", border = "white",
           xlab = "Pearson correlations", main = "", xlim = c(min(x) - 0.1, max(x) + 0.1))
    })

    # show query data
    output$query_res <- DT::renderDataTable({

      query_res <- query_res()
      wide_cols <- c('MOA', 'Target', 'Disease Area', 'Indication', 'Vendor', 'Catalog #', 'Vendor Name')
      # -1 needed with rownames = FALSE
      elipsis_targets <- which(colnames(query_res) %in% wide_cols) - 1

      DT::datatable(
        query_res,
        class = 'cell-border',
        rownames = FALSE,
        selection = 'none',
        escape = FALSE, # to allow HTML in table
        options = list(
          columnDefs = list(list(className = 'dt-nopad sim-cell', height=38, width=120, targets = 0),
                            list(targets = elipsis_targets, render = DT::JS(
                              "function(data, type, row, meta) {",
                              "return type === 'display' && data !== null && data.length > 17 ?",
                              "'<span title=\"' + data + '\">' + data.substr(0, 17) + '...</span>' : data;",
                              "}"))),
          ordering=FALSE,
          scrollX = TRUE,
          pageLength = 50,
          scrollY = TRUE,
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

#' Set up DT for explore_results
#'
#' Appends perturbation annotations, summarizes by perturbation, and adds HTML for
#' correlation plots and hyperlinks.
#'
#' @param query_res Named numeric vector returned from \code{\link{query_drugs}}.
#' @param study Character vector. Either \code{'CMAP02'} or \code{'L1000'}.
#'
#' @return \code{data.frame} of perturbation correlations and annotations.
#' @export
#'
#' @examples
study_table <- function(query_res, study) {
  query_res <- append_pdata(query_res, study)
  query_res <- summarize_compound(query_res)
  query_res <- add_table_html(query_res)
  return(query_res)
}


#' Summarize query results and annotations by perturbation
#'
#' Takes a \code{data.frame} with one row per signatures and summarizes to one row per compound.
#'
#' Variables related to individual signatures (cell line, dose, duration, and sample number) are
#' pasted together and added as a list to \code{'title'} column. Query correlation values are also added as a list to
#' the \code{'Correlation'} column.
#'
#' Clinical status is summarized by keeping the most advanced phase only (e.g. Launched > Phase 3). For all other variables,
#' all unique entries are paste together seperated by \code{'|'}.
#'
#' @param query_res \code{data.frame} of perturbation correlations and annotations returned by \code{\link{append_pdata}}.
#'
#' @return \code{data.frame} of perturbation correlations and annotations summarized by perturbation.
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#'
#'
summarize_compound <- function(query_res) {

  # group by compound
  query_res <- query_res %>%
    dplyr::group_by(Compound)

  # merge cell line, dose, and duration and samples for correlation titles
  query_res <- query_res %>%
    tidyr::unite(title, 'Cell Line', 'Dose', 'Duration', 'Samples(n)')

  query_title <- query_res %>%
    dplyr::select(title, Compound) %>%
    dplyr::summarize(title = I(list(title))) %>%
    dplyr::select(title)

  # put all correlations together in list
  # keep minimum correlation for sorting
  query_cors <- query_res %>%
    dplyr::select(Correlation, Compound) %>%
    dplyr::summarise(min_cor = min(Correlation), Correlation = I(list(Correlation))) %>%
    dplyr::select(Correlation, min_cor)

  # keep furthest clinical phase
  query_phase <- query_res %>%
    dplyr::select(`Clinical Phase`, Compound) %>%
    dplyr::summarise(`Clinical Phase` = max(`Clinical Phase`)) %>%
    dplyr::pull(`Clinical Phase`)

  # summarize rest
  query_rest <- query_res %>%
    dplyr::select(-Correlation, -`Clinical Phase`, -title) %>%
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

  query_res <- dplyr::bind_cols(query_cors, query_rest, query_title) %>%
    dplyr::arrange(min_cor) %>%
    dplyr::select(-min_cor)

  return(query_res)
}

#' Add linkout HTML
#'
#' Non \code{NA} values in \code{id_col} of \code{query_res} are inserted between
#' \code{pre_url} and \code{post_url} to form hyperlinks. Relevant HTML markup is also added.
#'
#'
#' @param query_res \code{data.frame} returned by \code{\link{summarize_compound}}.
#' @param id_col Character with column in \code{query_res} that contains ids to
#'   be inserted between \code{pre_url} and \code{post_url} to form the link. \code{NA}
#'   values will be ignored.
#' @param img_url Character with url to an image to display as hyperlink.
#' @param pre_url Character with url portion to paste before \code{id_col} column values.
#' @param post_url Character with url portion to paste before \code{id_col} column values.
#' @param title Character that will be added to hyperlink title attribute. Default is \code{id_col}.
#'
#' @return \code{query_res} with HTML for hyperlinks in \code{id_col}.
#' @export
#'
#' @examples
add_linkout <- function(query_res, id_col, img_url, pre_url, post_url = NULL, title = id_col) {

  ids <- query_res[[id_col]]
  have_ids <- !is.na(ids)
  query_res[[id_col]][have_ids] <- paste0('<a href="',
                                          pre_url, ids[have_ids], post_url,
                                          '" target="_blank" title="', paste('Go to', title), '">',
                                          '<img src="', img_url, '" height="22px" hspace="4px"></img>',
                                          # ids[have_ids],
                                          '</a>')

  return(query_res)
}

#' Add HTML to query results table
#'
#' @param query_res \code{data.frame} returned by \code{\link{summarize_compound}}
#'
#' @return \code{query_res} with pubchem cid links and correlation plot HTML.
#' @export
#'
#' @examples
add_table_html <- function(query_res) {

  pre_urls <- c('https://pubchem.ncbi.nlm.nih.gov/compound/',
                'http://sideeffects.embl.de/drugs/',
                'https://www.drugbank.ca/drugs/')

  img_urls <- c('https://pubchem.ncbi.nlm.nih.gov/pcfe/favicon/favicon.ico',
                'http://sideeffects.embl.de/media/images/EMBL_Logo.png',
                'https://www.drugbank.ca/favicons/favicon.ico')


  # add linkout to Pubchem, SIDER, and DrugBank
  query_res <- add_linkout(query_res, 'Pubchem CID', img_urls[1], pre_urls[1], title = 'Pubchem')
  query_res <- add_linkout(query_res, 'SIDER', img_urls[2], pre_urls[2])
  query_res <- add_linkout(query_res, 'DrugBank', img_urls[3], pre_urls[3])

  # merge linkouts into single column
  query_res <- merge_linkouts(query_res, c('Pubchem CID', 'DrugBank', 'SIDER'))

  # replace correlation with svg element
  cors <- query_res$Correlation

  # move titles to plots
  titles <- query_res$title
  query_res$title <- NULL

  cors_range <- range(unlist(cors))
  query_res$Correlation <- paste0('<svg class="simplot" width="180" height="38">
                            <line x1="90" x2="90" y1="0" y2="38" style="stroke: rgb(221, 221, 221); shape-rendering: crispEdges; stroke-width: 1px; stroke-dasharray: 3, 3;"></line>',
                            get_cors_html(cors, titles, cors_range),
                            '</svg>')

  return(query_res)
}

#' Merge columns with image links
#'
#' @param query_res \code{data.frame} after calling \code{\link{add_linkout}} to \code{cols}.
#' @param cols Character vector of columns in \code{query_res} that \code{\link{add_linkout}} has been called on.
#
#' @importFrom magrittr "%>%"
#'
#' @return \code{query_res} with column \code{'External Links'} formed from pasting \code{cols} together. \code{cols} are removed.
#' @export
#'
#' @examples
merge_linkouts <- function(query_res, cols) {

  # paste cols with non-NA values
  paste.na <- function(x) paste(x[!is.na(x)], collapse = '')
  new_vals <- apply(query_res[ ,cols] , 1, paste.na)

  # remove cols that pasted
  query_res <- query_res %>%
    tibble::add_column(`External Links` = new_vals, .before = cols[1]) %>%
    dplyr::select(-cols)

  return(query_res)
}

#' Get HTML for correlation values.
#'
#' @param cors List of numeric vectors of pearson correlations.
#' @param titles List of character vectors of treatment titles for pearson correlations (e.g. MCF7_1e-05M_6h_3).
#' @param cors_range Numeric vector of length two specifying the range of correlation values.
#'
#' @return Character vector of HTML markup for the title/circle/text for a correlation plot.
#' @export
#'
#' @examples
get_cors_html <- function(cors, titles, cors_range) {

  cors_html <- sapply(seq_along(cors), function(i) {
    x <- cors[[i]]
    xtitle <- titles[[i]]
    paste0('<g class="cor-point">
              <title>', xtitle, '</title>
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
