#' Logic for Drugs page
#' @export
#' @keywords internal
drugsPage <- function(input, output, session, new_anal, data_dir) {

  # the form area inputs/results
  form <- callModule(drugsForm, 'form',
                     new_anal = new_anal,
                     data_dir = data_dir)

  # the output table
  callModule(drugsTable, 'table',
             query_res = form$query_res,
             drug_study = form$drug_study,
             cells = form$cells,
             sort_by = form$sort_by,
             show_clinical = form$show_clinical)




}

# Logic for form on drugs page
#' @export
#' @keywords internal
drugsForm <- function(input, output, session, new_anal, data_dir) {

  querySignature <- callModule(querySignature, 'signature',
                               new_anal = new_anal,
                               data_dir = data_dir)

  drugStudy <- callModule(selectedDrugStudy, 'drug_study',
                          anal = querySignature$anal)

  advancedOptions <- callModule(advancedOptions, 'advanced',
                                cmap_res = querySignature$cmap_res,
                                l1000_res = querySignature$l1000_res,
                                drug_study = drugStudy$drug_study,
                                show_advanced = drugStudy$show_advanced)

  query_res <- reactive({
    drug_study <- drugStudy$drug_study()
    cmap_res <- querySignature$cmap_res()
    l1000_res <- querySignature$l1000_res()

    if (drug_study == 'CMAP02') return(cmap_res)
    if (drug_study == 'L1000') return(l1000_res)

    return(NULL)
  })




  return(list(
    query_res = query_res,
    drug_study = drugStudy$drug_study,
    cells = advancedOptions$cells,
    sort_by = advancedOptions$sort_by,
    show_clinical = drugStudy$show_clinical
  ))


}

#' Logic for query signature in drugsForm
#' @export
#' @keywords internal
querySignature <- function(input, output, session, new_anal, data_dir) {


  cmap_res <- reactiveVal()
  l1000_res <- reactiveVal()

  # reload query choices if new analysis
  anals <- reactive({
    new_anal()
    load_bulk_anals(data_dir)
  })

  observe({
    anals <- anals()
    req(anals)
    updateSelectizeInput(session, 'query', choices = rbind(rep(NA, 5), anals), server = TRUE)
  })

  anal <- reactive({
    row_num <- input$query
    anals <- anals()
    req(row_num, anals)

    anals[row_num, ]
  })


  # paths to analysis and drug query results
  res_paths <- reactive({
    anal <- anal()
    dataset_dir <- file.path(data_dir, 'bulk', anal$dataset_dir)
    anal_name <- anal$anal_name

    list(
      anal = file.path(dataset_dir, paste0('diff_expr_symbol_', anal_name, '.rds')),
      cmap = file.path(dataset_dir, paste0('cmap_res_', anal_name, '.rds')),
      l1000 = file.path(dataset_dir, paste0('l1000_res_', anal_name, '.rds'))
    )
  })


  # get saved cmap/l1000 query results
  observe({
    res_paths <- res_paths()

    # load if available
    if (file.exists(res_paths$cmap)) {
      cmap_res <- readRDS(res_paths$cmap)
      l1000_res <- readRDS(res_paths$l1000)


    } else {
      # otherwise run
      disable('query')

      progress <- Progress$new(session, min = 0, max = 4)
      progress$set(message = "Querying drugs", value = 1)
      on.exit(progress$close())

      cmap_path  <- system.file('extdata', 'cmap_es_ind.rds', package = 'drugseqr', mustWork = TRUE)
      l1000_path <- system.file('extdata', 'l1000_es.rds', package = 'drugseqr', mustWork = TRUE)

      cmap_es  <- readRDS(cmap_path)
      progress$inc(1)
      l1000_es <- readRDS(l1000_path)
      progress$inc(1)


      # get dprime effect size values for analysis
      anal <- readRDS(res_paths$anal)
      dprimes <- get_dprimes(anal)

      # get correlations between query and drug signatures
      cmap_res <- query_drugs(dprimes, cmap_es)
      l1000_res <- query_drugs(dprimes, l1000_es)
      progress$inc(1)

      saveRDS(cmap_res, res_paths$cmap)
      saveRDS(l1000_res, res_paths$l1000)
      enable('query')
    }

    cmap_res(cmap_res)
    l1000_res(l1000_res)
  })



  return(list(
    cmap_res = cmap_res,
    l1000_res = l1000_res,
    anal = anal
  ))

}

#' Logic for selected drug study in drugsForm
#' @export
#' @keywords internal
selectedDrugStudy <- function(input, output, session, anal) {


  drug_study <- reactive(input$study)

  # boolean for advanced options
  show_advanced <- reactive({
    input$advanced %% 2 != 0
  })

  observe({
    req(anal())
    updateSelectizeInput(session, 'study', choices = c('CMAP02', 'L1000'), selected = NULL)
  })

  # toggle for clinical status
  show_clinical <- reactive({
    input$clinical %% 2 != 0
  })
  observe({
    toggleClass('advanced', 'btn-primary', condition = show_advanced())
    toggleClass('clinical', 'btn-primary', condition = show_clinical())
  })



  return(list(
    drug_study = drug_study,
    show_clinical = show_clinical,
    show_advanced = show_advanced
  ))

}

#' Logic for advanced options in drugsForm
#' @export
#' @keywords internal
advancedOptions <- function(input, output, session, cmap_res, l1000_res, drug_study, show_advanced) {

  # available cell lines
  cmap_cells <- cell_info$cmap
  l1000_cells <- cell_info$l1000

  # update choices for cell lines based on selected study
  cell_choices <- shiny::reactive({
    req(drug_study())
    if (drug_study() == 'L1000') return(l1000_cells)
    else if (drug_study() == 'CMAP02') return(cmap_cells)
  })

  #  toggle  showing advanced options
  shiny::observe({
    toggle('advanced-panel', condition = show_advanced(), anim = TRUE)
  })

  # update choices for cell lines
  shiny::observe({
    shiny::updateSelectizeInput(session, 'cells', choices = cell_choices(), selected = NULL, server = TRUE)
  })

  return(list(
    cells = reactive(input$cells),
    sort_by = reactive(input$sort_by)
  ))

}

#' Logic for drug table
#' @export
#' @keywords internal
#' @importFrom magrittr "%>%"
drugsTable <- function(input, output, session, query_res, drug_study, cells, show_clinical, sort_by) {

  # will update with proxy to analent redraw
  dummy_table <- data.frame('Correlation' = NA,
                            'Compound' = NA,
                            'Clinical Phase' = NA,
                            'External Links' = NA,
                            'MOA' = NA,
                            'Target' = NA,
                            'Disease Area' = NA,
                            'Indication' = NA,
                            'Vendor' = NA,
                            'Catalog #' = NA,
                            'Vendor Name' = NA, check.names = FALSE)
  dummy_rendered <- reactiveVal(FALSE)

  # get either cmap or l1000 annotations
  drug_annot <- reactive({
    drug_study <- drug_study()
    req(drug_study)

    if (drug_study == 'CMAP02') return(cmap_annot)
    else if (drug_study == 'L1000') return(l1000_annot)
  })

  # add annotations to query result
  query_table_full <- reactive({
    query_res <- query_res()
    if (is.null(query_res)) return(NULL)

    drug_annot <- drug_annot()
    req(query_res, drug_annot)
    stopifnot(all.equal(drug_annot$title, names(query_res)))

    tibble::add_column(drug_annot, Correlation = query_res, .before=0)
  })

  # subset to selected cells, summarize by compound, and add html
  query_table_summarised <- reactive({
    query_table_full <- query_table_full()
    if (is.null(query_table_full)) return(NULL)

    query_table <- query_table_full %>%
      limit_cells(cells()) %>%
      summarize_compound() %>%
      add_table_html()
  })

  query_table_final <- reactive({
    query_table <- query_table_summarised()
    if (is.null(query_table)) return(NULL)
    sort_by <- sort_by()

    # subset by clinical phase
    if (show_clinical()) query_table <- dplyr::filter(query_table, !is.na(`Clinical Phase`))

    if (sort_by == 'avg_cor') {
      query_table$Correlation <- gsub('simplot', 'simplot show-meanline', query_table$Correlation)
    }

    # sort as desired
    dplyr::arrange(query_table, !!sym(sort_by)) %>%
      select(-min_cor, -avg_cor)
  })


  # show query data
  output$query_table <- DT::renderDataTable({
    # ellipses for wide columns
    wide_cols <- c('MOA', 'Target', 'Disease Area', 'Indication', 'Vendor', 'Catalog #', 'Vendor Name')
    # -1 needed with rownames = FALSE
    elipsis_targets <- which(colnames(dummy_table) %in% wide_cols) - 1
    dummy_rendered(TRUE)

    DT::datatable(
      dummy_table,
      class = 'cell-border',
      rownames = FALSE,
      selection = 'none',
      escape = FALSE, # to allow HTML in table
      options = list(
        columnDefs = list(list(className = 'dt-nopad sim-cell', height=38, width=120, targets = 0),
                          list(targets = c(4,5,6,7), render = DT::JS(
                            "function(data, type, row, meta) {",
                            "return type === 'display' && data !== null && data.length > 17 ?",
                            "'<span title=\"' + data + '\">' + data.substr(0, 17) + '...</span>' : data;",
                            "}"))),
        ordering=FALSE,
        scrollX = TRUE,
        pageLength = 20,
        paging = TRUE,
        bInfo = 0,
        dom = 'ftp'
      )
    )
  }, server = TRUE)

  # proxy used to replace data
  # low priority to make sure data has been rendered
  proxy <- DT::dataTableProxy("query_table")
  observe({
    req(dummy_rendered())
    query_table <- query_table_final()
    DT::replaceData(proxy, query_table, rownames = FALSE)
  })
}
