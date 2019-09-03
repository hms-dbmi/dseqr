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
             show_clinical = form$show_clinical,
             min_signatures = form$min_signatures)

}

#' Logic for custom query form on Drugs page
#' @export
#' @keywords internal
customQueryForm <- function(input, output, session, show_custom, is_custom, anal, new_anal, data_dir) {

  # setup choices and options
  choices <- data.frame(gene = c(genes$common, genes$cmap_only), stringsAsFactors = FALSE)
  choices$value <- choices$label <- choices$gene
  choices$cmap_only <- choices$gene %in% genes$cmap_only
  options <- list(render = I('{option: queryGenesOption, item: queryGenesItem}'))


  # show hide custom query signature stuff
  observe({
    shinyjs::toggle('custom_query_container', anim = TRUE, condition = show_custom())
  })


  # update genes to downregulate
  observe({

    # if current signature is custom, then load and show previously selected genes
    if (is_custom()) {
      fname <- paste0('query_genes_', anal()$anal_name, '.rds')
      query_genes <- readRDS(file.path(data_dir, 'custom_queries', fname))
      sel_up <- query_genes$up
      sel_dn <- query_genes$dn

    } else {
      sel_up <- sel_dn <- NULL
    }

    updateSelectizeInput(session, 'dn_genes', choices = choices, selected = sel_dn, options = options, server = TRUE)
    updateSelectizeInput(session, 'up_genes', choices = choices, selected = sel_up, options = options, server = TRUE)
  })

  res_paths <- reactive({
    custom_name <- input$custom_name
    custom_dir <- file.path(data_dir, 'custom_queries')
    if (!dir.exists(custom_dir)) dir.create(custom_dir)

    list(
      cmap = file.path(custom_dir, paste0('cmap_res_', custom_name, '.rds')),
      l1000_drugs = file.path(custom_dir, paste0('l1000_drugs_res_', custom_name, '.rds')),
      l1000_genes = file.path(custom_dir, paste0('l1000_genes_res_', custom_name, '.rds')),
      query_genes = file.path(custom_dir, paste0('query_genes_', custom_name, '.rds'))
    )
  })


  observeEvent(input$submit_custom, {

    error_msg <- validate_custom_query(dn_genes = input$dn_genes,
                                       up_genes = input$up_genes,
                                       custom_name = input$custom_name)

    if (is.null(error_msg)) {
      removeClass('validate', class = 'has-error')

      query_genes <- list(dn = input$dn_genes, up = input$up_genes)
      res <- run_custom_query(query_genes = query_genes,
                              res_paths = res_paths(),
                              session = session)

      new_anal(input$custom_name)


    } else {
      html('error_msg', html = error_msg)
      addClass('validate', class = 'has-error')
    }
  })
}


# Logic for form on drugs page
#' @export
#' @keywords internal
drugsForm <- function(input, output, session, new_anal, data_dir) {

  cmap_res <- reactiveVal()
  l1000_drugs_res <- reactiveVal()
  l1000_genes_res <- reactiveVal()

  querySignature <- callModule(querySignature, 'signature',
                               new_anal = new_anal,
                               data_dir = data_dir)

  # show/hide single cell stuff
  is_sc <- reactive({
    anal <- querySignature$anal()
    anal$type == 'Single Cell'
  })

  observe({
    shinyjs::toggle('sc_clusters_container', condition = is_sc())
  })

  input_ids <- c('run_comparison', 'selected_clusters')
  sc_inputs <- scSampleComparison(input, output, session,
                                  data_dir = data_dir,
                                  anal = querySignature$anal,
                                  is_sc = is_sc,
                                  input_ids = input_ids,
                                  with_drugs = TRUE)


  # if currently selected analysis is custom then show genes
  is_custom <- reactive({
    anal <- querySignature$anal()
    anal$type == 'Custom'
  })

  custom_query <- callModule(customQueryForm, 'custom-query',
                             show_custom = querySignature$show_custom,
                             is_custom = is_custom,
                             anal = querySignature$anal,
                             new_anal = new_anal,
                             data_dir = data_dir)


  # get saved cmap/l1000 query results
  observe({
    disable('signature')

    if (is_sc())  {
      res <- sc_inputs$results()

    } else if (is_custom()) {
      res_paths <- querySignature$res_paths()
      res <- load_custom_results(res_paths)

    } else {
      res_paths <- querySignature$res_paths()
      res <- run_drugs_comparison(res_paths, session)
    }


    cmap_res(res$cmap)
    l1000_drugs_res(res$l1000_drugs)
    l1000_genes_res(res$l1000_genes)

    enable('signature')
  })


  drugStudy <- callModule(selectedDrugStudy, 'drug_study',
                          anal = querySignature$anal)



  advancedOptions <- callModule(advancedOptions, 'advanced',
                                cmap_res = cmap_res,
                                l1000_res = l1000_res,
                                drug_study = drugStudy$drug_study,
                                show_advanced = drugStudy$show_advanced)

  query_res <- reactive({
    drug_study <- drugStudy$drug_study()
    cmap_res <- cmap_res()
    l1000_drugs_res <- l1000_drugs_res()
    l1000_genes_res <- l1000_genes_res()

    if (drug_study == 'CMAP02') return(cmap_res)
    else if (drug_study == 'L1000 Drugs') return(l1000_drugs_res)
    else if (drug_study == 'L1000 Genetic') return(l1000_genes_res)

    return(NULL)
  })




  return(list(
    query_res = query_res,
    drug_study = drugStudy$drug_study,
    cells = advancedOptions$cells,
    sort_by = advancedOptions$sort_by,
    show_clinical = drugStudy$show_clinical,
    min_signatures = advancedOptions$min_signatures
  ))


}


#' Logic for query signature in drugsForm
#' @export
#' @keywords internal
querySignature <- function(input, output, session, new_anal, data_dir) {


  # reload query choices if new analysis
  anals <- reactive({
    new_anal()

    scseq_anals <- load_scseq_anals(data_dir, with_type = TRUE)
    bulk_anals <- load_bulk_anals(data_dir, with_type = TRUE)
    custom_anals <- load_custom_anals(data_dir)

    anals <- rbind(bulk_anals, scseq_anals, custom_anals)
    anals$value <- seq_len(nrow(anals))

    return(anals)
  })

  observe({
    anals <- anals()
    req(anals)
    updateSelectizeInput(session, 'query', choices = anals, server = TRUE)
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
    req(anal)

    dataset_dir <- ifelse(anal$type == 'Custom',
                          file.path(data_dir, 'custom_queries'),
                          file.path(data_dir, 'bulk', anal$dataset_dir))

    anal_name <- anal$anal_name

    list(
      anal = file.path(dataset_dir, paste0('diff_expr_symbol_', anal_name, '.rds')),
      cmap = file.path(dataset_dir, paste0('cmap_res_', anal_name, '.rds')),
      l1000_drugs = file.path(dataset_dir, paste0('l1000_drugs_res_', anal_name, '.rds')),
      l1000_genes = file.path(dataset_dir, paste0('l1000_genes_res_', anal_name, '.rds'))
    )
  })


  show_custom <- reactive(input$show_custom %% 2 != 0)

  # show/hide integration form
  observe({
    toggleClass(id = "show_custom", 'btn-primary', condition = show_custom())
  })



  return(list(
    res_paths = res_paths,
    anal = anal,
    show_custom = show_custom
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
    choices <- data.frame(study = c('CMAP02', 'L1000', 'L1000'),
                          subset = c('drugs', 'drugs', 'genetic'),
                          value = c('CMAP02', 'L1000 Drugs', 'L1000 Genetic'),
                          stringsAsFactors = FALSE)

    updateSelectizeInput(session, 'study', choices = choices, selected = NULL, options = list(render = I('{option: studyOption, item: studyItem}')),  server = TRUE)
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

  # update choices for cell lines based on selected study
  cell_choices <- shiny::reactive({
    req(drug_study())
    get_cell_choices(drug_study())
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
    sort_by = reactive(input$sort_by),
    min_signatures = reactive(input$min_signatures)
  ))

}

#' Logic for drug table
#' @export
#' @keywords internal
#' @importFrom magrittr "%>%"
drugsTable <- function(input, output, session, query_res, drug_study, cells, show_clinical, sort_by, min_signatures) {
  drug_cols <- c('Correlation', 'Compound', 'Clinical Phase', 'External Links', 'MOA', 'Target', 'Disease Area', 'Indication', 'Vendor', 'Catalog #', 'Vendor Name', 'n')
  gene_cols <- c('Correlation', 'Compound', 'External Links', 'Description', 'n')

  runjs('initContextMenu();')

  dummy_rendered <- reactiveVal(FALSE)

  # get either cmap or l1000 annotations
  drug_annot <- reactive({
    drug_study <- drug_study()
    req(drug_study)

    # update globals when first use
    if (is.null(cmap_annot)) {
      cmap_annot <<- get_drugs_table('CMAP02')
      l1000_drugs_annot <<- get_drugs_table('L1000_drugs')
      l1000_genes_annot <<- get_drugs_table('L1000_genes')
    }

    if (drug_study == 'CMAP02') return(cmap_annot)
    else if (drug_study == 'L1000 Drugs') return(l1000_drugs_annot)
    else if (drug_study == 'L1000 Genetic') return(l1000_genes_annot)
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

    is_genetic <- grepl('Genetic', drug_study())
    cols <- if(is_genetic) gene_cols else drug_cols

    query_table <- query_table_full %>%
      limit_cells(cells()) %>%
      summarize_compound(is_genetic = is_genetic) %>%
      add_table_html() %>%
      select(cols, everything())
  })

  query_table_final <- reactive({
    query_table <- query_table_summarised()
    if (is.null(query_table)) return(NULL)
    sort_by <- sort_by()
    drug_study <- drug_study()

    query_table <- dplyr::filter(query_table, n >= min_signatures())

    # subset by clinical phase
    if (show_clinical() && 'Clinical Phase' %in% colnames(query_table))
      query_table <- dplyr::filter(query_table, !is.na(`Clinical Phase`))

    if (sort_by == 'avg_cor') {
      query_table$Correlation <- gsub('simplot', 'simplot show-meanline', query_table$Correlation)
    }

    # show largest absolute correlations first for genetic
    # as both directions are informative
    if (drug_study == 'L1000 Genetic') {
      query_table <- query_table %>%
        dplyr::mutate(min_cor = -pmax(abs(min_cor), abs(max_cor))) %>%
        dplyr::mutate(avg_cor = -abs(avg_cor))
    }

    # sort as desired
    query_table %>%
      dplyr::arrange(!!sym(sort_by)) %>%
      dplyr::select(-min_cor, -avg_cor, -max_cor)
  })

  # will update with proxy to prevent redraw
  dummy_table <- reactive({
    study <- drug_study()
    dummy_rendered(FALSE)
    cols <- if (study == 'L1000 Genetic') gene_cols else drug_cols

    data.frame(matrix(ncol = length(cols), dimnames = list(NULL, cols)), check.names = FALSE)
  })


  # show query data
  output$query_table <- DT::renderDataTable({
    # ellipses for wide columns
    dummy_table <- dummy_table()
    wide_cols <- c('MOA', 'Target', 'Disease Area', 'Indication', 'Vendor', 'Catalog #', 'Vendor Name')
    # -1 needed with rownames = FALSE
    elipsis_targets <- which(colnames(dummy_table) %in% wide_cols) - 1
    hide_target <- which(colnames(dummy_table) %in% c('n')) - 1

    dummy_rendered(TRUE)

    DT::datatable(
      dummy_table,
      class = 'cell-border',
      rownames = FALSE,
      selection = 'none',
      escape = FALSE, # to allow HTML in table
      options = list(
        columnDefs = list(list(className = 'dt-nopad sim-cell', height=38, width=120, targets = 0),
                          list(targets = hide_target, visible=FALSE),
                          list(targets = elipsis_targets, render = DT::JS(
                            "function(data, type, row, meta) {",
                            "return type === 'display' && data !== null && data.length > 17 ?",
                            "'<span title=\"' + data + '\">' + data.substr(0, 17) + '...</span>' : data;",
                            "}"))
        ),
        drawCallback = DT::JS('function(setting, json) { setupContextMenu(); }'),
        ordering = FALSE,
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
    runjs('setupContextMenu();')
  })
}

