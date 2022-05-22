#' Logic for Drugs Tab
#'
#' @inheritParams bulkPage
#' @inheritParams run_dseqr
#'
#' @return Called with \link[shiny]{callModule} to generate logic for
#'   single-cell tab.
#'
#' @export
drugsPage <- function(input, output, session, project_dir, pert_query_dir, pert_signature_dir, tx2gene_dir) {

  # the form area inputs/results
  form <- callModule(drugsForm, 'form',
                     project_dir = project_dir,
                     pert_query_dir = pert_query_dir,
                     pert_signature_dir = pert_signature_dir,
                     tx2gene_dir = tx2gene_dir)

  callModule(drugsGenesPlotly, 'genes',
             project_dir = project_dir,
             top_table = form$top_table,
             pert_signature = form$pert_signature,
             drug_study = form$drug_study)


  # the output table
  callModule(drugsTable, 'table',
             query_res = form$query_res,
             sorted_query = form$sorted_query,
             drug_study = form$drug_study,
             anal_name = form$anal_name,
             cells = form$cells,
             sort_by = form$sort_by,
             show_clinical = form$show_clinical,
             min_signatures = form$min_signatures,
             is_pert = form$is_pert,
             direction = form$direction)


}

# Logic for form on drugs page
#'
#' @keywords internal
#' @noRd
drugsForm <- function(input, output, session, project_dir, pert_query_dir, pert_signature_dir, tx2gene_dir) {

  new_custom <- reactiveVal()

  # dataset/analysis choices
  choices <- reactive({
    # reactive to new custom query, or bulk change (e.g. number of SVs)
    new_custom()

    project_dir <- project_dir()
    scseq_datasets <- load_scseq_datasets(project_dir)
    bulk_datasets <- load_bulk_datasets(project_dir, with.explore = TRUE)
    custom_anals <- load_custom_anals(project_dir)
    pert_anals <- load_pert_anals()

    choices <- rbind(bulk_datasets, scseq_datasets, custom_anals, pert_anals)
    choices$value <- seq_len(nrow(choices))+1
    choices$title <- choices$dataset_name
    choices$label <- stringr::str_trunc(choices$label, 35)
    choices <- rbind(NA, choices)

    return(choices)
  })

  # the selected dataset/analysis results
  selectedAnal <- callModule(selectedAnal, 'drugs',
                             choices = choices,
                             project_dir = project_dir,
                             tx2gene_dir = tx2gene_dir,
                             new_custom = new_custom,
                             pert_query_dir = pert_query_dir)





  # if currently selected analysis is custom then show genes for query
  custom_query <- callModule(customQueryForm, 'custom-query',
                             show_custom = selectedAnal$show_custom,
                             is_custom = selectedAnal$is_custom,
                             anal_name = selectedAnal$name,
                             new_custom = new_custom,
                             project_dir = project_dir)


  drugStudy <- callModule(selectedDrugStudy, 'drug_study',
                          drug_queries = selectedAnal$drug_queries,
                          is_pert = selectedAnal$is_pert)


  advancedOptions <- callModule(advancedOptions, 'advanced',
                                drug_study = drugStudy$drug_study,
                                show_advanced = drugStudy$show_advanced)

  pertSignature <- callModule(selectedPertSignature, 'genes',
                              project_dir = project_dir,
                              query_res = drugStudy$query_res,
                              query_type = drugStudy$query_type,
                              pert_signature_dir = pert_signature_dir)

  # show
  have_queries <- reactive(isTruthy(selectedAnal$drug_queries()))
  have_query <- reactive(isTruthy(drugStudy$query_res()))

  observe(toggle('drug_study_container', condition = have_queries()))
  observe(toggle('pert_signature_container', condition = have_query() & !selectedAnal$is_pert()))

  return(list(
    top_table = selectedAnal$top_table,
    is_pert = selectedAnal$is_pert,
    anal_name = selectedAnal$name,
    query_res = drugStudy$query_res,
    drug_study = drugStudy$drug_study,
    show_clinical = drugStudy$show_clinical,
    direction = drugStudy$direction,
    cells = advancedOptions$cells,
    sort_by = advancedOptions$sort_by,
    min_signatures = advancedOptions$min_signatures,
    pert_signature = pertSignature$pert_signature,
    sorted_query = pertSignature$sorted_query
  ))
}

#' Logic for custom query form on Drugs page
#'
#' @keywords internal
#' @noRd
customQueryForm <- function(input, output, session, show_custom, is_custom, anal_name, new_custom, project_dir) {
  input_ids <- c('custom_name', 'click_custom')

  # show hide custom query signature stuff
  observe({
    shinyjs::toggle('custom_query_container', anim = TRUE, condition = show_custom())
  })


  res_paths <- reactive({
    custom_name <- input$custom_name
    custom_dir <- file.path(project_dir(), 'custom_queries')
    if (!dir.exists(custom_dir)) dir.create(custom_dir)

    res_paths <- get_drug_paths(custom_dir, custom_name)
    res_paths$query_genes <- file.path(custom_dir, 'drugs', paste0('query_genes_', custom_name, '.qs'))
    return(res_paths)
  })


  # Custom upload
  error_msg <- reactiveVal()


  observe({
    msg <- error_msg()
    html('error_msg', html = msg)
    shinyjs::toggleClass('validate', 'has-error', condition = isTruthy(msg))
  })

  observeEvent(input$click_custom, {
    if (!isTruthy(input$custom_name)) {
      error_msg('Need name for query')
      return(NULL)
    }
    error_msg(NULL)
    shinyjs::click('up_custom')
  })


  observeEvent(input$up_custom, {


    infile <- input$up_custom
    if (!isTruthy(infile)){
      msg <- NULL

    } else {
      top_table <- tryCatch(
        utils::read.csv(infile$datapath, check.names = FALSE, stringsAsFactors = FALSE, row.names = 1),
        error = function(e) return(NULL))

      msg <- validate_up_custom(top_table, input$custom_name)

      if (is.null(msg)) {
        shinyjs::removeClass('validate', class = 'has-error')

        # slowest part is loading es
        disableAll(input_ids)

        progress <- Progress$new(session, min = 0, max = 2)
        progress$set(message = "Querying drugs", value = 1)
        on.exit(progress$close())

        top_table <- format_up_custom(top_table)
        res_paths <- res_paths()
        es <- load_drug_es()
        progress$inc(1)

        run_drug_queries(top_table, res_paths, es, ngenes = nrow(top_table))
        progress$inc(1)

        qs::qsave(top_table, res_paths$query_genes)
        new_custom(input$custom_name)
        enableAll(input_ids)
      }
    }

    error_msg(msg)
  })

}

#' Logic for selected drug study in drugsForm
#'
#' @keywords internal
#' @noRd
selectedDrugStudy <- function(input, output, session, drug_queries, is_pert) {


  # toggle display of perturbation study inputs
  have_queries <- reactive(isTruthy(drug_queries()))

  # show advanced options,  clinical only results, and gene plots
  show_advanced <- reactive(input$advanced %% 2 != 0)
  show_clinical <- reactive(input$clinical %% 2 != 1)
  show_genes <- reactive(input$show_genes %% 2 != 0)

  observe({
    shinyjs::toggleClass('advanced', 'btn-primary', condition = show_advanced())
    shinyjs::toggleClass('clinical', 'btn-primary', condition = show_clinical())
    shinyjs::toggleClass('show_genes', 'btn-primary', condition = show_genes())
  })

  # show pert study choices
  choices <- data.frame(study = c(NA, 'CMAP02', 'L1000', 'L1000'),
                        subset = c(NA, 'drugs', 'drugs', 'genetic'),
                        value = c(NA, 'CMAP02', 'L1000 Drugs', 'L1000 Genetic'),
                        stringsAsFactors = FALSE)



  observe({
    updateSelectizeInput(session, 'study', choices = choices,
                         options = list(render = I('{option: studyOption, item: studyItem}')), server = TRUE)
  })


  # toggle correlation direction filter and icon
  direction <- reactive({
    switch((input$direction %% 3) + 1,
           'both',
           'similar',
           'opposing')
  })

  direction_icon <- reactive({
    switch(direction(),
           'both' = 'arrows-alt-v',
           'similar' = 'chevron-up',
           'opposing' = 'chevron-down')
  })

  observe(updateActionButton(session, 'direction', icon = icon(direction_icon(), 'fa-fw')))

  # disable buttons based on genetic/pert
  is_genetic <- reactive(input$study == 'L1000 Genetic')

  # sort by absolute if either is genetic or is CMAP/L1000 pert
  observe({
    toggle('direction-parent', condition = is_genetic() | is_pert())
    toggle('clinical-parent', condition = !is_genetic())
  })

  # order pert choices same as in drugs
  query_type <- reactive(switch(input$study,
                                'CMAP02' = 'cmap',
                                'L1000 Drugs' = 'l1000_drugs',
                                'L1000 Genetic' = 'l1000_genes'))

  # query result for selected perturbation study
  query_res <- reactive({
    type <- query_type()
    drug_queries <- drug_queries()
    if (is.null(drug_queries) | is.null(type)) return(NULL)
    drug_queries[[type]]
  })



  return(list(
    drug_study = reactive(input$study),
    query_res = query_res,
    query_type = query_type,
    show_clinical = show_clinical,
    show_advanced = show_advanced,
    show_genes = show_genes,
    direction = direction
  ))

}

#' Logic for advanced options in drugsForm
#'
#' @keywords internal
#' @noRd
advancedOptions <- function(input, output, session, drug_study, show_advanced) {

  # update choices for cell lines based on selected study
  cell_choices <- reactive({
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

#' Logic for showing pert selection for gene plot
#'
#' @keywords internal
#' @noRd
selectedPertSignature <- function(input, output, session, project_dir, query_res, query_type, pert_signature_dir) {
  pert_options <- list(render = I('{option: pertOptions, item: pertItem}'))

  sorted_query <- reactive({
    q <- query_res()
    type <- query_type()
    if (is.null(type) | is.null(q)) return(NULL)
    if(type == 'l1000_genes') q[order(abs(q), decreasing = TRUE)] else sort(q)
  })

  pert_choices <- reactive({
    query <- sorted_query()
    req(query)

    pert_choices <- data.frame(
      label = gsub('_-700-666.0_\\d+h$', '', names(query)),
      value = names(query),
      cor = format(round(query, digits = 3)), stringsAsFactors = FALSE
    )
    return(pert_choices)
  })

  # update pert signature choices
  observe({
    updateSelectizeInput(session,
                         'pert',
                         choices = rbind(rep(NA, 3), pert_choices()),
                         options = pert_options,
                         server = TRUE)
  })

  # load signature for pert
  # allow empty signature if just want to show query genes
  pert_signature <- reactive({
    pert <- input$pert
    if (is.null(pert) || pert == '') return(NULL)
    load_pert_signature(pert, query_type(), pert_signature_dir)
  })



  return(list(
    pert_signature = pert_signature,
    sorted_query = sorted_query
  ))


}

#' Logic for drug table
#'
#' @keywords internal
#' @importFrom magrittr "%>%"
#' @noRd
drugsTable <- function(input, output, session, query_res, sorted_query, drug_study, anal_name, cells, show_clinical, sort_by, min_signatures, is_pert, direction, ntop = 1500) {
  pert_options <- list(render = I('{option: pertOptions, item: pertItem}'))
  cmap_annot <- l1000_drugs_annot <- l1000_genes_annot <- NULL

  dummy_rendered <- reactiveVal(FALSE)

  # get either cmap or l1000 annotations

  drug_annot <- reactive({
    drug_study <- drug_study()
    req(drug_study)

    # update globals when first use
    if (drug_study == 'CMAP02') {
      if (is.null(cmap_annot)) cmap_annot <- get_drugs_table('CMAP02')
      return(cmap_annot)

    } else if (drug_study == 'L1000 Drugs') {
      if (is.null(l1000_drugs_annot)) l1000_drugs_annot <- get_drugs_table('L1000_drugs')
      return(l1000_drugs_annot)

    } else if (drug_study == 'L1000 Genetic') {
      if (is.null(l1000_genes_annot)) l1000_genes_annot <- get_drugs_table('L1000_genes')
      return(l1000_genes_annot)
    }
  })

  # add annotations to query result
  query_table_annot <- reactive({
    query_res <- query_res()
    if (is.null(query_res)) return(NULL)

    drug_annot <- drug_annot()
    req(drug_annot)

    annot_query_res(query_res, drug_annot)
  })


  # sort by absolute if either is genetic or is CMAP/L1000 pert
  is_genetic <- reactive(drug_study() == 'L1000 Genetic')
  sort_abs <- reactive(is_genetic() | is_pert())

  # subset to selected cells, summarize by compound, and add html
  query_table_summarised <- reactive({
    query_table_annot <- query_table_annot()
    if (is.null(query_table_annot)) return(NULL)

    summarise_query_table(query_table_annot,
                          is_genetic = is_genetic(),
                          cells = cells(),
                          sort_abs = sort_abs(),
                          ntop = ntop)
  })



  # build up to final table
  # ---

  # subset by min signatures
  query_table_nsig <- reactive({
    query_table_summarised <- query_table_summarised()
    if (is.null(query_table_summarised)) return(NULL)

    filter_nsig(query_table_summarised, min_signatures())
  })

  # subset by clinical phase
  query_table_clin <- reactive({
    query_table_nsig <- query_table_nsig()
    if (is.null(query_table_nsig)) return(NULL)

    filter_clinical(query_table_nsig, show_clinical())
  })

  # final sorting/filtering
  query_table_final <- reactive({
    query_table_clin <- query_table_clin()
    if (is.null(query_table_clin)) return(NULL)

    sort_query_table_clin(query_table_clin(),
                          sort_by=sort_by(),
                          sort_abs=sort_abs(),
                          direction=direction(),
                          drug_study=drug_study())
  })

  # will update with proxy to prevent redraw
  dummy_table <- reactive({
    study <- drug_study()
    req(study, query_res())
    dummy_rendered(FALSE)
    cols <- get_query_cols(study == 'L1000 Genetic')

    data.frame(matrix(ncol = length(cols), dimnames = list(NULL, cols)), check.names = FALSE)
  })


  # query table for downloading
  query_table_dl <- reactive({

    q <- sorted_query()
    if (is.null(q)) return(NULL)
    data.frame(correlation = q, signature = names(q), row.names = NULL)
  })

  # show query data
  output$query_table <- DT::renderDataTable({
    # ellipses for wide columns
    dummy_table <- dummy_table()
    wide_cols <- c('MOA', 'Target', 'Disease Area', 'Indication', 'Vendor', 'Catalog #', 'Vendor Name')
    # -1 needed with rownames = FALSE
    elipsis_targets <- which(colnames(dummy_table) %in% wide_cols) - 1
    hide_target <- which(colnames(dummy_table) %in% c('Vendor', 'Catalog #', 'Vendor Name')) - 1

    # don't show column visibility button for genetic
    dom <- ifelse(is_genetic(), 'ftp', 'Bfrtip')

    dummy_rendered(TRUE)

    DT::datatable(
      dummy_table,
      class = 'cell-border',
      rownames = FALSE,
      selection = 'none',
      escape = FALSE, # to allow HTML in table
      options = list(
        columnDefs = list(list(className = 'dt-nopad sim-cell', height=38, width=120, targets = 1),
                          list(className = 'dt-rank', targets = 0, width=50),
                          list(targets = hide_target, visible=FALSE),
                          list(targets = elipsis_targets, render = DT::JS(
                            "function(data, type, row, meta) {",
                            "return type === 'display' && data !== null && data.length > 27 ?",
                            "'<span title=\"' + data + '\">' + data.substr(0, 27) + '...</span>' : data;",
                            "}"))
        ),
        drawCallback = DT::JS('function(setting, json) { setupContextMenu(); }'),
        ordering = FALSE,
        scrollX = TRUE,
        pageLength = 50,
        paging = TRUE,
        bInfo = 0,
        dom = dom
      )
    )
  },
  server = TRUE)

  # proxy used to replace data
  # low priority to make sure data has been rendered
  proxy <- DT::dataTableProxy("query_table")
  observe({
    req(dummy_rendered())
    query_table <- query_table_final()
    DT::replaceData(proxy, query_table, rownames = FALSE)
  })


  # download link for table
  dl_ready <- reactive(isTruthy(query_table_dl()))

  observe(toggle('dl_drugs', condition = dl_ready()))

  output$dl_drugs <- downloadHandler(
    filename = function() {
      drug_study <- drug_study()
      anal_name <- anal_name()
      req(anal_name)
      drug_study <- tolower(drug_study)
      drug_study <- gsub(' ', '_', drug_study)
      paste0(anal_name, '_', drug_study, '.csv')
    },
    content = function(con) {utils::write.csv(query_table_dl(), con, row.names = FALSE)}
  )

}

#' Logic query/drug genes plotly
#'
#' @keywords internal
#' @importFrom magrittr "%>%"
#' @noRd
drugsGenesPlotly <- function(input, output, session, project_dir, top_table, drug_study, pert_signature) {

  path_id <- reactive({
    study <- drug_study()
    if (study == 'CMAP02') return('CMAP02')
    if (grepl('^L1000', study)) return('L1000')
  })

  # the gene plot
  pl <- reactive({
    top_table <- top_table()
    path_id <- path_id()
    if (is.null(path_id) | is.null(top_table)) return(NULL)

    pert_signature <- pert_signature()
    path_df <- get_path_df(top_table, path_id, pert_signature)
    if (nrow(path_df) == 0) return(NULL)

    plot_dprimes(path_df)

  })

  observe(toggle('container', condition = isTruthy(pl())))

  output$plotly <- snapshotPreprocessOutput(
    plotly::renderPlotly({
      pl()
    }),
    function(value) { 'genes_plotly' }
  )
}


#' Logic for selected dataset/analysis in Drugs tab
#'
#' @keywords internal
#' @noRd
selectedAnal <- function(input, output, session, project_dir, choices, new_custom, tx2gene_dir, pert_query_dir = NULL) {
  options <- list(render = I('{option: querySignatureOptions, item: querySignatureItem}'),
                  searchField = c('dataset_name', 'label'))


  # update dataset/analysis choice
  observe({
    choices <- choices()
    req(choices)
    updateSelectizeInput(session, 'query', choices = choices, server = TRUE, options = options)
  })

  observe({
    new_custom <- new_custom()
    req(new_custom)

    choices <- isolate(choices())
    new_custom <- which(choices$label == new_custom)
    updateSelectizeInput(session, 'query', choices = choices, server = TRUE, options = options, selected = new_custom)
  })

  # right click load signature logic for Drugs tab
  runjs(paste0('initContextMenu("', session$ns('pert_query_name_load'), '", "', session$ns('pert_query_name_show'), '");'))
  observe({
    sel <- input$pert_query_name_load
    req(sel)

    choices <- choices()
    sel_idx <- which(choices$label == sel)

    updateSelectizeInput(session, 'query', choices = choices, selected = sel_idx, server = TRUE)
  })

  # the selected dataset/analysis
  sel <- reactive({
    row_num <- input$query
    choices <- choices()
    req(choices)
    if (!isTruthy(row_num)) return(list(group = ''))
    req(row_num, choices)

    choices[row_num, ]
  })

  # the type of dataset/analysis
  is_bulk <- reactive(sel()$group == 'Bulk Data')
  is_sc <- reactive(sel()$group == 'Single Cell')
  is_custom <- reactive(sel()$group == 'Custom')
  is_pert <- reactive(sel()$group == 'CMAP02/L1000 Perturbations')


  observe({
    shinyjs::toggle('sc_clusters_container', condition = is_sc())
    shinyjs::toggle('bulk_groups_container', condition = is_bulk())
  })


  # show/hide custom signature inputs in Drugs tab
  show_custom <- reactive(input$show_custom %% 2 != 0)

  observe({
    shinyjs::toggleClass(id = "show_custom", 'btn-primary', condition = show_custom())
  })

  sel_name <- reactive(sel()$dataset_name)

  dataset_dir <- reactive({
    if (!isTruthy(input$query)) return(NULL)
    dataset_dir <- sel()$dataset_dir
    if (is.na(dataset_dir)) return(NULL)

    file.path(project_dir(), dataset_dir)
  })


  # Bulk analysis
  # ---
  eset  <- reactive(qread.safe(file.path(dataset_dir(), 'eset.qs')))
  pdata <- reactive(qread.safe(file.path(dataset_dir(), 'pdata_explore.qs')))
  svobj <- reactive(qread.safe(file.path(dataset_dir(), 'svobj.qs')))

  numsv_path <- reactive(file.path(dataset_dir(), 'numsv.qs'))
  numsv <- reactive({
    numsv_path <- numsv_path()
    if (!file.exists(numsv_path)) qs::qsave(0, numsv_path)
    qs::qread(numsv_path)
  })

  bulk_eset <- exploreEset(eset = eset,
                           dataset_dir = dataset_dir,
                           explore_pdata = pdata,
                           numsv = numsv,
                           svobj = svobj)

  bulkAnal <- callModule(bulkAnal, 'bulk',
                         pdata = pdata,
                         eset = bulk_eset,
                         svobj = svobj,
                         numsv = numsv,
                         dataset_dir = dataset_dir,
                         dataset_name = dataset_name,
                         is_bulk = is_bulk)




  # Single cell analysis
  # ---

  resoln <- reactive({
    dataset_dir <- dataset_dir()
    req(dataset_dir)
    resoln_path <- file.path(dataset_dir, 'resoln.qs')
    qread.safe(resoln_path)
  })

  resoln_dir <- reactive({
    resoln <- resoln()
    if (is.null(resoln)) return(NULL)
    file.path(dataset_dir(), get_resoln_dir(resoln))
  })


  dataset_name <- reactive(sel()$dataset_name)



  scSampleGroups <- callModule(scSampleGroups, 'sample_groups',
                               dataset_dir = dataset_dir,
                               resoln_dir = resoln_dir,
                               dataset_name = dataset_name)

  scSampleClusters <- callModule(scSampleClusters, 'sample_clusters',
                                 input_scseq = scSampleGroups$scseq,
                                 meta = scSampleGroups$meta,
                                 lm_fit = scSampleGroups$lm_fit,
                                 groups = scSampleGroups$groups,
                                 dataset_dir = dataset_dir,
                                 resoln_dir = resoln_dir,
                                 tx2gene_dir = tx2gene_dir,
                                 resoln = resoln,
                                 dataset_name = dataset_name,
                                 page = 'drugs')



  # results for selected type (bulk, sc, custom, pert)
  drug_queries <- reactive({
    sel_name <- sel_name()
    if (is_sc()) {
      drug_queries <- scSampleClusters$drug_queries()

    } else if (is_bulk()) {
      drug_queries <- bulkAnal$drug_queries()

    } else if (is_custom()) {
      custom_dir <- file.path(project_dir(), 'custom_queries')
      drug_paths <- get_drug_paths(custom_dir, sel_name)
      drug_queries <- lapply(drug_paths, function(x) if (file.exists(x)) qs::qread(x))

    } else if (is_pert()) {
      drug_paths <- get_drug_paths(pert_query_dir, fs::path_sanitize(sel_name), ftype = '.rds')
      sapply(drug_paths, dl_pert_result)
      drug_queries <- lapply(drug_paths, readRDS)

    } else {
      drug_queries <- NULL
    }

    return(drug_queries)
  })

  top_table <- reactive({
    if (is_sc()) {
      top_table <- scSampleClusters$top_table()[[1]]

    } else if (is_bulk()) {
      top_table <- bulkAnal$top_table()

    } else if (is_custom()) {
      fname <- paste0('query_genes_', sel_name(), '.qs')
      top_table <- qs::qread(file.path(project_dir(), 'custom_queries', 'drugs', fname))

    } else {
      top_table <- NULL
    }
    return(top_table)
  })


  path_res <- reactive({
    if (is_sc()) {
      path_res <- scSampleClusters$path_res()

    } else if (is_bulk()) {
      path_res <- bulkAnal$path_res()

    } else {
      path_res <- NULL
    }
    return(path_res)
  })

  anal_name <- reactive({
    sel_name <- sel_name()
    if (is_sc()) {
      anal_name <- paste(sel_name, scSampleClusters$annot_clusters(), sep = '_')

    } else if (is_bulk()) {
      anal_name <- paste(sel_name, bulkAnal$name(), sep = '_')

    } else {
      anal_name <- sel_name
    }
    return(anal_name)
  })


  return(list(
    name = anal_name,
    bulk_eset = bulk_eset,
    contrast_groups = bulkAnal$contrast_groups,
    top_table = top_table,
    path_res = path_res,
    show_custom = show_custom,
    drug_queries = drug_queries,
    is_bulk = is_bulk,
    is_sc = is_sc,
    is_custom = is_custom,
    is_pert = is_pert
  ))

}
