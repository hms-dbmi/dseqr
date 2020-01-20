#' Logic for Drugs page
#' @export
#' @keywords internal
drugsPage <- function(input, output, session, new_bulk, data_dir, pert_query_dir, pert_signature_dir) {

  # the form area inputs/results
  form <- callModule(drugsForm, 'form',
                     new_bulk = new_bulk,
                     data_dir = data_dir,
                     pert_query_dir = pert_query_dir,
                     pert_signature_dir = pert_signature_dir)

  callModule(drugsGenesPlotly, 'genes',
             data_dir = data_dir,
             top_table = form$top_table,
             ambient = form$ambient,
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
#' @export
#' @keywords internal
drugsForm <- function(input, output, session, data_dir, new_bulk, pert_query_dir, pert_signature_dir) {

  # TODO: trigger new_custom after running custom query
  new_custom <- reactiveVal()

  # dataset/analysis choices
  choices <- reactive({
    # reactive to new datasets, new custom query, or bulk change (e.g. number of SVs)
    new_bulk()
    new_custom()
    scseq_datasets <- load_scseq_datasets(data_dir)
    bulk_datasets <- load_bulk_datasets(data_dir)
    custom_anals <- load_custom_anals(data_dir)
    pert_anals <- load_pert_anals()

    choices <- rbind(bulk_datasets, scseq_datasets, custom_anals, pert_anals)
    choices$value <- seq_len(nrow(choices))+1
    choices <- rbind(rep(NA, 5), choices)

    return(choices)
  })

  # the selected dataset/analysis results
  selectedAnal <- callModule(selectedAnal, 'drugs',
                             choices = choices,
                             data_dir = data_dir,
                             pert_query_dir = pert_query_dir)





  # if currently selected analysis is custom then show genes for query
  custom_query <- callModule(customQueryForm, 'custom-query',
                             show_custom = selectedAnal$show_custom,
                             is_custom = selectedAnal$is_custom,
                             anal_name = selectedAnal$anal_name,
                             new_anal = new_anal,
                             data_dir = data_dir)


  drugStudy <- callModule(selectedDrugStudy, 'drug_study',
                          drug_queries = selectedAnal$drug_queries,
                          is_pert = selectedAnal$is_pert)


  advancedOptions <- callModule(advancedOptions, 'advanced',
                                drug_study = drugStudy$drug_study,
                                show_advanced = drugStudy$show_advanced)

  pertSignature <- callModule(selectedPertSignature, 'genes',
                              data_dir = data_dir,
                              query_res = drugStudy$query_res,
                              query_type = drugStudy$query_type,
                              pert_signature_dir = pert_signature_dir)

  # show
  have_queries <- reactive(isTruthy(selectedAnal$drug_queries()))
  have_query <- reactive(isTruthy(drugStudy$query_res()))

  observe(toggle('drug_study_container', condition = have_queries()))
  observe(toggle('pert_signature_container', condition = have_query()))

  return(list(
    top_table = selectedAnal$top_table,
    ambient = selectedAnal$ambient,
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
#' @export
#' @keywords internal
customQueryForm <- function(input, output, session, show_custom, is_custom, anal_name, new_anal, data_dir) {

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
      fname <- paste0('query_genes_', anal_name(), '.rds')
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

    res_paths <- get_drug_paths(custom_dir, custom_name)
    res_paths$query_genes <- file.path(custom_dir, paste0('query_genes_', custom_name, '.rds'))
  })


  observeEvent(input$submit_custom, {

    error_msg <- validate_custom_query(dn_genes = input$dn_genes,
                                       up_genes = input$up_genes,
                                       custom_name = input$custom_name)

    if (is.null(error_msg)) {
      shinyjs::removeClass('validate', class = 'has-error')

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

#' Logic for selected drug study in drugsForm
#' @export
#' @keywords internal
selectedDrugStudy <- function(input, output, session, drug_queries, is_pert) {


  # toggle display of perturbation study inputs
  have_queries <- reactive(isTruthy(drug_queries()))

  # show advanced options,  clinical only results, and gene plots
  show_advanced <- reactive(input$advanced %% 2 != 0)
  show_clinical <- reactive(input$clinical %% 2 != 0)
  show_genes <- reactive(input$show_genes %% 2 != 0)

  observe({
    toggleClass('advanced', 'btn-primary', condition = show_advanced())
    toggleClass('clinical', 'btn-primary', condition = show_clinical())
    toggleClass('show_genes', 'btn-primary', condition = show_genes())
  })

  # show pert study choices
  choices <- data.frame(study = c('CMAP02', 'L1000', 'L1000'),
                        subset = c('drugs', 'drugs', 'genetic'),
                        value = c('CMAP02', 'L1000 Drugs', 'L1000 Genetic'),
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
#' @export
#' @keywords internal
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
#' @export
#' @keywords internal
selectedPertSignature <- function(input, output, session, data_dir, query_res, query_type, pert_signature_dir) {
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
#' @export
#' @keywords internal
#' @importFrom magrittr "%>%"
drugsTable <- function(input, output, session, query_res, sorted_query, drug_study, anal_name, cells, show_clinical, sort_by, min_signatures, is_pert, direction) {
  drug_cols <- c('Rank', 'Correlation', 'Compound', 'Clinical Phase', 'External Links', 'MOA', 'Target', 'Disease Area', 'Indication', 'Vendor', 'Catalog #', 'Vendor Name')
  gene_cols <- c('Rank', 'Correlation', 'Compound', 'External Links', 'Description')
  pert_options <- list(render = I('{option: pertOptions, item: pertItem}'))

  dummy_rendered <- reactiveVal(FALSE)

  # get either cmap or l1000 annotations

  drug_annot <- reactive({
    drug_study <- drug_study()
    req(drug_study)

    # update globals when first use
    if (drug_study == 'CMAP02') {
      if (is.null(cmap_annot)) cmap_annot <<- get_drugs_table('CMAP02')
      return(cmap_annot)

    } else if (drug_study == 'L1000 Drugs') {
      if (is.null(l1000_drugs_annot)) l1000_drugs_annot <<- get_drugs_table('L1000_drugs')
      return(l1000_drugs_annot)

    } else if (drug_study == 'L1000 Genetic') {
      if (is.null(l1000_genes_annot)) l1000_genes_annot <<- get_drugs_table('L1000_genes')
      return(l1000_genes_annot)
    }
  })

  # add annotations to query result
  query_table_annot <- reactive({
    query_res <- query_res()
    req(query_res)

    drug_annot <- drug_annot()
    req(query_res, drug_annot)
    drug_annot <- drug_annot[drug_annot$title %in% names(query_res), ]
    stopifnot(all.equal(drug_annot$title, names(query_res)))

    tibble::add_column(drug_annot,
                       Rank = NA,
                       Correlation = query_res,
                       .before=0)
  })


  # sort by absolute if either is genetic or is CMAP/L1000 pert
  is_genetic <- reactive(drug_study() == 'L1000 Genetic')
  sort_abs <- reactive(is_genetic() | is_pert())

  # subset to selected cells, summarize by compound, and add html
  query_table_summarised <- reactive({
    cols <- if(is_genetic()) gene_cols else drug_cols

    query_table_annot() %>%
      limit_cells(cells()) %>%
      summarize_compound(is_genetic = sort_abs()) %>%
      add_table_html() %>%
      select(cols, everything())
  })



  # build up to final table
  # ---

  # subset by min signatures
  query_table_nsig <- reactive(
    query_table_summarised() %>% dplyr::filter(n >= min_signatures())
  )

  # subset by clinical phase
  query_table_clin <- reactive({
    query_table_nsig() %>% {
      if (show_clinical() && 'Clinical Phase' %in% colnames(.))
        dplyr::filter(., !is.na(`Clinical Phase`))
      else .
    }
  })

  # final sorting/filtering
  query_table_final <- reactive({
    sort_by <- sort_by()
    sort_abs <- sort_abs()
    direction <- direction()
    drug_study <- drug_study()
    q <- query_table_clin()

    if (sort_by == 'avg_cor') q$Correlation <- gsub('simplot', 'simplot show-meanline', q$Correlation)

    # show largest absolute correlations first for genetic and pert queries
    # as both directions are informative
    if (sort_abs) {

      # filter none, opposing, or similar signatures based on direction toggle
      if (sort_by == 'avg_cor') {
        is.sim <- q$avg_cor > 0

      } else {
        mm <- q[, c('min_cor', 'max_cor')]
        mcol <- max.col(abs(mm), ties.method = 'last')
        is.sim <- mcol == 2 & mm[, 2] > 0
      }

      q <- switch(direction, 'both' = q, 'similar' = q[is.sim, ], 'opposing' = q[!is.sim, ]) %>%
        dplyr::mutate(min_cor = -pmax(abs(min_cor), abs(max_cor))) %>%
        dplyr::mutate(avg_cor = -abs(avg_cor))
    }

    # indicate total number of unique perts in title for rank
    rank_title <- switch(drug_study,
                         'CMAP02' = 'out of 1,309',
                         'L1000 Genetic' = 'out of 6,943',
                         'L1000 Drugs' = 'out of 19,360')

    # sort as desired then add rank
    q %>%
      dplyr::arrange(!!sym(sort_by)) %>%
      dplyr::select(-min_cor, -avg_cor, -max_cor, -n) %>%
      dplyr::mutate(Rank = paste0('<span class="rank-label label label-default" title="', rank_title, '">', 1:nrow(q), '</span>'))

  })

  # will update with proxy to prevent redraw
  dummy_table <- reactive({
    study <- drug_study()
    req(study, query_res())
    dummy_rendered(FALSE)
    cols <- if (study != 'L1000 Genetic') drug_cols else gene_cols

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
    content = function(con) {write.csv(query_table_dl(), con, row.names = FALSE)}
  )

}

#' Logic query/drug genes plotly
#' @export
#' @keywords internal
#' @importFrom magrittr "%>%"
drugsGenesPlotly <- function(input, output, session, data_dir, top_table, ambient, drug_study, pert_signature) {


  path_id <- reactive({
    study <- drug_study()
    if (study == 'CMAP02') return('Query genes - CMAP02')
    if (grepl('^L1000', study)) return('Query genes - L1000')
  })




  # the gene plot
  pl <- reactive({

    top_table <- top_table()
    ambient <- ambient()
    path_id <- path_id()

    if (is.null(path_id) | is.null(top_table)) return(NULL)

    pert_signature <- pert_signature()
    path_df <- get_path_df(top_table, path_id, pert_signature, ambient = ambient)

    dprimesPlotly(path_df)

  })

  observe(toggle('container', condition = isTruthy(pl())))

  output$plotly <- snapshotPreprocessOutput(
    plotly::renderPlotly({
      pl()
    }),
    function(value) { 'genes_plotly' }
  )
}

#' Generate plotly of dprimes values for Drugs and Pathways tab
#'
#' @param path_df result of \link{get_path_df}.
#'
#' @return
#' @export
#'
#' @keywords internal
dprimesPlotly <- function(path_df) {

  # so that still shows hover if no sd
  path_df$sd[is.na(path_df$sd)] <- 'NA'

  # 30 pixels width per gene in pathway
  ngenes <- length(unique(path_df$Gene))
  plot_height <- max(400, ngenes*25 + 125)


  (pl <- plotly::plot_ly(data = path_df,
                         y = ~Gene,
                         x = ~Dprime,
                         text = ~Gene,
                         customdata = apply(path_df, 1, as.list),
                         type = 'scatter',
                         mode = 'markers',
                         height = plot_height,
                         marker = list(size = 5, color = path_df$color),
                         error_x = ~list(array = sd, color = '#000000', thickness = 0.5, width = 0),
                         hoverlabel = list(bgcolor = '#000000', align = 'left'),
                         hovertemplate = paste0(
                           '<span style="color: crimson; font-weight: bold; text-align: left;">Gene</span>: %{text}<br>',
                           '<span style="color: crimson; font-weight: bold; text-align: left;">Description</span>: %{customdata.description}<br>',
                           '<span style="color: crimson; font-weight: bold; text-align: left;">Dprime</span>: %{x:.2f}<br>',
                           '<span style="color: crimson; font-weight: bold; text-align: left;">SD</span>: %{customdata.sd:.2f}',
                           '<extra></extra>')
  ) %>%
      plotly::config(displayModeBar = FALSE) %>%
      plotly::layout(hoverdistance = -1,
                     hovermode = 'y',
                     margin = list(t = 65, r = 20, l = 0, pad = 10),
                     title = list(text = 'Standardized Effect Size for Query Genes', y = 1, x = 0),
                     xaxis = list(fixedrange = TRUE, rangemode = "tozero", side = 'top', title = '', tickfont = list(size = 12)),
                     yaxis = list(fixedrange = TRUE,
                                  title = '',
                                  range = c(ngenes, -1),
                                  tickmode = 'array',
                                  tickvals = 0:ngenes,
                                  ticktext = ~Link,
                                  tickfont = list(size = 12)),
                     autosize = TRUE))


  # add arrow to show drug effect
  if ('dprime_sum' %in% colnames(path_df))
    pl <- pl %>%
    plotly::add_annotations(x = ~dprime_sum,
                            y = ~Gene,
                            xref = "x", yref = "y",
                            axref = "x", ayref = "y",
                            text = "",
                            showarrow = TRUE,
                            arrowcolor = ~arrow_color,
                            arrowwidth = 1,
                            ay = ~Gene,
                            ax = ~Dprime)

  return(pl)
}



#' Logic for selected dataset/analysis in Drugs and Pathways tabs
#' @export
#' @keywords internal
selectedAnal <- function(input, output, session, data_dir, choices, pert_query_dir = NULL) {

  # update dataset/analysis choice
  observe({
    choices <- choices()
    req(choices)
    updateSelectizeInput(session, 'query', choices = choices, server = TRUE, options = list(render = I('{item: querySignatureItem}')))
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
    if (!isTruthy(row_num)) return(list(type = ''))
    req(row_num, choices)

    choices[row_num, ]
  })

  # the type of dataset/analysis
  is_bulk <- reactive(sel()$type == 'Bulk Data')
  is_sc <- reactive(sel()$type == 'Single Cell')
  is_custom <- reactive(sel()$type == 'Custom')
  is_pert <- reactive(sel()$type == 'CMAP02/L1000 Perturbations')


  observe({
    shinyjs::toggle('sc_clusters_container', condition = is_sc())
    shinyjs::toggle('bulk_groups_container', condition = is_bulk())
  })


  # show/hide custom signature inputs in Drugs tab
  show_custom <- reactive(input$show_custom %% 2 != 0)

  observe({
    toggleClass(id = "show_custom", 'btn-primary', condition = show_custom())
  })

  sel_name <- reactive(sel()$label)

  dataset_dir <- reactive({
    if (!isTruthy(input$query)) return(NULL)
    file.path(data_dir, sel()$dataset_dir)
  })


  # Bulk analysis
  # ---
  eset  <- reactive(readRDS(file.path(dataset_dir(), 'eset.rds')))
  pdata <- reactive(readRDS(file.path(dataset_dir(), 'pdata_explore.rds')))
  svobj <- reactive(readRDS(file.path(dataset_dir(), 'svobj.rds')))

  numsv_path <- reactive(file.path(dataset_dir(), 'numsv.rds'))
  numsv <- reactive({
    numsv_path <- numsv_path()
    if (!file.exists(numsv_path)) saveRDS(0, numsv_path)
    readRDS(numsv_path)
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
                         is_bulk = is_bulk)



  # Single cell analysis
  # ---
  scSampleComparison <- callModule(scSampleComparison, 'sc',
                                   is_sc = is_sc,
                                   dataset_dir = dataset_dir)


  # results for selected type (bulk, sc, custom, pert)
  drug_queries <- reactive({
    sel_name <- sel_name()
    if (is_sc()) {
      drug_queries <- scSampleComparison$drug_queries()

    } else if (is_bulk()) {
      drug_queries <- bulkAnal$drug_queries()

    } else if (is_custom()) {
      drug_paths <- get_drug_paths(sel$dataset_dir, sel_name)
      drug_queries <- lapply(drug_paths, readRDS)

    } else if (is_pert()) {
      drug_paths <- get_drug_paths(pert_query_dir, fs::path_sanitize(sel_name))
      sapply(drug_paths, dl_pert_result)
      drug_queries <- lapply(drug_paths, readRDS)

    } else {
      drug_queries <- NULL
    }

    return(drug_queries)
  })

  top_table <- reactive({
    if (is_sc()) {
      top_table <- scSampleComparison$top_table()

    } else if (is_bulk()) {
      top_table <- bulkAnal$top_table()

    } else {
      top_table <- NULL
    }
    return(top_table)
  })

  ambient <- reactive({
    if (is_sc()) {
      ambient <- scSampleComparison$ambient()

    } else if (is_bulk()) {
      ambient <- NULL
    }
  })

  path_res <- reactive({
    if (is_sc()) {
      path_res <- scSampleComparison$path_res()

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
      anal_name <- paste(sel_name, scSampleComparison$name(), sep = '_')

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
    ambient = ambient,
    path_res = path_res,
    show_custom = show_custom,
    drug_queries = drug_queries,
    is_bulk = is_bulk,
    is_sc = is_sc,
    is_custom = is_custom,
    is_pert = is_pert
  ))

}

