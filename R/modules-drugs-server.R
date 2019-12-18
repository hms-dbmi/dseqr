#' Logic for Drugs page
#' @export
#' @keywords internal
drugsPage <- function(input, output, session, new_dataset, bulk_changed, data_dir, pert_query_dir, pert_signature_dir) {

  # the form area inputs/results
  # form <- callModule(drugsForm, 'form',
  #                    new_dataset = new_dataset,
  #                    bulk_changed = bulk_changed,
  #                    data_dir = data_dir,
  #                    pert_query_dir = pert_query_dir,
  #                    pert_signature_dir = pert_signature_dir)

  # callModule(drugsGenesPlotly, 'genes',
  #            data_dir = data_dir,
  #            is_sc = form$is_sc,
  #            anal = form$anal,
  #            pert_signature = form$pert_signature,
  #            sc_inputs = form$sc_inputs,
  #            show_genes = form$show_genes,
  #            drug_study = form$drug_study)


  # the output table
  # drugs_table <- callModule(drugsTable, 'table',
  #                           data_dir = data_dir,
  #                           query_res = form$query_res,
  #                           drug_study = form$drug_study,
  #                           anal = form$anal,
  #                           cells = form$cells,
  #                           sort_by = form$sort_by,
  #                           show_clinical = form$show_clinical,
  #                           min_signatures = form$min_signatures,
  #                           is_pert = form$is_pert,
  #                           is_sc = form$is_sc,
  #                           sc_inputs = form$sc_inputs,
  #                           direction = form$direction)
  #
  #
  # # download link for table
  # output$dl_drugs <- downloadHandler(
  #   filename = function() {
  #     drug_study <- form$drug_study()
  #     anal_name <- form$anal()$anal_name
  #     drug_study <- tolower(drug_study)
  #     drug_study <- gsub(' ', '_', drug_study)
  #     paste0(anal_name, '_', drug_study, '.csv')
  #   },
  #   content = function(con) {write.csv(drugs_table$query_table_dl(), con, row.names = FALSE)}
  # )

}

# Logic for form on drugs page
#' @export
#' @keywords internal
drugsForm <- function(input, output, session, data_dir, new_dataset, bulk_changed, pert_query_dir, pert_signature_dir) {

  cmap_res <- reactiveVal()
  l1000_drugs_res <- reactiveVal()
  l1000_genes_res <- reactiveVal()

  # TODO: trigger new_custom after running custom query
  new_custom <- reactiveVal()

  # dataset/analysis choices
  choices <- reactive({
    # reactive to new datasets, new custom query, or bulk change (e.g. number of SVs)
    new_dataset()
    new_custom()
    bulk_changed()
    scseq_datasets <- load_scseq_datasets(data_dir)
    bulk_datasets <- load_bulk_datasets(data_dir)
    custom_anals <- load_custom_anals(data_dir)
    pert_anals <- load_pert_anals()

    choices <- rbind(bulk_datasets, scseq_datasets, custom_anals, pert_anals)
    choices$value <- seq_len(nrow(choices))

    return(choices)
  })

  # the selected dataset/query signaure
  selectedAnal <- callModule(selectedAnal, 'anal',
                             choices = choices,
                             data_dir = data_dir,
                             pert_query_dir = pert_query_dir)


  # end of filenames for drug queries
  fname_end <- reactive({
    type <- selectedAnal$type()
    sel <- selectedAnal$sel()
    if (type$is_bulk) {
      fname_end <- paste0(selectedAnal$bulk_name(), '_', selectedAnal$numsv(), 'svs')

    } else if (type$is_sc) {
      #TODO
      fname_end <- selectedAnal$sc_name()

    } else {
      fname_end <- fs::path_sanitize(sel$label)

    }
    fname_end <- paste0(fname_end, '.rds')
    return(fname_end)
  })


  # path to query results
  res_paths <- reactive({
    dataset_dir <- selectedAnal$dataset_dir()
    fname_end <- fname_end()

    list(
      cmap = file.path(dataset_dir, paste0('cmap_res_', fname_end)),
      l1000_drugs = file.path(dataset_dir, paste0('l1000_drugs_res_', fname_end)),
      l1000_genes = file.path(dataset_dir, paste0('l1000_genes_res_', fname_end))
    )
  })

  observe({
    print(res_paths())
  })

  # sc_inputs <- scSampleComparison(input, output, session,
  #                                 data_dir = data_dir,
  #                                 anal = querySignature$anal,
  #                                 is_sc = is_sc,
  #                                 input_ids = input_ids,
  #                                 with_drugs = TRUE)
  #
  #
  #
  # # if currently selected analysis is custom then show genes
  #
  # custom_query <- callModule(customQueryForm, 'custom-query',
  #                            show_custom = querySignature$show_custom,
  #                            is_custom = is_custom,
  #                            anal = querySignature$anal,
  #                            new_anal = new_anal,
  #                            data_dir = data_dir)
  #
  #
  # # get saved cmap/l1000 query results
  # observe({
  #   disable('signature')
  #   is_sc <- querySignature$is_sc()
  #   is_pert <- querySignature$is_pert()
  #   is_custom <- querySignature$is_custom()
  #
  #   if (is_sc)  {
  #     res <- sc_inputs$results()
  #
  #   } else if (is_custom || is_pert) {
  #     res_paths <- querySignature$res_paths()
  #     res <- load_custom_results(res_paths, is_pert = is_pert)
  #
  #   } else {
  #     res <- bulk_inputs$results()
  #   }
  #
  #   cmap_res(res$cmap)
  #   l1000_drugs_res(res$l1000_drugs)
  #   l1000_genes_res(res$l1000_genes)
  #
  #   enable('signature')
  # })
  #
  #
  # drugStudy <- callModule(selectedDrugStudy, 'drug_study',
  #                         anal = querySignature$anal,
  #                         is_pert = is_pert)
  #
  #
  #
  # advancedOptions <- callModule(advancedOptions, 'advanced',
  #                               cmap_res = cmap_res,
  #                               l1000_res = l1000_res,
  #                               drug_study = drugStudy$drug_study,
  #                               show_advanced = drugStudy$show_advanced)
  #
  # pert_signature <- callModule(selectedPertSignature, 'genes',
  #                              data_dir = data_dir,
  #                              pert_signature_dir = pert_signature_dir,
  #                              show_genes = drugStudy$show_genes,
  #                              drug_study = drugStudy$drug_study,
  #                              is_sc = is_sc,
  #                              is_bulk = is_bulk,
  #                              sc_inputs = sc_inputs,
  #                              bulk_inputs = bulk_inputs,
  #                              anal = querySignature$anal)
  #
  #
  # query_res <- reactive({
  #   drug_study <- drugStudy$drug_study()
  #   cmap_res <- cmap_res()
  #   l1000_drugs_res <- l1000_drugs_res()
  #   l1000_genes_res <- l1000_genes_res()
  #
  #   if (drug_study == 'CMAP02') return(cmap_res)
  #   else if (drug_study == 'L1000 Drugs') return(l1000_drugs_res)
  #   else if (drug_study == 'L1000 Genetic') return(l1000_genes_res)
  #
  #   return(NULL)
  # })
  #
  #
  #
  #
  # return(list(
  #   query_res = reactiveVal(),
  #   drug_study = drugStudy$drug_study,
  #   anal = querySignature$anal,
  #   cells = advancedOptions$cells,
  #   sort_by = advancedOptions$sort_by,
  #   show_clinical = drugStudy$show_clinical,
  #   show_genes = drugStudy$show_genes,
  #   min_signatures = advancedOptions$min_signatures,
  #   is_pert = querySignature$is_pert,
  #   is_sc = querySignature$is_sc,
  #   sc_inputs = sc_inputs,
  #   direction = drugStudy$direction,
  #   pert_signature = pert_signature
  # ))


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
    sel_idx <- which(choices$value == sel)

    updateSelectizeInput(session, 'query', choices = choices, selected = sel_idx, server = TRUE)
  })

  # the selected dataset/analysis
  sel <- reactive({
    row_num <- input$query
    choices <- choices()
    req(row_num, choices)

    choices[row_num, ]
  })

  # the type of dataset/analysis
  type <- reactive({
    type <- sel()$type
    list(
      is_bulk = type == 'Bulk Data',
      is_sc = type == 'Single Cell',
      is_custom = type == 'Custom',
      is_pert = type == 'CMAP02/L1000 Perturbations'
    )
  })

  observe({
    shinyjs::toggle('sc_clusters_container', condition = type()$is_sc)
    shinyjs::toggle('bulk_groups_container', condition = type()$is_bulk)
  })


  # show/hide custom signature inputs in Drugs tab
  show_custom <- reactive(input$show_custom %% 2 != 0)

  observe({
    toggleClass(id = "show_custom", 'btn-primary', condition = show_custom())
  })

  # Bulk analysis
  # ---
  fastq_dir <- reactive({
    req(type()$is_bulk)
    file.path(data_dir, sel()$dataset_dir)
  })

  eset  <- reactive(readRDS(file.path(fastq_dir(), 'eset.rds')))
  pdata <- reactive(readRDS(file.path(fastq_dir(), 'pdata_explore.rds')))
  numsv <- reactive(readRDS(file.path(fastq_dir(), 'numsv.rds')))
  svobj <- reactive(readRDS(file.path(fastq_dir(), 'svobj.rds')))


  explore_eset <- exploreEset(eset = eset,
                              fastq_dir = fastq_dir,
                              explore_pdata = pdata,
                              numsv = numsv,
                              svobj = svobj)

  bulkAnal <- callModule(bulkAnal, 'drugs',
                         pdata = pdata,
                         eset = explore_eset,
                         svobj = svobj,
                         numsv = numsv,
                         fastq_dir = fastq_dir,
                         type = type)

  # Single cell analysis
  # ---
  sc_name <- reactive({
    req(type()$is_sc)
    sel()$value
  })
  scAnal <- callModule(scAnal, 'drugs',
                       data_dir,
                       anal_name = sc_name
  )

  # TODO: get scAnal

  # folder for selected dataset
  dataset_dir <- reactive({
    sel <- sel()
    type <- type()

    if (type$is_pert) {
      dataset_dir <- pert_query_dir

    } else {
      dataset_dir <- file.path(data_dir, sel$dataset_dir)
    }

    return(dataset_dir)
  })

  return(list(
    lm_fit = bulkAnal$lm_fit,
    is_lmfit = bulkAnal$is_lmfit,
    bulk_name = bulkAnal$name,
    sc_name = scAnal$name,
    numsv = numsv,
    dataset_dir = dataset_dir,
    sel = sel,
    type = type
  ))

}

scAnal <- function(input, output, session, data_dir, anal_name, with_drugs = FALSE, with_path = FALSE) {
  contrast_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))

  # differential expression, pathway analyses, and drug queries
  results <- reactiveVal()

  sc_dir <- reactive({
    file.path(data_dir, 'single-cell')
  })

  scseq <- reactive({
    scseq_path <- scseq_part_path(sc_dir(), anal_name(), 'scseq')
    readRDS(scseq_path)
  })

  # TODO: update if annotation change from Single Cell tab
  annot <- reactive({
    annot_path <- scseq_part_path(sc_dir(), anal_name(), 'annot')
    readRDS(annot_path)
  })


  cluster_choices <- reactive({
    get_cluster_choices(clusters = annot(), anal_name(), sc_dir(), sample_comparison = TRUE)
  })


  # update UI for contrast/cluster choices
  observeEvent(cluster_choices(), {
    updateSelectizeInput(session, 'selected_clusters',
                         choices = cluster_choices(),
                         options = contrast_options, server = TRUE)
  })

  # reset results if change selected clusters
  observeEvent(input$selected_clusters, {
    result(NULL)
  })


  observeEvent(input$run_comparison, {

    selected_clusters <- input$selected_clusters
    req(selected_clusters)

    toggleAll(input_ids)

    scseq <- scseq()
    sc_dir <- sc_dir()
    anal_name <- anal_name()


    res <- run_comparison(scseq,
                          selected_clusters = selected_clusters,
                          sc_dir = sc_dir,
                          anal_name = anal_name,
                          session = session,
                          with_path = with_path,
                          with_drugs = with_drugs)


    toggleAll(input_ids)
    results(res)
  })
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
selectedDrugStudy <- function(input, output, session, anal, is_pert) {

  drug_study <- reactive(input$study)

  # boolean for advanced options
  show_advanced <- reactive({
    input$advanced %% 2 != 0
  })

  show_genes <- reactive({
    input$show_genes %% 2 != 0
  })

  observe({
    req(anal())
    choices <- data.frame(study = c('CMAP02', 'L1000', 'L1000'),
                          subset = c('drugs', 'drugs', 'genetic'),
                          value = c('CMAP02', 'L1000 Drugs', 'L1000 Genetic'),
                          stringsAsFactors = FALSE)

    prev_selected <- isolate(input$study)
    if (prev_selected == '') prev_selected <- NULL

    updateSelectizeInput(session, 'study', choices = choices, selected = prev_selected,
                         options = list(render = I('{option: studyOption, item: studyItem}')), server = TRUE)
  })

  # toggle for clinical status
  show_clinical <- reactive({
    input$clinical %% 2 != 0
  })


  observe({
    toggleClass('advanced', 'btn-primary', condition = show_advanced())
    toggleClass('clinical', 'btn-primary', condition = show_clinical())
    toggleClass('show_genes', 'btn-primary', condition = show_genes())
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
  is_genetic <- reactive(drug_study() == 'L1000 Genetic')

  # sort by absolute if either is genetic or is CMAP/L1000 pert
  observe(toggle('direction-parent', condition = is_genetic() | is_pert()))
  observe(toggle('clinical-parent', condition = !is_genetic()))


  return(list(
    drug_study = drug_study,
    show_clinical = show_clinical,
    show_advanced = show_advanced,
    show_genes = show_genes,
    direction = direction
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

#' Logic for advanced options in drugsForm
#' @export
#' @keywords internal
selectedPertSignature <- function(input, output, session, data_dir, pert_signature_dir, drug_study, show_genes, sc_inputs, bulk_inputs, is_sc, is_bulk, anal) {


  # toggle showing genes plotly inputs
  observe({
    toggle('pert_container', condition = show_genes())
  })

  # file paths to pathway, analysis, and drug query results
  custom_fpaths <- reactive({
    anal <- anal()
    is_sc <- is_sc()
    is_bulk <- is_bulk()

    if (is_sc | is_bulk) {
      return(NULL)
    }

    dataset_dir <-  file.path(data_dir, anal$dataset_dir)
    anal_name <- anal$anal_name

    list(
      anal = file.path(dataset_dir, paste0('diff_expr_symbol_', anal_name, '.rds')),
      cmap = file.path(dataset_dir, paste0('cmap_res_', anal_name, '.rds')),
      l1000_drugs = file.path(dataset_dir, paste0('l1000_drugs_res_', anal_name, '.rds')),
      l1000_genes = file.path(dataset_dir, paste0('l1000_genes_res_', anal_name, '.rds'))
    )
  })

  # load analysis results
  diffs <- reactive({
    fpaths <- custom_fpaths()
    if (is_sc()) return (sc_inputs$results())
    if (is_bulk()) return (bulk_inputs$results())

    list(anal = readRDS(fpaths$anal))
  })

  # load drug/genetic query results
  queries <- reactive({
    fpaths <- custom_fpaths()

    if (is_sc()) {
      sc_res <- sc_inputs$res()
      req(sc_res)

      cmap <- sc_res$cmap
      l1000_genes <- sc_res$l1000_genes
      l1000_drugs <- sc_res$l1000_drugs

    } else if (is_bulk()) {
      bulk_res <- bulk_inputs$res()
      req(bulk_res)

      cmap <- bulk_res$cmap
      l1000_genes <- bulk_res$l1000_genes
      l1000_drugs <- bulk_res$l1000_drugs

    } else {
      res <- run_drugs_comparison(fpaths, session)
      cmap <- res$cmap
      l1000_genes <- res$l1000_genes
      l1000_drugs <- res$l1000_drugs
    }

    # order pert choices same as in drugs
    list(
      cmap = sort(cmap),
      l1000_drugs = sort(l1000_drugs),
      l1000_genes = l1000_genes[order(abs(l1000_genes), decreasing = TRUE)]
    )

  })

  pert_type <- reactive({
    drug_study <- drug_study()
    req(drug_study)
    if (drug_study == 'CMAP02') return('cmap')
    ifelse(drug_study == 'L1000 Drugs', 'l1000_drugs', 'l1000_genes')
  })

  pert <- callModule(drugsPert,
                     'pert',
                     queries = queries,
                     pert_type = pert_type,
                     pert_signature_dir = pert_signature_dir)


  return(pert$signature)


}


#' Logic for drug table
#' @export
#' @keywords internal
#' @importFrom magrittr "%>%"
drugsTable <- function(input, output, session, query_res, drug_study, cells, show_clinical, sort_by, min_signatures, is_pert, direction) {
  drug_cols <- c('Rank', 'Correlation', 'Compound', 'Clinical Phase', 'External Links', 'MOA', 'Target', 'Disease Area', 'Indication', 'Vendor', 'Catalog #', 'Vendor Name')
  gene_cols <- c('Rank', 'Correlation', 'Compound', 'External Links', 'Description')
  pert_options <- list(render = I('{option: pertOptions, item: pertItem}'))



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
    drug_annot <- drug_annot[drug_annot$title %in% names(query_res), ]
    stopifnot(all.equal(drug_annot$title, names(query_res)))

    tibble::add_column(drug_annot,
                       Rank = NA,
                       Correlation = query_res,
                       .before=0)
  })

  is_genetic <- reactive({
    drug_study() == 'L1000 Genetic'
  })

  # sort by absolute if either is genetic or is CMAP/L1000 pert
  sort_abs <- reactive({
    is_genetic() | is_pert()
  })

  # subset to selected cells, summarize by compound, and add html
  query_table_summarised <- reactive({
    query_table_full <- query_table_full()
    if (is.null(query_table_full)) return(NULL)

    cols <- if(is_genetic()) gene_cols else drug_cols

    query_table <- query_table_full %>%
      limit_cells(cells()) %>%
      summarize_compound(is_genetic = sort_abs()) %>%
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

    # show largest absolute correlations first for genetic and pert queries
    # as both directions are informative
    if (sort_abs()) {

      # filter none, opposing, or similar signatures based on direction toggle
      if (sort_by == 'avg_cor') {
        is.similar <- query_table$avg_cor > 0
      } else {
        mm <- query_table[, c('min_cor', 'max_cor')]
        mcol <- max.col(abs(mm), ties.method = 'last')
        is.similar <- mcol == 2 & mm[, 2] > 0
      }

      query_table <- switch(direction(),
                            'both' = query_table,
                            'similar' = query_table[is.similar, ],
                            'opposing' = query_table[!is.similar, ])

      query_table <- query_table %>%
        dplyr::mutate(min_cor = -pmax(abs(min_cor), abs(max_cor))) %>%
        dplyr::mutate(avg_cor = -abs(avg_cor))
    }

    # indicate total number of unique perts in title for rank
    rank_title <- switch(drug_study(),
                         'CMAP02' = 'out of 1,309',
                         'L1000 Genetic' = 'out of 6,943',
                         'L1000 Drugs' = 'out of 19,360' )

    # sort as desired then add rank
    query_table <- query_table %>%
      dplyr::arrange(!!sym(sort_by)) %>%
      dplyr::select(-min_cor, -avg_cor, -max_cor, -n) %>%
      dplyr::mutate(Rank = paste0('<span class="rank-label label label-default" title="', rank_title, '">', 1:nrow(query_table), '</span>'))

    return(query_table)
  })

  # will update with proxy to prevent redraw
  dummy_table <- reactive({
    study <- drug_study()
    dummy_rendered(FALSE)
    cols <- if (study == 'L1000 Genetic') gene_cols else drug_cols

    data.frame(matrix(ncol = length(cols), dimnames = list(NULL, cols)), check.names = FALSE)
  })

  sorted_query <- reactive({
    query_res <- query_res()
    if (sort_abs()) {
      query_res <- query_res[order(abs(query_res), decreasing = TRUE)]
    } else {
      query_res <- sort(query_res)
    }
    return(query_res)
  })

  # query table for downloading
  query_table_dl <- reactive({
    query_res <- sorted_query()
    data.frame(correlation = query_res, signature = names(query_res), row.names = NULL)
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

  return(list(
    query_table_dl = query_table_dl
  ))
}

#' Logic query/drug genes plotly
#' @export
#' @keywords internal
#' @importFrom magrittr "%>%"
drugsGenesPlotly <- function(input, output, session, data_dir, anal, is_sc, sc_inputs, drug_study, show_genes, pert_signature) {

  #  toggle  showing genes plotly
  shiny::observe({
    toggle('container', condition = show_genes(), anim = TRUE)
  })

  # load pathway and analysis results
  diffs <- reactive({
    anal <- anal()
    if (is_sc()) return (sc_inputs$results())

    dataset_dir <-  file.path(data_dir, anal$dataset_dir)
    anal_name <- anal$anal_name

    list(
      anal = readRDS(file.path(dataset_dir, paste0('diff_expr_symbol_', anal_name, '.rds')))
    )
  })

  path_id <- reactive({
    study <- drug_study()
    if (study == 'CMAP02') return('Query genes - CMAP02')
    if (grepl('^L1000', study)) return('Query genes - L1000')
  })


  # the gene plot
  pl <- reactive({

    diffs <- diffs()
    path_id <- path_id()
    anal <- diffs$anal

    req(path_id, anal)

    pert_signature <- pert_signature()
    path_df <- get_path_df(anal, path_id, pert_signature)

    # so that still shows hover if no sd
    path_df$sd[is.na(path_df$sd)] <- 'NA'

    # 30 pixels width per gene in pathway
    ngenes <- length(unique(path_df$Gene))
    plot_width <- max(400, ngenes*25 + 125)

    pl <- plotly::plot_ly(data = path_df,
                          y = ~Dprime,
                          x = ~Gene,
                          text = ~Gene,
                          customdata = apply(path_df, 1, as.list),
                          type = 'scatter',
                          mode = 'markers',
                          width = plot_width,
                          height = 550,
                          marker = list(size = 5, color = path_df$color),
                          error_y = ~list(array = sd, color = '#000000', thickness = 0.5, width = 0),
                          hoverlabel = list(bgcolor = '#000000', align = 'left'),
                          hovertemplate = paste0(
                            '<span style="color: crimson; font-weight: bold; text-align: left;">Gene</span>: %{text}<br>',
                            '<span style="color: crimson; font-weight: bold; text-align: left;">Description</span>: %{customdata.description}<br>',
                            '<span style="color: crimson; font-weight: bold; text-align: left;">Dprime</span>: %{y:.2f}<br>',
                            '<span style="color: crimson; font-weight: bold; text-align: left;">SD</span>: %{customdata.sd:.2f}',
                            '<extra></extra>')
    ) %>%
      plotly::config(displayModeBar = FALSE) %>%
      plotly::layout(hoverdistance = -1,
                     hovermode = 'x',
                     yaxis = list(fixedrange = TRUE, rangemode = "tozero"),
                     xaxis = list(fixedrange = TRUE,
                                  range = c(-2, ngenes + 1),
                                  tickmode = 'array',
                                  tickvals = 0:ngenes,
                                  ticktext = ~Link,
                                  tickangle = -45),
                     autosize = FALSE)


    # add arrow to show drug effect
    if ('dprime_sum' %in% colnames(path_df))
      pl <- pl %>%
      plotly::add_annotations(x = ~Gene,
                              y = ~dprime_sum,
                              xref = "x", yref = "y",
                              axref = "x", ayref = "y",
                              text = "",
                              showarrow = TRUE,
                              arrowcolor = ~arrow_color,
                              arrowwidth = 1,
                              ax = ~Gene,
                              ay = ~Dprime)

    return(pl)

  })

  output$plotly <- snapshotPreprocessOutput(
    plotly::renderPlotly({
      pl()
    }),
    function(value) { 'genes_plotly' }
  )
}

#' Logic for perturbation selection in Pathways tab
#' @export
#' @keywords internal
drugsPert <- function(input, output, session, pert_type, queries, pert_signature_dir) {
  pert_options <- list(render = I('{option: pertOptions, item: pertItem}'))




  sorted_query <- reactive({

    queries <- queries()
    pert_type <- pert_type()
    req(pert_type)

    query_res <- queries[[pert_type]]

    query_res <- data.frame(
      label = gsub('_-700-666.0_\\d+h$', '', names(query_res)),
      value = names(query_res),
      cor = format(round(query_res, digits = 3)), stringsAsFactors = FALSE
    )
    return(query_res)
  })


  # update pert signature choices
  observe({
    updateSelectizeInput(session,
                         'pert',
                         choices = rbind(rep(NA, 3), sorted_query()),
                         options = pert_options,
                         server = TRUE)
  })

  # load signature for pert
  pert_signature <- reactive({
    pert <- input$pert
    pert_type <- pert_type()
    if (pert == '' | is.null(pert)) return(NULL)
    load_pert_signature(pert, pert_type, pert_signature_dir)
  })


  return(list(
    name = reactive(input$pert),
    signature = pert_signature
  ))

}
