
#' Logic for Pathways tab
#' @export
#' @keywords internal
pathPage <- function(input, output, session, new_anal, data_dir) {
  form <- callModule(pathForm, 'form',
                     new_anal = new_anal,
                     data_dir)

  observe({
    toggle('pert-legend', condition = isTruthy(form$pert_signature()))
  })


  # the gene plot
  pl <- reactive({

    diffs <- form$diffs()
    path_id <- form$pathway()
    anal <- diffs$anal

    req(path_id, anal)

    pert_signature <- form$pert_signature()
    path_df <- get_path_df(anal, path_id, pert_signature)

    # 30 pixels width per gene in pathway
    ngenes <- length(unique(path_df$Gene))
    plot_width <- max(400, ngenes*25 + 125)

    plotly::plot_ly(data = path_df,
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
      plotly::layout(yaxis = list(fixedrange = TRUE, rangemode = "tozero"),
                     xaxis = list(fixedrange = TRUE,
                                  range = c(-2, ngenes + 1),
                                  tickmode = 'array',
                                  tickvals = 0:ngenes,
                                  ticktext = ~Link,
                                  tickangle = -45),
                     autosize = FALSE)

  })

  # inside observe to allow dynamic width
  output$path_plot <- plotly::renderPlotly({
    pl()
  })


}


#' Logic for form in Pathways tab
#' @export
#' @keywords internal
pathForm <- function(input, output, session, new_anal, data_dir) {


  # reload analysis choices if new analysis
  anals <- reactive({
    new_anal()
    scseq_anals <- load_scseq_anals(data_dir, with_type = TRUE)
    bulk_anals <- load_bulk_anals(data_dir, with_type = TRUE)
    bulk_anals <- bulk_anals[bulk_anals$dataset_name != 'raw', ]

    anals <- rbind(bulk_anals, scseq_anals)
    anals$value <- seq_len(nrow(anals))

    return(anals)
  })


  # update analysis choices
  observe({
    anals <- anals()
    req(anals)
    updateSelectizeInput(session, 'anal', choices = rbind(rep(NA, 5), anals), options = list(render = I('{item: querySignatureItem}')), server = TRUE)
  })

  # get directory/name info about analysis
  anal <- reactive({
    row_num <- input$anal
    anals <- anals()
    req(row_num, anals)

    anals[row_num, ]
  })

  # show/hide perturbation signature selection
  is.cmap2 <- reactive(input$pathway == 'Query genes - CMAP02')
  is.l1000 <- reactive(input$pathway == 'Query genes - L1000')

  observe({
    toggle('l1000_pert_container', condition = is.l1000())
    toggle('cmap2_pert_container', condition = is.cmap2())
  })


  # show/hide single cell stuff
  is_sc <- reactive({
    anal <- anal()
    anal$type == 'Single Cell'
  })

  observe({
    shinyjs::toggle('sc_clusters_container', condition = is_sc())
  })

  # get single-cell data

  # inputs/buttons that can access/disable
  input_ids <- c('anal', 'kegg', 'pathway', 'run_comparison', 'selected_clusters')
  sc_inputs <- scSampleComparison(input, output, session,
                                  data_dir = data_dir,
                                  anal = anal,
                                  is_sc = is_sc,
                                  input_ids = input_ids,
                                  with_path = TRUE)


  # file paths to pathway, analysis, and drug query results
  fpaths <- reactive({
    anal <- anal()
    is_sc <- is_sc()

    if (is_sc()) {
      return(NULL)
    }

    dataset_dir <-  file.path(data_dir, anal$dataset_dir)
    anal_name <- anal$anal_name

    # for compatability with previous versions
    diff_path <- file.path(dataset_dir, paste0('diff_path_kegg_', anal_name, '.rds'))
    diff_path_old <- file.path(dataset_dir, paste0('diff_path_', anal_name, '.rds'))

    if(file.exists(diff_path_old)) file.rename(diff_path_old, diff_path)

    list(
      diff_anal = file.path(dataset_dir, paste0('diff_expr_symbol_', anal_name, '.rds')),
      diff_path = file.path(dataset_dir, paste0('diff_path_kegg_', anal_name, '.rds')),
      cmap = file.path(dataset_dir, paste0('cmap_res_', anal_name, '.rds')),
      l1000_drugs = file.path(dataset_dir, paste0('l1000_drugs_res_', anal_name, '.rds')),
      l1000_genes = file.path(dataset_dir, paste0('l1000_genes_res_', anal_name, '.rds'))
    )
  })

  # load pathway and analysis results
  diffs <- reactive({
    fpaths <- fpaths()
    if (is_sc()) return (sc_inputs$results())

    list(
      path = readRDS(fpaths$diff_path),
      anal = readRDS(fpaths$diff_anal)
    )
  })

  # load drug/genetic query results
  queries <- reactive({
    fpaths <- fpaths()

    # TODO: implement for single cell
    if (is_sc()) return(NULL)

    # order pert choices same as in drugs
    l1000_genes <- readRDS(fpaths$l1000_genes)
    l1000_genes <- l1000_genes[order(abs(l1000_genes), decreasing = TRUE)]

    list(
      cmap = sort(readRDS(fpaths$cmap)),
      l1000_drugs = sort(readRDS(fpaths$l1000_drugs)),
      l1000_genes = l1000_genes
    )

  })

  cmap2Pert <- callModule(pathPert, 'cmap2', type = 'CMAP02', queries = queries)
  l1000Pert <- callModule(pathPert, 'l1000', type = 'L1000', queries = queries)


  # get gex signature from drug/genetic query selection
  pert_signature <- reactive({
    if (!is.l1000() & !is.cmap2()) return(NULL)

    if (is.cmap2()) res <- cmap2Pert
    if (is.l1000()) res <- l1000Pert

    return(res$signature())
  })

  path_choices <- reactive({
    diffs <- diffs()
    res <- diffs$path$res

    if (is.null(res)) return(NULL)

    # for showing top up/down regulated
    all_choices <-  data.frame(
      name = c('used for CMAP02 drug queries', 'used for L1000 drug/genetic queries'),
      value = c('Query genes - CMAP02', 'Query genes - L1000'),
      label = c('Query genes - CMAP02', 'Query genes - L1000'),
      fdr = c(NA, NA),
      stringsAsFactors = FALSE
    )

    path_choices <- data.frame(
      name = res$Name,
      value = res$ID,
      label = res$Name,
      fdr = format.pval(res$FDRpadog, eps = 0.001, digits = 2),
      stringsAsFactors = FALSE)

    path_choices <- rbind(all_choices, path_choices)

    return(path_choices)
  })

  # update pathway dropdown
  observe({
    updateSelectizeInput(session, 'pathway',
                         choices = path_choices(),
                         options = list(render= I('{option: pathOptions, item: pathItem}')),
                         server = TRUE)
  })


  # open KEGG when click button
  observeEvent(input$kegg, {
    path_id <- input$pathway
    req(path_id)
    kegg_link <- paste0('https://www.genome.jp/kegg-bin/show_pathway?map', path_id)
    runjs(paste0("window.open('", kegg_link, "')"))
  })

  observe({
    toggleState('kegg', condition = input$pathway != 'all')
  })

  return(list(
    diffs = diffs,
    pathway = reactive(input$pathway),
    pert_signature = pert_signature
  ))
}

#' Logic for perturbation selection in Pathways tab
#' @export
#' @keywords internal
pathPert <- function(input, output, session, type, queries) {
  pert_options <- list(render = I('{option: pertOptions, item: pertItem}'))

  # update toggle for l1000 drugs/genes
  l1000_type <- reactive({
    req(!is.null(input$l1000_type))
    ifelse(input$l1000_type %% 2 == 0, 'drugs', 'genes')
  })

  observe({
    icon_name <- ifelse(l1000_type() == 'drugs', 'pills', 'dna')
    updateActionButton(session, 'l1000_type', icon = icon(icon_name, 'fa-fw'))
  })

  pert_type <- reactive({
    if (type == 'CMAP02') return('cmap')

    l1000_type <- l1000_type()
    if (type == 'L1000') return(paste0('l1000_', l1000_type))
  })


  sorted_query <- reactive({

    queries <- queries()
    pert_type <- pert_type()
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
    load_pert_signature(pert, pert_type)
  })


  return(list(
    name = reactive(input$pert),
    signature = pert_signature
  ))

}


#' Logic for single cell clusters selector for pathForm
#' @export
#' @keywords internal
scSampleComparison <- function(input, output, session, data_dir, anal, is_sc, input_ids, with_path = FALSE, with_drugs = FALSE) {
  contrast_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))

  # differentialexpression, pathway analyses, and drug queries
  results <- reactiveVal()

  sc_dir <- reactive({
    req(is_sc())
    file.path(data_dir, 'single-cell')
  })

  scseq <- reactive({
    scseq_path <- scseq_part_path(sc_dir(), anal()$anal_name, 'scseq')
    readRDS(scseq_path)
  })

  # TODO: update if annotation change from Single Cell tab
  annot <- reactive({
    annot_path <- scseq_part_path(sc_dir(), anal()$anal_name, 'annot')
    readRDS(annot_path)
  })


  cluster_choices <- reactive({
    # value is original cluster number so that saved pathway analysis name
    # isn't affected by updates to cluster annotation
    scseq <- scseq()
    value <- levels(Seurat::Idents(scseq))
    get_cluster_choices(clusters = annot(), scseq = scseq, value = value, sample_comparison = TRUE)
  })




  # update UI for contrast/cluster choices
  observeEvent(cluster_choices(), {
    updateSelectizeInput(session, 'selected_clusters',
                         choices = cluster_choices(),
                         options = contrast_options, server = TRUE)
  })


  observeEvent(input$run_comparison, {

    selected_clusters <- input$selected_clusters
    req(selected_clusters)

    toggleAll(input_ids)

    scseq <- scseq()
    sc_dir <- sc_dir()
    anal_name <- anal()$anal_name


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


  return(list(
    results = results,
    clusters = reactive(input$selected_clusters)
  ))


}
