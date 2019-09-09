
#' Logic for Pathways tab
#' @export
#' @keywords internal
pathPage <- function(input, output, session, new_anal, data_dir) {
  form <- callModule(pathForm, 'form',
                     new_anal = new_anal,
                     data_dir)

  observe({
    toggle('l1000-label', condition = isTruthy(form$pathway()))
  })



  # the gene plot
  pl <- reactive({

    diffs <- form$diffs()
    path_id <- form$pathway()
    path_genes <- form$custom_path_genes()
    anal <- diffs$anal

    req(path_id, anal)

    path_df <- get_path_df(anal, path_id, path_genes)

    # 30 pixels width per gene in pathway
    plot_width <- max(400, nrow(path_df)*25 + 125)

    pt.color <- ifelse(path_df$Gene %in% genes$common, 'red', 'black')

    plotly::plot_ly(data = path_df,
                    y = ~Dprime,
                    x = ~Gene,
                    text = ~Gene,
                    customdata = apply(path_df, 1, as.list),
                    type = 'scatter',
                    mode = 'markers',
                    width = plot_width,
                    height = 550,
                    marker = list(size = 5, color = pt.color),
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
                                  range = c(-2, nrow(path_df) + 1),
                                  tickmode = 'array',
                                  tickvals = 0:nrow(path_df),
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


  # file paths to pathway and analysis results
  fpaths <- reactive({
    anal <- anal()
    is_sc <- is_sc()

    if (is_sc()) {
      return(NULL)
    }

    dataset_dir <-  file.path(data_dir, 'bulk', anal$dataset_dir)
    anal_name <- anal$anal_name

    list(
      diff_anal = file.path(dataset_dir, paste0('diff_expr_symbol_', anal_name, '.rds')),
      diff_path = file.path(dataset_dir, paste0('diff_path_', anal_name, '.rds'))
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

  path_choices <- reactive({
    diffs <- diffs()
    res <- diffs$path$res

    if (is.null(res)) return(NULL)

    # for showing top up/down regulated
    all_choices <-  data.frame(
      name = 'used for drug/genetic queries',
      value = 'Drug and genetic query genes',
      label = 'Drug and genetic query genes',
      fdr = NA,
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

  # show custom pathway inputs
  show_custom <- reactive(input$show_custom %% 2 != 0)

  observe({
    shinyjs::toggle('custom_path_container', anim = TRUE, condition = show_custom())
    toggleClass(id = "show_custom", 'btn-primary', condition = show_custom())
  })

  # choices for custom pathway genes
  gene_choices <- reactive({
    diffs <- diffs()
    req(diffs())
    gene <- row.names(diffs$anal$top_table)
    description <- tx2gene$description[match(gene, tx2gene$gene_name)]

    data.frame(gene, description, value = gene, label = gene, stringsAsFactors = FALSE)
  })

  observe({
    updateSelectizeInput(session, 'custom_path_genes', choices = gene_choices(), options = list(render = I('{option: pathGene, item: pathGene}')), server = TRUE)
  })

  return(list(
    diffs = diffs,
    show_custom = show_custom,
    custom_path_genes = reactive(input$custom_path_genes),
    pathway = reactive(input$pathway)
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
