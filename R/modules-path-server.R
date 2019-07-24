
#' Logic for Pathways tab
#' @export
#' @keywords internal
pathPage <- function(input, output, session, new_anal, data_dir) {
  form <- callModule(pathForm, 'form',
                     new_anal = new_anal,
                     data_dir)



  # the gene plot
  pl <- reactive({

    diffs <- form$diffs()
    path_id <- form$pathway()
    show_up <- form$show_up()
    anal <- diffs$anal

    req(path_id, anal)

    if (path_id == 'all') {
      path_df <- get_all_df(anal, show_up)
    } else {
      path_df <- get_path_df(path_id, anal)
    }

    # 30 pixels width per gene in pathway
    plot_width <- max(400, nrow(path_df)*25 + 125)


    plotly::plot_ly(data = path_df,
                    y = ~Dprime,
                    x = ~Gene,
                    text = ~Gene,
                    customdata = ~sd,
                    type = 'scatter',
                    mode = 'markers',
                    width = plot_width,
                    height = 550,
                    marker = list(size = 5, color = '#000000'),
                    error_y = ~list(array = sd, color = '#000000', thickness = 0.5, width = 0),
                    hovertemplate = paste0(
                      '<span style="color: crimson; font-weight: bold;">Gene</span>: %{text}<br>',
                      '<span style="color: crimson; font-weight: bold;">Dprime</span>: %{y:.2f}<br>',
                      '<span style="color: crimson; font-weight: bold;">SD</span>: %{customdata:.2f}',
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
    bulk_anals <- load_bulk_anals(data_dir, with_type = TRUE)
    scseq_anals <- load_scseq_anals(data_dir, with_type = TRUE)

    anals <- rbind(bulk_anals, scseq_anals)
    anals$value <- seq_len(nrow(anals))
    return(anals)
  })


  # update analysis choices
  observe({
    anals <- anals()
    req(anals)
    updateSelectizeInput(session, 'anal', choices = rbind(rep(NA, 5), anals), server = TRUE)
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
  sc_inputs <- scPathClusters(input, output, session,
                              data_dir = data_dir,
                              anal = anal,
                              is_sc = is_sc)



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
    if (is_sc()) return (sc_inputs$path_diffs())


    list(
      path = readRDS(fpaths$diff_path),
      anal = readRDS(fpaths$diff_anal)
    )
  })

  show_up <- reactive({
    input$show_up %% 2 == 0
  })

  observe({
    name <- ifelse(show_up(),  'chevron-up', 'chevron-down')
    updateActionButton(session, 'show_up', icon = icon(name, 'fa-fw'))
  })

  path_directions <- reactive({
    diffs <- diffs()
    res <- diffs$path$res

    if (is.null(res)) return(NULL)

    path_directions <- get_path_directions(diffs$anal$top_table)
    path_directions[res$ID, ]
  })

  path_choices <- reactive({
    diffs <- diffs()
    res <- diffs$path$res
    directions <- path_directions()

    if (is.null(res)) return(NULL)

    # for showing top up/down regulated
    all_choices <-  data.frame(
      name = 'all',
      value = 'all',
      label = 'all',
      direction_label = c('Mostly Up', 'Mostly Down'),
      is.up = c(TRUE, FALSE),
      fdr = c(NA, NA),
      stringsAsFactors = FALSE
    )

    path_choices <- data.frame(
      name = res$Name,
      value = res$ID,
      label = res$Name,
      direction_label = directions$label,
      is.up = directions$is.up,
      fdr = format.pval(res$Ppadog, eps = 0.001, digits = 2),
      stringsAsFactors = FALSE)

    path_choices <- rbind(all_choices, path_choices)

    filter <- if (show_up()) path_choices$is.up else !path_choices$is.up

    return(path_choices[filter, ])
  })

  # update pathway dropdown
  observe({
    updateSelectizeInput(session, 'pathway',
                         choices = path_choices(),
                         options = list(render= I('{option: pathOptions}')),
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
    show_up = show_up,
    pathway = reactive(input$pathway)
  ))
}

#' Logic for single cell clusters selector for pathForm
#' @export
#' @keywords internal
scPathClusters <- function(input, output, session, data_dir, anal, is_sc) {
  # inputs/buttons that can access/disable
  input_ids <- c('anal', 'kegg', 'pathway', 'run_comparison', 'selected_clusters', 'show_up')

  contrast_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))

  # differential and pathway analyses
  path_diffs <- reactiveVal()

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
    get_cluster_choices(clusters = annot(), scseq = scseq, value = value)
  })




  # update UI for contrast/cluster choices
  observeEvent(cluster_choices(), {
    updateSelectizeInput(session, 'selected_clusters',
                         choices = cluster_choices(),
                         options = contrast_options, server = TRUE)
  })


  observeEvent(input$run_comparison, {
    toggleAll(input_ids)


    # run differential expression analysis ----
    # set idents to ctrl and test
    scseq <- scseq()
    clusters <- input$selected_clusters
    prev_anal <- diff_expr_scseq(scseq, clusters = clusters)
    prev_anal <- ambient.omit(scseq, anal = prev_anal)

    # run pathway analysis
    panal <- diff_path_scseq(scseq,
                             prev_anal = prev_anal,
                             data_dir = sc_dir(),
                             anal_name = anal()$anal_name,
                             clusters = clusters)

    toggleAll(input_ids)
    path_diffs(list(anal = prev_anal, path = panal))
  })


  return(list(
    path_diffs = path_diffs,
    clusters = reactive(input$selected_clusters)
  ))


}
