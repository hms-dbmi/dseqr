#' Logic for Pathways tab
#' @export
#' @keywords internal
pathPage <- function(input, output, session, new_dataset, data_dir) {

  form <- callModule(pathForm, 'form',
                     new_dataset = new_dataset,
                     data_dir = data_dir)


  # the gene plot
  gene_pl <- reactive({

    top_table <- form$top_table()
    path_res <- form$path_res()
    path_id <- form$path_id()

    req(path_id, top_table, path_res)
    path_df <- get_path_df(top_table, path_id)

    dprimesPlotly(path_df)
  })

  output$path_plot <- snapshotPreprocessOutput(
    plotly::renderPlotly({
      gene_pl()
    }),
    function(value) { 'path_plotly' }
  )

  # heatmap plot
  heat_pl <- reactive({
    path_id <- form$path_id()
    req(path_id)

    eset <- form$bulk_eset()
    top_table <- form$top_table()
    contrast_groups <- form$contrast_groups()
    get_pathway_heatmap(eset, top_table, contrast_groups, path_id)
  })

  heat_width <- reactive({
    eset <- form$bulk_eset()
    pdata <- Biobase::pData(eset)
    contrast_groups <- form$contrast_groups()

    (sum(pdata$group %in% contrast_groups)*20) + 140
  })

  observe({

    output$heatmap <- renderPlot({
      req(form$is_bulk())
      heat_pl()
    }, width = heat_width())

  })


}



#' Get heatmap of genes in a GO pathways
#'
#' @param eset ExpressionSet
#' @param top_table limma topTable
#' @param contrast_groups groups to include in heatmap
#' @param path_id GO pathway id string
#'
#' @return \code{pheatmap} object.
#' @export
get_pathway_heatmap <- function(eset, top_table, contrast_groups, path_id) {

  # subset to pathway genes
  adj <- Biobase::assayDataElement(eset, 'adjusted')
  path_genes <- names(gslist.go[[path_id]])

  # subset to max 1000 top genes
  sig.genes <- head(row.names(top_table[path_genes, ]), 1000)
  adj <- adj[row.names(adj) %in% sig.genes, ]

  # get group colors and levels
  pdata <- Biobase::pData(eset)
  group_levels <- get_group_levels(pdata)
  group_colors <- get_group_colors(group_levels)

  # subset eset to selected groups
  in.sel <- pdata$group %in% contrast_groups
  adj <- adj[, in.sel]
  if ((nrow(adj) %% 2) == 1) adj <- rbind(adj, rnorm(ncol(adj)))

  # cluster rows on correlation
  dist.rows <- as.dist(1 - cor(t(adj)))
  annot <- data.frame(Group = pdata$group[in.sel], row.names = colnames(adj))
  annotation_colors <- list(Group = group_colors)
  names(annotation_colors$Group) <- group_levels
  annotation_colors$Group <- annotation_colors$Group[contrast_groups]

  pheatmap::pheatmap(adj,
                     color = gplots::greenred(10),
                     clustering_distance_rows = dist.rows,
                     show_colnames = TRUE,
                     show_rownames = FALSE,
                     border_color = NA,
                     scale = 'row',
                     treeheight_row = 0,
                     annotation_colors = annotation_colors,
                     annotation_names_col = FALSE,
                     annotation_col = annot,na_col = 'black',
                     cellwidth=18, fontsize = 12)


}


#' Logic for form in Pathways tab
#' @export
#' @keywords internal
pathForm <- function(input, output, session, new_dataset, data_dir) {

  # dataset/analysis choices
  choices <- reactive({
    # reactive to new datasets, new custom query, or bulk change (e.g. number of SVs)
    new_dataset()
    scseq_datasets <- load_scseq_datasets(data_dir)
    bulk_datasets <- load_bulk_datasets(data_dir)

    choices <- rbind(bulk_datasets, scseq_datasets)
    choices$value <- seq_len(nrow(choices))+1
    choices <- rbind(rep(NA, 5), choices)

    return(choices)
  })

  # the selected dataset/analysis results
  selectedAnal <- callModule(selectedAnal, 'path',
                             choices = choices,
                             data_dir = data_dir)


  # the pathway direction
  path_sort <- reactive({
    switch((input$direction %% 3) + 1,
           NULL,
           'up',
           'down')
  })

  sort_icon <- reactive({
    switch((input$direction %% 3) + 1,
           'arrows-alt-v',
           'chevron-up',
           'chevron-down')
  })

  observe(updateActionButton(session, 'direction', icon = icon(sort_icon(), 'fa-fw')))


  # the pathway choices
  path_choices <- reactive({
    res <- selectedAnal$path_res()

    if (is.null(res)) return(NULL)
    res <- limma::topGO(res, sort = path_sort(), ontology = 'BP', number = Inf)

    # limma authors recommend ignoring unadjusted p-values > 10^-5
    path_choices <- data.frame(
      name = res$Term,
      value = row.names(res),
      label = res$Term,
      dirLabel = ifelse(res$P.Up < res$P.Down, 'Up', 'Down'),
      ignore = pmin(res$P.Up, res$P.Down) > 10^-5,
      stringsAsFactors = FALSE)

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
  # observeEvent(input$kegg, {
  #   path_id <- input$pathway
  #   req(path_id)
  #   kegg_link <- paste0('https://www.genome.jp/kegg-bin/show_pathway?map', path_id)
  #   runjs(paste0("window.open('", kegg_link, "')"))
  # })

  return(list(
    bulk_eset = selectedAnal$bulk_eset,
    contrast_groups = selectedAnal$contrast_groups,
    top_table = selectedAnal$top_table,
    path_res = selectedAnal$path_res,
    is_bulk = selectedAnal$is_bulk,
    path_id = reactive(input$pathway)
  ))
}
