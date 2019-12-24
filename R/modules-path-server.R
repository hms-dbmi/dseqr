#' Logic for Pathways tab
#' @export
#' @keywords internal
pathPage <- function(input, output, session, new_dataset, data_dir) {

  form <- callModule(pathForm, 'form',
                     new_dataset = new_dataset,
                     data_dir = data_dir)




  # the gene plot
  pl <- reactive({

    top_table <- form$top_table()
    path_res <- form$path_res()
    path_id <- form$path_id()

    req(path_id, top_table, path_res)
    path_df <- get_path_df(top_table, path_id)

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

  output$path_plot <- snapshotPreprocessOutput(
    plotly::renderPlotly({
      pl()
    }),
    function(value) { 'path_plotly' }
  )


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
    res <- topGO(res, sort = path_sort(), ontology = 'BP', number = Inf)

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
    top_table = selectedAnal$top_table,
    path_res = selectedAnal$path_res,
    path_id = reactive(input$pathway)
  ))
}
