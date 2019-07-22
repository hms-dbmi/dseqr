

get_path_df <- function(path_id, anal) {
  # add dprimes and vardprime values
  anal <- add_es(anal)
  top_table <- anal$top_table

  path_enids <- gslist[[path_id]]
  path_genes <- names(path_enids)

  # subset top table to genes in the pathway
  top_table <- top_table[row.names(top_table) %in% path_genes, ]

  path_df <- data.frame(
    Gene = row.names(top_table),
    Dprime = top_table$dprime,
    Vardprime = top_table$vardprime,
    Link = paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", row.names(top_table), "'>", row.names(top_table), "</a>"), stringsAsFactors = FALSE
  )

  path_df <- path_df %>%
    arrange(desc(abs(Dprime))) %>%
    mutate(Gene = factor(Gene, levels = Gene))

  return(path_df)
}

pathPage <- function(input, output, session, new_anal) {
  form <- callModule(pathForm, 'form',
                     new_anal = new_anal)



  # the gene plot
  pl <- reactive({

    diffs <- form$diffs()
    path_id <- form$pathway()
    anal <- diffs$anal

    req(path_id, anal)

    path_df <- get_path_df(path_id, anal)
    # 30 pixels width per gene in pathway
    plot_width <- nrow(path_df)*25


    plotly::plot_ly(data = path_df,
                    y = ~Dprime,
                    x = ~Gene,
                    text = ~Gene,
                    type = 'scatter',
                    mode = 'markers',
                    width = plot_width,
                    height = 550,
                    marker = list(size = 5, color = '#000000'),
                    error_y = ~list(array = Vardprime, color = '#000000', thickness = 0.5, width = 0),
                    hovertemplate = paste0(
                      '<b>Gene</b>: %{text}<br>',
                      '<b>Dprime</b>: %{y:.2f}',
                      '<extra></extra>')
                    ) %>%
      plotly::config(displayModeBar = FALSE) %>%
      plotly::layout(yaxis = list(fixedrange = TRUE),
                     xaxis = list(fixedrange = TRUE,
                                  range = c(-2, nrow(path_df) + 1),
                                  tickmode = 'array',
                                  tickvals = 0:nrow(path_df),
                                  ticktext = ~Link,
                                  tickangle = -45),
                     autosize = FALSE)

  })

  # inside observe to allow dynamic width
  output$path_plot <- renderPlotly({
    pl()
  })


}


pathForm <- function(input, output, session, new_anal) {

  # reload analysis choices if new analysis
  anals <- reactive({
    new_anal()
    load_bulk_anals(data_dir)
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


  # file paths to pathway and analysis results
  fpaths <- reactive({
    anal <- anal()
    dataset_dir <- file.path(data_dir, 'bulk', anal$dataset_dir)
    anal_name <- anal$anal_name

    list(
      diff_anal = file.path(dataset_dir, paste0('diff_expr_symbol_', anal_name, '.rds')),
      diff_path = file.path(dataset_dir, paste0('diff_path_', anal_name, '.rds'))
    )
  })

  # load pathway and analysis results
  diffs <- reactive({
    fpaths <- fpaths()

    list(
      path = readRDS(fpaths$diff_path),
      anal = readRDS(fpaths$diff_anal)
    )
  })

  # update pathway dropdown
  observe({
    diffs <- diffs()
    res <- diffs$path$res

    path_choices <- data.frame(
      name = res$Name,
      value = res$ID,
      label = res$Name,
      fdr = format.pval(res$Ppadog, eps = 0.001, digits = 2),
      stringsAsFactors = FALSE)

    updateSelectizeInput(session, 'pathway',
                         choices = path_choices,
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

  return(list(
    diffs = diffs,
    pathway = reactive(input$pathway)
  ))



}

get_path_choices <- function(res) {
  data.frame(
    name = res$Name,
    value = res$ID,
    fdr = res$Ppadog
  )
}

server <- function(input, output, session) {

  # get arguments from calling function
  # defaults for server
  # base directory contains data_dir folder
  data_dir <- getShinyOption('data_dir', '/srv/shiny-server/drugseqr/data_dir')

  # for testing don't seem to be able to pass arguments as options
  if (isTRUE(getOption('shiny.testmode'))) {

    # reset data for testing
    data_dir <- 'tests/data/test'
    static_dir <- 'tests/data/static'
    unlink(data_dir, recursive = TRUE)
    dir.create(data_dir)
    file.copy(list.files(static_dir, full.names = TRUE), data_dir, recursive = TRUE)
  }

  sc_dir <- file.path(data_dir, 'single-cell')
  bulk_dir <- file.path(data_dir, 'bulk')

  dir.create(sc_dir, showWarnings = FALSE)
  dir.create(bulk_dir, showWarnings = FALSE)


  # single cell analysis and options
  scPage <- callModule(scPage, 'sc',
                       sc_dir = sc_dir)

  dsPage <- callModule(dsPage, 'datasets',
                       data_dir = data_dir)


  drugsPage <- callModule(drugsPage, 'drug',
                          new_anal = dsPage$new_anal,
                          data_dir = data_dir)

  pathPage <- callModule(pathPage, 'pathways',
                         new_anal = dsPage$new_anal)



}
