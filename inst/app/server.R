



pathPage <- function(input, output, session, new_anal) {
  form <- callModule(pathForm, 'form',
                     new_anal = new_anal)

  output$path_plot <- shiny::renderPlot({

    df <- data.frame(gene = c('CRX1', 'CRX2', 'CRX3'),
                     dprime = c(5, -4, 9.2),
                     sd = c(0.1, 1, 0.5)
    )

    ggplot2::ggplot(data = df, ggplot2::aes_string(x = 'gene', y = 'dprime')) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=dprime-sd, ymax=dprime+sd), width = 0.03, size = 0.15) +
      ggplot2::ylab("Dprime") +
      ggplot2::xlab("") +
      ggplot2::geom_hline(yintercept = 0, color = 'black', linetype = 2, size = .2) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust=0.5),
                     legend.title = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_line(size = 0.2),
                     panel.border = ggplot2::element_rect(size = 0.05), text = element_text(size = 14.5))



  })
}


pathForm <- function(input, output, session, new_anal) {

  # reload query choices if new analysis
  anals <- reactive({
    new_anal()
    load_bulk_anals(data_dir)
  })

  observe({
    anals <- anals()
    req(anals)
    updateSelectizeInput(session, 'anal', choices = rbind(rep(NA, 5), anals), server = TRUE)
  })

  anal <- reactive({
    row_num <- input$anal
    anals <- anals()
    req(row_num, anals)

    anals[row_num, ]
  })


  # paths to analysis and drug query results
  fpaths <- reactive({
    anal <- anal()
    dataset_dir <- file.path(data_dir, 'bulk', anal$dataset_dir)
    anal_name <- anal$anal_name

    list(
      diff_anal = file.path(dataset_dir, paste0('diff_expr_symbol_', anal_name, '.rds')),
      diff_path = file.path(dataset_dir, paste0('diff_path_', anal_name, '.rds'))
    )
  })

  # the results
  diffs <- reactive({
    fpaths <- fpaths()

    list(
      path = readRDS(fpaths$diff_path),
      anal = readRDS(fpaths$diff_anal)
    )
  })

  observe({
    diffs <- diffs()
    res <- diffs$path$res

    path_choices <- data.frame(
      name = res$Name,
      value = res$ID,
      label = res$Name,
      fdr = format.pval(res$Ppadog, eps = 0.001, digits = 2),
      stringsAsFactors = FALSE)

    updateSelectizeInput(session, 'path_name',
                         choices = path_choices,
                         options = list(render= I('{option: pathOptions}')),
                         server = TRUE)
  })



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
