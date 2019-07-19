

pathPage <- function(input, output, session) {

  output$path_plot <- plotly::renderPlotly({

    df <- data.frame(gene = c('CRX1', 'CRX2', 'CRX3'),
                     dprime = c(5, -4, 9.2)
                     )

    g <-  ggplot2::ggplot(data = df,
                          ggplot2::aes_string(x = 'gene', y = 'dprime')) +
      ggplot2::geom_point() +
      ggplot2::ylab("Dprime") +
      ggplot2::xlab("") +
      ggplot2::geom_hline(yintercept = 0, colour = '#999999') +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust=0.5),
                     legend.title = ggplot2::element_blank())


    pl <- plotly::plotly_build(g)

    pl
  })
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

  pathPage <- callModule(pathPage, 'pathways')



}
