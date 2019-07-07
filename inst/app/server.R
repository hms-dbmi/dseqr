

server <- function(input, output, session) {

  # get arguments from calling function
  # defaults for server
  data_dir <- getShinyOption('data_dir', '/srv/shiny-server/drugseqr/scseq/sjia')

  # for testing don't seem to be able to pass arguments as options
  if (isTRUE(getOption('shiny.testmode'))) {

    # reset data for testing
    data_dir <- 'tests/data/test'
    static_dir <- 'tests/data/static'
    unlink(data_dir, recursive = TRUE)
    dir.create(data_dir)
    file.copy(list.files(static_dir, full.names = TRUE), data_dir, recursive = TRUE)
  }


  # single cell analysis and options
  scPage <- callModule(scPage, 'sc',
                       data_dir = data_dir)

}
