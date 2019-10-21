server <- function(input, output, session) {

  # get arguments from calling function
  # defaults for server

  # base directory contains data_dir folder
  data_dir <- getShinyOption('data_dir', 'data_dir')

  # path where drug gene queries will be stored
  pert_query_dir <- getShinyOption(
    'pert_query_dir',
    '/srv/shiny-server/pert_query_dir')

  # path where kallisto index is downloaded and stored
  indices_dir <- getShinyOption(
    'indices_dir',
    '/srv/shiny-server/indices')

  if (!dir.exists(pert_query_dir)) dir.create(pert_query_dir)

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
                       data_dir = data_dir,
                       indices_dir = indices_dir)


  drugsPage <- callModule(drugsPage, 'drug',
                          new_anal = dsPage$new_anal,
                          data_dir = data_dir,
                          pert_query_dir = pert_query_dir)

  pathPage <- callModule(pathPage, 'pathways',
                         new_anal = dsPage$new_anal,
                         data_dir = data_dir)



}
