server <- function(input, output, session) {

  # get arguments from calling function
  # defaults for testing

  # base directory contains data_dir folder
  data_dir <- getShinyOption('data_dir', 'tests/data/test/example')

  # path where pert queries will be stored
  pert_query_dir <- getShinyOption('pert_query_dir', '/srv/drugseqr/pert_query_dir')

  # path where pert signatures will be stored
  pert_signature_dir <- getShinyOption('pert_signature_dir', '/srv/drugseqr/pert_signature_dir')

  # path where kallisto index is downloaded and stored
  indices_dir <- getShinyOption('indices_dir', '/srv/drugseqr/indices')

  if (!dir.exists(pert_query_dir)) dir.create(pert_query_dir)
  if (!dir.exists(pert_signature_dir)) dir.create(pert_signature_dir)

  # for testing don't seem to be able to pass arguments as options
  if (isTRUE(getOption('shiny.testmode'))) {
    # reset data for testing
    data_dir <- 'tests/data/test/example'
    static_dir <- 'tests/data/static/example'
    unlink(data_dir, recursive = TRUE)
    dir.create(data_dir)
    file.copy(list.files(static_dir, full.names = TRUE), data_dir, recursive = TRUE)
  }

  sc_dir <- file.path(data_dir, 'single-cell')
  bulk_dir <- file.path(data_dir, 'bulk')

  dir.create(sc_dir, showWarnings = FALSE)
  dir.create(bulk_dir, showWarnings = FALSE)



  bulkPage <- callModule(bulkPage, 'bulk',
                         data_dir = data_dir,
                         sc_dir = sc_dir,
                         bulk_dir = bulk_dir,
                         indices_dir = indices_dir)

  scPage <- callModule(scPage, 'sc',
                       sc_dir = sc_dir,
                       indices_dir = indices_dir)

  # TODO: get new_dataset from bulkPage and scPage
  new_dataset <- reactiveVal()
  bulk_changed <- reactiveVal()

  drugsPage <- callModule(drugsPage, 'drug',
                          data_dir = data_dir,
                          new_dataset = new_dataset,
                          bulk_changed = bulk_changed,
                          pert_query_dir = pert_query_dir,
                          pert_signature_dir = pert_signature_dir)


}
