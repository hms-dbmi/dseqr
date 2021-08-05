server <- function(input, output, session) {

    # get arguments from calling function
    # defaults for testing

    # base directory contains data_dir folder
    data_dir <- getShinyOption('data_dir', 'tests/data/test/example')

    # path where pert queries will be stored
    pert_query_dir <- getShinyOption('pert_query_dir', '/srv/dseqr/pert_query_dir')

    # path where gene set data is stored
    gs_dir <- getShinyOption('gs_dir', '/srv/dseqr/gs_dir')

    # path where pert signatures will be stored
    pert_signature_dir <- getShinyOption('pert_signature_dir', '/srv/dseqr/pert_signature_dir')

    # path where kallisto index is downloaded and stored
    indices_dir <- getShinyOption('indices_dir', '/srv/dseqr/indices')

    # path where save tx2genes
    tx2gene_dir <- getShinyOption('tx2gene_dir', '/srv/dseqr/tx2gene')

    is_example <- getShinyOption('is_example', FALSE)

    if (!dir.exists(pert_query_dir)) dir.create(pert_query_dir)
    if (!dir.exists(pert_signature_dir)) dir.create(pert_signature_dir)
    if (!dir.exists(indices_dir)) dir.create(indices_dir)
    if (!dir.exists(tx2gene_dir)) dir.create(tx2gene_dir)

    # for testing don't seem to be able to pass arguments as options
    if (isTRUE(getOption('shiny.testmode'))) {
        # reset data for testing
        data_dir <- 'tests/data/test/example'
        static_dir <- 'tests/data/static/example'
        unlink(data_dir, recursive = TRUE)
        dir.create(data_dir, recursive = TRUE)
        file.copy(list.files(static_dir, full.names = TRUE), data_dir, recursive = TRUE)
    }


    sc_dir <- file.path(data_dir, 'single-cell')
    bulk_dir <- file.path(data_dir, 'bulk')

    dir.create(sc_dir, showWarnings = FALSE)
    dir.create(bulk_dir, showWarnings = FALSE)

    # hide tour button for docs page
    observe(toggleClass('start_tour', 'invisible', condition = input$tabs == 'Docs'))


    # rintrojs
    observeEvent(input$start_tour, {
        if (input$tabs == 'Single Cell') {
            steps <- utils::read.csv('www/sc_intro.csv', stringsAsFactors = FALSE)
        } else if (input$tabs == 'Bulk Data') {
            steps <- utils::read.csv('www/bulk_intro.csv', stringsAsFactors = FALSE)
        } else if (input$tabs == 'Drugs') {
            steps <- utils::read.csv('www/drugs_intro.csv', stringsAsFactors = FALSE)
        } else {
            print(input$tabs)
            return(NULL)
        }

        rintrojs::introjs(session,
                          options = list(
                              showStepNumbers = 'false',
                              steps = steps))

    })

    feedback_counter <- reactiveVal(0)
    observeEvent(input$submit_feedback, {
        params <- rev(strsplit(data_dir, '/')[[1]])
        dataset <- params[1]
        user <- params[2]

        url <- readRDS(system.file('extdata/slack.rds', package = 'dseqr'))

        httr::POST(url = url,
                   httr::add_headers('Content-Type' = 'application/json'),
                   body = sprintf('{"text": "%s \n user: %s, app: %s"}', input$feedback, user, dataset))

        updateTextAreaInput(session, 'feedback', value = '', placeholder = 'Thank you!')
        shinyjs::delay(3000, updateTextAreaInput(session, 'feedback', placeholder = ''))
    })



    observe({
        toggle('add_dataset', condition = input$tabs == 'Single Cell')
        toggle('remove_dataset', condition = input$tabs == 'Single Cell')
    })


    is_mobile <- reactive(input$is_mobile)

    bulkPage <- callModule(bulkPage, 'bulk',
                           data_dir = data_dir,
                           sc_dir = sc_dir,
                           bulk_dir = bulk_dir,
                           gs_dir = gs_dir,
                           indices_dir = indices_dir)

    if (is_example) {
        add_sc <- reactiveVal()
        remove_sc <- reactiveVal()

    } else {
        add_sc <- reactive(input$add_dataset)
        remove_sc <- reactive(input$remove_dataset)
    }


    scPage <- callModule(scPage, 'sc',
                         sc_dir = sc_dir,
                         indices_dir = indices_dir,
                         tx2gene_dir = tx2gene_dir,
                         gs_dir = gs_dir,
                         is_mobile = is_mobile,
                         add_sc = add_sc,
                         remove_sc = remove_sc)

    drugsPage <- callModule(drugsPage, 'drug',
                            data_dir = data_dir,
                            new_bulk = bulkPage$new_dataset,
                            pert_query_dir = pert_query_dir,
                            pert_signature_dir = pert_signature_dir)


}
