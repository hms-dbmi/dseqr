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

    is_local <- getShinyOption('is_local', TRUE)

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
    observe(toggleClass('start_tour', 'invisible', condition = input$tab == 'Docs'))


    # rintrojs
    observeEvent(input$start_tour, {
        if (input$tab == 'Single Cell') {
            steps <- utils::read.csv('www/sc_intro.csv', stringsAsFactors = FALSE)

        } else if (input$tab == 'Bulk Data') {
            steps <- utils::read.csv('www/bulk_intro.csv', stringsAsFactors = FALSE)

        } else if (input$tab == 'Drugs') {
            steps <- utils::read.csv('www/drugs_intro.csv', stringsAsFactors = FALSE)

        } else {
            print(input$tab)
            return(NULL)
        }

        rintrojs::introjs(session,
                          options = list(
                              showStepNumbers = 'false',
                              steps = steps))

    })

    observe({
        if (length(list.dirs(sc_dir, recursive = FALSE))) return(NULL)
        # show hints and add dataset modal if no datasets
        shinyjs::click('add_dataset')
        rintrojs::hintjs(session,
                         options = list(hints =
                                            data.frame(
                                                element = '#docs-link',
                                                hint = 'Read the docs for all the details.',
                                                hintPosition = 'middle-middle')))
    })

    # open modal selectors
    observeEvent(input$feedback, {
        showModal(feedbackModal(session))
    })


    observeEvent(input$submit_feedback, {

        user <- Sys.getenv('SHINYPROXY_USERNAME', 'localhost')

        project <- rev(strsplit(data_dir, '/')[[1]])[1]
        project <- ifelse(project == user, 'private', project)

        slack <- readRDS(system.file('extdata/slack.rds', package = 'dseqr'))

        httr::POST(url = slack$feedback,
                   httr::add_headers('Content-Type' = 'application/json'),
                   body = sprintf(
                       '{"text": "🧑 💬 \n\n>_%s_ \n\n *project*: %s \n *user*: %s"}',
                       input$feedback_text,
                       project,
                       user
                   ))

        updateTextAreaInput(session, 'feedback_text', value = '', placeholder = 'Thank you!')
        shinyjs::delay(1000, {updateTextAreaInput(session, 'feedback', placeholder = ''); removeModal()})
    })



    observe({
        toggle('datasets_dropdown', condition = input$tab != 'Drugs')
        toggle('integrate_dataset', condition = input$tab == 'Single Cell')
        toggle('download_dataset', condition = input$tab == 'Single Cell')
    })


    is_mobile <- reactive(input$is_mobile)


    add_sc <- reactiveVal(NULL)
    remove_sc <- reactiveVal(NULL)
    integrate_sc <- reactiveVal(NULL)
    export_sc <- reactiveVal(NULL)

    add_bulk <- reactiveVal(NULL)
    remove_bulk <- reactiveVal(NULL)

    increment <- function(rval) {
        curr <- rval()
        if (is.null(curr)) curr <- 0
        rval(curr+1)
    }

    observeEvent(input$add_dataset, {
        req(!is_example)
        if (input$tab == 'Single Cell') increment(add_sc)
        if (input$tab == 'Bulk Data') increment(add_bulk)
    })

    observeEvent(input$remove_dataset, {
        req(!is_example)
        if (input$tab == 'Single Cell') increment(remove_sc)
        if (input$tab == 'Bulk Data') increment(remove_bulk)
    })

    observeEvent(input$integrate_dataset, {
        req(!is_example)
        if (input$tab == 'Single Cell') increment(integrate_sc)
    })

    observeEvent(input$export_dataset, {
        if (input$tab == 'Single Cell') increment(export_sc)
    })


    # login notification
    observe({
        req(!is_local)
        user <- Sys.getenv('SHINYPROXY_USERNAME', 'localhost')

        project <- rev(strsplit(data_dir, '/')[[1]])[1]
        project <- ifelse(project == user, 'private', project)

        slack <- readRDS(system.file('extdata/slack.rds', package = 'dseqr'))

        httr::POST(
            url = slack$logins,
            httr::add_headers('Content-Type' = 'application/json'),
            body = sprintf('{"text": "⭐⭐⭐ \n\n *project*: %s \n *user*: %s 🧑"}', project, user)
        )
    })

    #  call each page module only on first tab visit
    pages <- reactiveValues()
    tabs <- reactiveValues()

    observeEvent(input$tab, {
        if (input$tab == 'Single Cell') tabs$sc <- TRUE
        if (input$tab == 'Bulk Data') tabs$bulk <- TRUE
        if (input$tab == 'Drugs') tabs$drugs <- TRUE
    })

    observeEvent(tabs$sc, {

        pages$scPage <- callModule(
            scPage, 'sc',
            sc_dir = sc_dir,
            indices_dir = indices_dir,
            tx2gene_dir = tx2gene_dir,
            gs_dir = gs_dir,
            is_mobile = is_mobile,
            add_sc = add_sc,
            remove_sc = remove_sc,
            integrate_sc = integrate_sc,
            export_sc = export_sc)

    }, once = TRUE)


    observeEvent(tabs$bulk, {

        pages$bulkPage <- callModule(
            bulkPage, 'bulk',
            data_dir = data_dir,
            sc_dir = sc_dir,
            bulk_dir = bulk_dir,
            tx2gene_dir = tx2gene_dir,
            gs_dir = gs_dir,
            indices_dir = indices_dir,
            add_bulk = add_bulk,
            remove_bulk = remove_bulk)


    }, once = TRUE)

    observeEvent(tabs$drugs, {

        pages$drugsPage <- callModule(
            drugsPage, 'drug',
            data_dir = data_dir,
            pert_query_dir = pert_query_dir,
            pert_signature_dir = pert_signature_dir,
            tx2gene_dir = tx2gene_dir)

    }, once = TRUE)
}
