#---- server

# modal to manage projects
projectModal <- function(session, choices, selected, options) {


  modalDialog(
    DT::dataTableOutput(session$ns('projects_table'), width = '100%'),
    tags$div(
      id = session$ns('validate_projects'),
      br(),
      tags$span(class = 'help-block', id = session$ns('error_msg_projects'))
    ),
    title = "Manage Your Projects",
    size = 'l',
    footer = tagList(
      actionButton(session$ns("add_project"),
                   "Add"),
      actionButton(session$ns("open_project"),
                   "Open Selected",
                   class='btn-success'),
      tags$div(class='pull-left', modalButton("Cancel"))
    ),
    easyClose = TRUE
  )
}


# modal to delete dataset
deleteProjectModal <- function(session, project) {

  modalDialog(
    tags$div(id=session$ns('confirm_delete_container'),
             tags$div(class='alert alert-danger',
                      'Delete ', tags$b(project), ' and all its datasets?',
                      br(), br(),
                      tags$b('This action cannot be undone.')
             ),
             br(),
             textInput(session$ns('confirm_delete'), HTML(paste0('<span> Type <i><span style="color: gray">',project, '</span></i> to confirm:</span>')),
                       placeholder = project, width = '100%'
             ),
    ),
    title = 'Delete Project',
    size = 'm',
    footer = tagList(
      actionButton(session$ns("delete_project"), "Delete Project"),
      tags$div(class='pull-left', modalButton("Cancel"))
    ),
    easyClose = FALSE,
  )
}

get_num_sc_datasets <- function(project, user_dir) {
  if (project == '') return(0)
  sc_dir <- file.path(user_dir, project, 'single-cell')
  dataset_names <- list.files(sc_dir)
  sum(check_has_scseq(dataset_names, sc_dir))
}

get_num_bulk_datasets <- function(project, user_dir) {
  if (project == '') return(0)
  bulk_dir <- file.path(user_dir, project, 'bulk')
  length(list.dirs(bulk_dir, recursive = FALSE))
}

server <- function(input, output, session) {


  # get arguments from calling function
  # defaults for testing
  # shiny::shinyOptions don't make it through

  # base directory contains data_dir folder
  app_name <- getShinyOption('app_name', 'test_user')
  data_dir <- getShinyOption('data_dir', 'tests/testthat/test_data_dir')

  # path where pert queries will be stored
  pert_query_dir <- getShinyOption('pert_query_dir', file.path(data_dir, '.pert_query_dir'))

  # path where pert signatures will be stored
  pert_signature_dir <- getShinyOption('pert_signature_dir', file.path(data_dir, '.pert_signature_dir'))

  # path where kallisto index is downloaded and stored
  indices_dir <- getShinyOption('indices_dir', file.path(data_dir, '.indices_dir'))

  # path where save tx2genes
  tx2gene_dir <- getShinyOption('tx2gene_dir', file.path(data_dir, '.tx2gene_dir'))

  # path where gene set data is stored
  gs_dir <- getShinyOption('gs_dir', file.path(data_dir, '.gs_dir'))

  is_example <- getShinyOption('is_example', FALSE)
  is_local <- getShinyOption('is_local', TRUE)

  # reset testing data
  if (isTRUE(getOption('shiny.testmode'))) {
    unlink(data_dir, recursive = TRUE)
    dir.create(data_dir, recursive = TRUE)
  }

  # ensure various directories exist
  app_dirs <- c(pert_query_dir, pert_signature_dir, indices_dir, tx2gene_dir, gs_dir)
  for (dir in app_dirs) dir.create(dir, showWarnings = FALSE)

  # on remote: send errors to slack
  if (!is_local) {
    options(shiny.error = function() {
      observe({
        user_name <- user_name()
        send_slack_error(app_name, user_name)
      })
    })
  }


  # hide tour button for docs page
  observe(shinyjs::toggleClass('start_tour', 'invisible', condition = input$tab == 'Docs'))

  user_name <- reactive({
    if (app_name == 'example' || is_local) return(app_name)

    # app_name is 'private'
    user_name <- session$request$HTTP_X_SP_USERID
    return(user_name)
  })

  user_dir <- reactive({
    user_name <- user_name()
    user_dir <- file.path(data_dir, user_name)

    if (!dir_exists(user_dir))
      init_dseqr(user_name, data_dir)

    return(user_dir)
  })

  # rintrojs
  observeEvent(input$start_tour, {
    if (input$tab == 'Bulk Data')
      steps <- utils::read.csv('www/bulk_intro.csv', stringsAsFactors = FALSE)
    else if (input$tab == 'Drugs')
      steps <- utils::read.csv('www/drugs_intro.csv', stringsAsFactors = FALSE)

    rintrojs::introjs(session, options = list(showStepNumbers = 'false', steps = steps))
  })

  observeEvent(input$tour_sc_clusters, {
    steps <- utils::read.csv('www/sc_intro_clusters.csv', stringsAsFactors = FALSE)
    rintrojs::introjs(session, options = list(showStepNumbers = 'false', steps = steps))
  })

  observeEvent(input$tour_sc_samples, {
    steps <- utils::read.csv('www/sc_intro_samples.csv', stringsAsFactors = FALSE)
    steps$step <- seq_len(nrow(steps))
    rintrojs::introjs(session, options = list(showStepNumbers = 'false', steps = steps))
  })

  # customize dataset management dropdown for each tab
  observe({
    toggle('tour_dropdown', condition = input$tab == 'Single Cell')
    toggle('start_tour_container', condition = input$tab != 'Single Cell')

    toggle('add_dataset', condition = input$tab != 'Drugs')
    toggle('integrate_dataset', condition = input$tab == 'Single Cell')

    toggle('export_dataset', condition = input$tab == 'Single Cell')
    toggle('remove_dataset', condition = input$tab != 'Drugs')

    toggle('dataset_management_sep1', condition = input$tab != 'Drugs')
    toggle('dataset_management_sep2', condition = input$tab != 'Drugs')
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

  # selecting the project
  project_choices <- reactiveVal()

  observe({
    choices <- list.dirs(user_dir(), recursive = FALSE, full.names = FALSE)
    project_choices(choices)
  })

  observeEvent(input$select_project, {
    req(!is_example)
    showModal(projectModal(session))
  })


  projects_table <- reactive({
    projects <- project_choices()
    project <- project()
    req(project)

    nsc <- sapply(projects, get_num_sc_datasets, user_dir())
    nbulk <- sapply(projects, get_num_bulk_datasets, user_dir())

    df <- data.frame(
      ` ` = getDeleteRowButtons(session, length(projects), title = 'Delete project'),
      'Project' = projects,
      'Single Cell Datasets' = nsc,
      'Bulk Datasets' = nbulk,
      selected = ifelse(projects == project, 'hl', 'other'),
      check.names = FALSE,
      row.names = NULL
    )

    return(df)
  })

  output$projects_table <- DT::renderDataTable({
    dt <- isolate(projects_table())

    DT::datatable(
      dt,
      class = 'cell-border',
      rownames = FALSE,
      escape = FALSE, # to allow HTML in table
      selection = list(mode = 'single'),
      editable = list(target = "cell", disable = list(columns = c(0, 2, 3))),
      options = list(
        scrollX = TRUE,
        ordering = FALSE,
        dom = 't',
        paging = FALSE,
        columnDefs = list(list(visible = FALSE, targets = 4))
      )) %>%
      DT::formatStyle(
        "selected",
        target = "row",
        backgroundColor = DT::styleEqual(c('hl', 'other'), values = c('#FFFFED', 'white'))
      )
  })


  # add row to projects table
  proxy <- DT::dataTableProxy('projects_table')
  observeEvent(input$add_project, {
    project_choices(c(project_choices(), ''))
  })

  # update projects table
  observe({
    editn()
    DT::replaceData(proxy, projects_table(), rownames = FALSE)
  })


  validate_open_project <- function(choices, row) {
    if (!length(row)) return('Select a project row')
    if (choices[row] == '') return('Add project name (double click cell to edit)')

    return(NULL)
  }

  error_msg <- reactiveVal()

  observe({
    msg <- error_msg()
    toggleClass('validate_projects', class = 'has-error', condition = !is.null(msg))
    # show error message
    html('error_msg_projects', html = msg)
  })

  # open selected project
  observeEvent(input$open_project, {
    row <- input$projects_table_rows_selected
    choices <- project_choices()

    msg <- validate_open_project(choices, row)
    error_msg(msg)
    if (!is.null(msg)) {
      return(NULL)
    }

    selected <- choices[row]
    project(selected)
    removeModal()
  })

  validate_edit_project_name <- function(choices, prev, new) {
    msg <- NULL
    if (new %in% choices) msg <- 'Project name already exists'
    return(msg)

  }

  delete_candidate <- reactiveVal()

  validate_delete_project <- function(df, row, sel) {
    if (nrow(df) == 1) return('Need at least one project')
    return(NULL)
  }

  # delete project
  observeEvent(input$delete_row, {

    row <- as.numeric(strsplit(input$delete_row, "_")[[1]][2])
    df <- projects_table()
    sel <- df$Project[row]

    msg <- validate_delete_project(df, row, sel)
    error_msg(msg)
    if (!is.null(msg)) return(NULL)

    # removing empty
    if (sel == '') {
      choices <- project_choices()[-row]
      project_choices(choices)
      return(NULL)
    }

    delete_candidate(sel)

    showModal(deleteProjectModal(session, project = sel))
  })

  allow_delete <- reactive(input$confirm_delete == delete_candidate())

  observe({
    shinyjs::toggleState('delete_project', condition = allow_delete())
    shinyjs::toggleClass('delete_project', class = 'btn-danger', condition = allow_delete())
  })

  shinyjs::toggleState('select_project', condition = !is_example)

  observeEvent(input$delete_project, {
    # remove from choices
    del <- delete_candidate()
    choices <- project_choices()
    choices <- choices[choices != del]
    project_choices(choices)

    # removing selected
    if (del == project()) {
      project(choices[1])
    }

    # remove data
    unlink(file.path(user_dir(), del), recursive = TRUE)
    removeModal()
  })

  # edit selected project name
  editn <- reactiveVal(0)
  observeEvent(input$projects_table_cell_edit, {
    editn(editn()+1)
    info <- input$projects_table_cell_edit
    choices <- project_choices()

    prev <- choices[info$row]
    new <- info$value
    msg <- validate_edit_project_name(choices, prev, new)
    error_msg(msg)

    if (!is.null(msg)) return(NULL)


    choices[info$row] <- new
    project_choices(choices)

    new_dir <- file.path(user_dir(), new)

    if (prev == "") {
      dir.create(new_dir, showWarnings = FALSE)
      return(NULL)
    }

    prev_dir <- file.path(user_dir(), prev)

    if (dir_exists(prev_dir))
      file.rename(prev_dir, new_dir)

    if (prev == project()) project(new)

  })

  observe({
    project <- project()
    req(project)
    qs::qsave(project, prev_path())
  })


  project_dir <- reactive(file.path(user_dir(), project()))
  prev_path <- reactive(file.path(user_dir(), 'prev_project.qs'))

  project <- reactiveVal()

  observe({
    project(qread.safe(prev_path(), .nofile = 'default'))
  })

  sc_dir <- reactive({
    project <- project()
    req(project)

    sc_dir <- file.path(user_dir(), project, 'single-cell')
    dir.create(sc_dir, showWarnings = FALSE)
    return(sc_dir)
  })

  bulk_dir <- reactive({
    project <- project()
    req(project)

    bulk_dir <- file.path(user_dir(), project, 'bulk')
    dir.create(bulk_dir, showWarnings = FALSE)
    return(bulk_dir)
  })


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
    if (user == 'alexvpickering@gmail.com') return(NULL)

    project <- project()
    slack <- readRDS(system.file('extdata/slack.rds', package = 'dseqr'))

    workspace <- basename(user_dir())
    workspace <- ifelse(workspace == user, 'private', workspace)

    httr::POST(
      url = slack$logins,
      httr::add_headers('Content-Type' = 'application/json'),
      body = sprintf('{"text": "\U2B50\U2B50\U2B50 \n\n *workspace*: %s \n *project*: %s \n *user*: %s \U1F9D1"}', workspace, project, user)
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
      project_dir = project_dir,
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
      project_dir = project_dir,
      pert_query_dir = pert_query_dir,
      pert_signature_dir = pert_signature_dir,
      tx2gene_dir = tx2gene_dir)

  }, once = TRUE)
}
