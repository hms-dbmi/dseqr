
datasetsPageUI <- function(id, tab, active) {
  ns <- NS(id)
  withTags({
    tabPane(tab, active,
            div(class = 'row',
                div(class = 'col-sm-6',
                    datasetsFormInput(ns('form'))
                )
            ),
            hr(),
            div(dsPairsTable(ns('pairs')))
    )
  })
}

datasetsFormInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(class = "well-form well-bg",
        dsAnalysisInput(ns('selected_anal')),
        dsEndTypeInput(ns('end_type'))
    )
  })
}

dsAnalysisInput <- function(id) {
  ns <- NS(id)

  withTags({
    textInputWithButtons(ns('anal_name'), 'Analysis name:',
                         button(id = ns('anal_dir'), type = 'button', class="shinyDirectories btn btn-default action-button shiny-bound-input",
                                `data-title` = 'Folder with fastq.gz files',
                                title = 'Select folder with fastq.gz files',
                                i(class = 'far fa-folder')
                         )
    )
  })
}

dsEndTypeInput <- function(id) {
  ns <- NS(id)

  selectizeInput(ns('end_type'),
                 'Confirm end-type:',
                 choices = NULL, width = '100%')
}

dsPairsTable <- function(id) {
  ns <- NS(id)

  tagList(
    dsLabelRowsUI(ns('label_rows')),
    DT::dataTableOutput(ns("pdata"))
  )
}

dsLabelRowsUI <- function(id) {
  ns <- NS(id)
  withTags({
    div(class = 'well-form no-well',
         div(class = 'btn-group btn-group-justified',
             div(class = 'btn-group',
                 button(type = 'button',
                        class = 'btn btn-default dropdown-toggle',
                        `data-toggle`='dropdown',
                        `aria-haspopup`='true',
                        `aria-expanded`='false',
                        span('Label Selected Rows')
                 ),
                 ul(class = 'dropdown-menu', style = 'width: 100%;',
                    li(class="dropdown-header", 'Files'),
                    dropdownMenuButton(ns('pair'), 'Paired'),
                    dropdownMenuButton(ns('rep'), 'Replicates'),
                    li(role = 'separator', class='divider'),
                    li(class="dropdown-header", 'Group'),
                    dropdownMenuButton(ns('test'), 'Test'),
                    dropdownMenuButton(ns('ctrl'), 'Control'),
                    li(role = 'separator', class='divider'),
                    dropdownMenuButton(ns('reset'), 'Reset Labels')
                 )
             )
         )

    )
  })
}

textInputWithButtons <- function(id, label, ...) {
  tags$div(class = 'form-group selectize-fh',
           tags$label(class = 'control-label', `for` = id, label),
           tags$div(class = 'input-group',
                    tags$input(id = id, type = 'text', class = 'form-control shiny-bound-input', value = '', placeholder = ''),
                    tags$span(class = 'input-group-btn', ...)
           )
  )
}

dropdownMenuButton <- function(id, label) {
  tags$li(
    tags$a(
      id = id, role = 'button', class = 'action-button shiny-bound-input', label
    )
  )
}


tabs <- c('Datasets', 'Single Cell', 'Bulk', 'Pathways', 'Drugs')
active <- 'Datasets'

bootstrapPage(
  useShinyjs(),
  includeScript(path = 'www/contrasts.js'),
  includeScript(path = 'www/toggleClinicalTitle.js'),
  includeCSS(path = 'www/custom.css'),
  includeCSS(path = 'www/drugs.css'),
  navbarUI(tabs, active),
  fluidPage(
    # make sure css/js loaded from packages where not using functions (not using default)
    tags$div(style = 'display: none;', selectizeInput('blah1', label = NULL, choices = '')),
    tags$div(style = 'display: none;', shinyFiles::shinyDirButton("blah2", title='', label='', icon=icon('plus'))),

    tags$div(class = "tab-content", `data-tabsetid` = "tabset", id = "tabs",
             # single cell tab
             scPageUI("sc", tab = 'Single Cell', active),
             bulkPageUI("bulk", tab = 'Bulk', active),
             drugsPageUI("drug", tab = 'Drugs', active),
             datasetsPageUI('datasets', tab = 'Datasets', active)
    )
  )
)
