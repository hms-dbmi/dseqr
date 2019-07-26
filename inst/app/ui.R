


tabs <- c('Datasets', 'Single Cell', 'Pathways', 'Drugs')
active <- 'Single Cell'


bootstrapPage(
  useShinyjs(),
  includeScript(path = 'www/contrasts.js'),
  includeScript(path = 'www/cellOptions.js'),
  includeScript(path = 'www/pathOptions.js'),
  includeScript(path = 'www/toggleClinicalTitle.js'),
  includeCSS(path = 'www/custom.css'),
  includeCSS(path = 'www/drugs.css'),
  includeCSS(path = 'www/pathways.css'),
  navbarUI(tabs, active),
  fluidPage(
    # make sure css/js loaded from packages where not using functions (not using default)
    tags$div(style = 'display: none;', selectizeInput('blah1', label = NULL, choices = '')),
    tags$div(style = 'display: none;', shinyFiles::shinyDirButton("blah2", title='', label='', icon=icon('plus'))),

    tags$div(class = "tab-content", `data-tabsetid` = "tabset", id = "tabs",
             # single cell tab
             scPageUI("sc", tab = 'Single Cell', active),
             drugsPageUI("drug", tab = 'Drugs', active),
             dsPageUI('datasets', tab = 'Datasets', active),
             pathPageUI('pathways', tab = 'Pathways', active)
    )
  )
)

