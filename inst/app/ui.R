tabs <- c('Bulk Data', 'Single Cell', 'Pathways', 'Drugs')
active <- 'Bulk Data'

bootstrapPage(
  useShinyjs(),
  # scrollspy for docs tab
  extendShinyjs(text = "shinyjs.init = function() {$('body').scrollspy({ target: '.bs-docs-sidenav', offset: 60 });}"),
  includeScript(path = 'www/renderSelectize.js'),
  includeScript(path = 'www/toggleClinicalTitle.js'),
  includeScript(path = 'www/contextMenu.js'),
  includeScript(path = 'www/anchor-polyfill.js'),
  includeCSS(path = 'www/custom.css'),
  includeCSS(path = 'www/bs-docs.css'),
  includeCSS(path = 'www/drugs.css'),
  includeCSS(path = 'www/pathways.css'),
  navbarUI(tabs, active),
  fluidPage(
    # make sure css/js loaded from packages where not using functions (not using default)
    tags$div(style = 'display: none;', selectizeInput('blah1', label = NULL, choices = '')),
    tags$div(style = 'display: none;', shinyFiles::shinyDirButton("blah2", title='', label='', icon=icon('plus'))),

    tags$div(class = "tab-content", `data-tabsetid` = "tabset", id = "tabs",
             bulkPageUI('bulk', tab = 'Bulk Data', active),
             scPageUI("sc", tab = 'Single Cell', active),
             pathPageUI('pathways', tab = 'Pathways', active),
             drugsPageUI("drug", tab = 'Drugs', active),
             docsPageUI('docs', tab = 'docs', active)

    )
  )
)



