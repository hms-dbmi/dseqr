tabs <- getShinyOption('tabs')
active <- tabs[1]

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
    tags$div(class = "tab-content", `data-tabsetid` = "tabset", id = "tabs",
             scPageUI("sc", tab = 'Single Cell', active),
             bulkPageUI('bulk', tab = 'Bulk Data', active),
             drugsPageUI("drug", tab = 'Drugs', active)
             # docsPageUI('docs', tab = 'docs', active)

    )
  )
)



