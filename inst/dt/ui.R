
bootstrapPage(
  useShinyjs(),
  includeScript(path = 'www/contextMenu.js'),
  includeCSS(path = 'www/custom.css'),
  includeCSS(path = 'www/drugs.css'),
  fluidPage(
    tags$div(style = 'padding-top: 15px; padding-bottom: 15px;',
             drugsTableOutput('table')

    )
  )
)



