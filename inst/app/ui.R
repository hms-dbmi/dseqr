bootstrapPage(
  useShinyjs(),
  includeScript(path = 'www/contrasts.js'),
  includeCSS(path = 'www/custom.css'),
  htmlTemplate("navbar.html"),
  fluidPage(
    # make sure selectize loaded (not using default)
    tags$div(style = 'display: none', selectizeInput('blah1', label = NULL, choices = '')),

    # THE TABS ----
    tags$div(class = "tab-content", `data-tabsetid` = "1823", id = "tabs",
             # single cell tab
             scPageUI("sc")
    )
  )
)
