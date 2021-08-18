tabs <- getShinyOption('tabs', c('Single Cell', 'Bulk Data', 'Drugs'))
data_dir <- getShinyOption('data_dir')
logout_url <- getShinyOption('logout_url')
is_example <- getShinyOption('is_example')
active <- tabs[1]


bootstrapPage(
  useShinyjs(),
  rintrojs::introjsUI(),
  # scrollspy for docs tab
  extendShinyjs(text = "shinyjs.init = function() {$('body').scrollspy({ target: '.bs-docs-sidenav', offset: 60 });}", functions = 'init'),
  includeScript(path = 'www/renderSelectize.js'),
  includeScript(path = 'www/progressBinding.js'),
  includeScript(path = 'www/isMobile.js'),
  includeScript(path = 'www/toggleClinicalTitle.js'),
  includeScript(path = 'www/contextMenu.js'),
  includeScript(path = 'www/anchor-polyfill.js'),
  includeCSS(path = 'www/custom.css'),
  includeCSS(path = 'www/bs-docs.css'),
  includeCSS(path = 'www/drugs.css'),
  includeCSS(path = 'www/pathways.css'),
  if (!is.null(logout_url)) tags$head(includeHTML("www/google-analytics.js")),
  tags$head(HTML("<title>Dseqr</title>"),
            tags$link(rel = "icon", type = "image/png", href = "https://raw.githubusercontent.com/hms-dbmi/dseqr.sp/master/favicon.png")),
  navbarUI(tabs, active, logout_url),
  navbar2UI(is_example),

  fluidPage(
    tags$div(
      class = "tab-content shiny-bound-input", `data-tabsetid` = "tabset", id = "tabs",

      # tabs
      scPageUI("sc", tab = 'Single Cell', active),
      bulkPageUI('bulk', tab = 'Bulk Data', active),
      drugsPageUI("drug", tab = 'Drugs', active),
      docsPageUI('docs', tab = 'Docs', active)
    )
  )
)



