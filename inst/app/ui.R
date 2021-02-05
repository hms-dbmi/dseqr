tabs <- getShinyOption('tabs')
data_dir <- getShinyOption('data_dir')
with_logout <- grepl('^/srv/drugseqr', data_dir)
active <- tabs[1]


bootstrapPage(
  useShinyjs(),
  rintrojs::introjsUI(),
  # scrollspy for docs tab
  extendShinyjs(text = "shinyjs.init = function() {$('body').scrollspy({ target: '.bs-docs-sidenav', offset: 60 });}", functions = 'init'),
  includeScript(path = 'www/renderSelectize.js'),
  includeScript(path = 'www/toggleClinicalTitle.js'),
  includeScript(path = 'www/contextMenu.js'),
  includeScript(path = 'www/anchor-polyfill.js'),
  includeCSS(path = 'www/custom.css'),
  includeCSS(path = 'www/bs-docs.css'),
  includeCSS(path = 'www/drugs.css'),
  includeCSS(path = 'www/pathways.css'),
  tags$head(HTML("<title>drugseqr</title>"),
            tags$link(rel = "icon", type = "image/png", href = "https://raw.githubusercontent.com/hms-dbmi/drugseqr.sp/master/favicon.png")),
  navbarUI(tabs, active, with_logout),
  fluidPage(
    tags$div(class = "row", style='height: 0px',
             tags$div(
               class = "lds-ripple col-lg-1 col-centered",
               id = 'loading_page',
               tags$div(),
               tags$div()),
    ),
    tags$div(style='display: none',
      id = 'main_page',
      tags$div(
        class = "tab-content shiny-bound-input", `data-tabsetid` = "tabset", id = "tabs",

        # rintrojs stuff
        span(id = 'start_tour', class='action-button shiny-bound-input btn-intro-icon btn-intro', icon('info', 'fa-fw fa-w-6')),
        shinyBS::bsTooltip(id = 'start_tour', title = 'Tour this page', placement = 'right', options = list(container = 'body')),

        # tabs
        scPageUI("sc", tab = 'Single Cell', active),
        bulkPageUI('bulk', tab = 'Bulk Data', active),
        drugsPageUI("drug", tab = 'Drugs', active),
        docsPageUI('docs', tab = 'Docs', active)

      )
    )


  )
)



