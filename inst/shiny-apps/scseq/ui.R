## ui.R ##
bootstrapPage(
  shinyjs::useShinyjs(),
  shiny::includeScript(path = 'www/contrasts.js'),
  shiny::includeCSS(path = 'www/custom.css'),
  shiny::htmlTemplate("navbar.html"),
  shiny::fluidPage(
    shiny::tags$div(class="tab-content", `data-tabsetid`="1823",
     # template for single cell
     source('www/templates/single-cell.R')$value,
     source('www/templates/contrasts.R')$value

    )
  )
)
