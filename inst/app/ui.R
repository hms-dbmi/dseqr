# defaults for testing
# shiny::shinyOptions don't make it through

tabs <- getShinyOption('tabs', c('Single Cell', 'Bulk Data', 'Drugs'))
logout_url <- getShinyOption('logout_url', NULL)
is_local <- getShinyOption('is_local', TRUE)
is_example <- getShinyOption('is_example', FALSE)
active <- tabs[1]

remoteDeps <- list()
if (!is_local) {

  selectizeDep <- htmltools::htmlDependency(
    'selectize', '0.12.4',
    src = c(href = 'https://d174upwcmdw9dj.cloudfront.net/'),
    stylesheet = 'selectize.bootstrap3.min.css',
    script = c('selectize.min.js', 'selectize-plugin-a11y.js')
  )

  shinypanelDep <- htmltools::htmlDependency(
    'shinypanel', '0.1.5',
    src = c(href = 'https://d174upwcmdw9dj.cloudfront.net/'),
    stylesheet = c('shinypanel.css')
  )

  dseqrDep <- htmltools::htmlDependency(
    'dseqr', '0.30.2',
    src = c(href = 'https://d174upwcmdw9dj.cloudfront.net/'),
    stylesheet = c('custom.css')
  )

  htmlwidgetsDep <- htmltools::htmlDependency(
    'htmlwidgets', '1.5.3',
    src = c(href = "https://d174upwcmdw9dj.cloudfront.net/"),
    script = c("htmlwidgets.js")
  )

  pickerDep <- htmltools::htmlDependency(
    'picker', '0.2.5',
    src = c(href = "https://d174upwcmdw9dj.cloudfront.net"),
    script = c("picker4.min.js")
  )

  deckglDep <- htmltools::htmlDependency(
    'deck.gl', '8.1.4',
    src = c(href = "https://d174upwcmdw9dj.cloudfront.net/"),
    script = c("deckgl.min.js")
  )

  bootstrapDep <- htmltools::htmlDependency(
    "bootstrap", "3.4.1",
    src = c(href = "https://d174upwcmdw9dj.cloudfront.net/"),
    script = "bootstrap.min.js", stylesheet = "bootstrap.min.css"
  )

  fontawesomeDep <- htmltools::htmlDependency(
    'font-awesome', '5.13.0',
    src = c(href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.13.0/css/"),
    stylesheet = c("all.min.css", "v4-shims.min.css")
  )

  bootstrapSelectDep <- htmltools::htmlDependency(
    'bootstrap-select', '1.13.18',
    src = c(href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap-select/1.13.18/"),
    stylesheet = c("css/bootstrap-select.min.css"),
    script = c("js/bootstrap-select.min.js")
  )

  remoteDeps <- list(
    selectizeDep,
    bootstrapDep,
    fontawesomeDep,
    shinypanelDep,
    dseqrDep,
    htmlwidgetsDep,
    deckglDep,
    pickerDep,
    bootstrapSelectDep
  )
}

checkBulkFileName <- 'function checkBulkFileName(fieldObj, shinyId, testString) {
    var fileName  = fieldObj.value;
    var fileBase = fileName.split(/[\\\\/]/).pop();

    const re = new RegExp(testString);

    console.log(fileBase)
    if (!re.test(fileBase)) {
        fieldObj.value = "";
        Shiny.setInputValue(shinyId, "", {priority: "event"})
        return false;
    }
    return true;
}'



checkSingleCellFileName <- 'function checkSingleCellFileName(fieldObj, shinyId) {
    var fileName  = fieldObj.value;
    var fileBase = fileName.split(/[\\\\/]/).pop();

    if (!/[.]h(df)?5(.gz)?$|[.]tsv(.gz)?$|[.]mtx(.gz)?$|[.]rds$|[.]qs$/.test(fileBase)) {
        fieldObj.value = "";
        Shiny.setInputValue(shinyId, "", {priority: "event"})
        return false;
    }
    return true;
}'

bootstrapPage(
  remoteDeps,
  if (!is_local) includeHTML("www/gtm.html"),
  if (!is_local) tags$head(includeHTML("www/gtag.html")),
  useShinyjs(),
  rintrojs::introjsUI(),
  includeScript(path = 'www/renderSelectize.js'),
  includeScript(path = 'www/isMobile.js'),
  includeScript(path = 'www/contextMenu.js'),
  tags$script(checkBulkFileName),
  tags$script(checkSingleCellFileName),
  includeCSS(path = 'www/custom.css'),
  includeCSS(path = 'www/drugs.css'),
  includeCSS(path = 'www/pathways.css'),
  tags$head(HTML("<title>Dseqr</title>"),
            tags$link(rel = "icon", type = "image/png", href = "https://raw.githubusercontent.com/hms-dbmi/dseqr.sp/master/favicon.png")),


  navbarUI(tabs, active, logout_url),
  navbar2UI(is_example),

  fluidPage(
    tags$div(
      class = "tab-content shiny-bound-input", `data-tabsetid` = "tabset",

      # tabs
      scPageUI("sc", tab = 'Single Cell', active),
      bulkPageUI('bulk', tab = 'Bulk Data', active),
      drugsPageUI("drug", tab = 'Drugs', active)
    )
  )
)



