tabs <- getShinyOption('tabs', c('Single Cell', 'Bulk Data', 'Drugs'))
data_dir <- getShinyOption('data_dir')
logout_url <- getShinyOption('logout_url')
is_local <- getShinyOption('is_local')
is_example <- getShinyOption('is_example')
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
    'dseqr', '0.20.27',
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
    script = c("picker3.min.js")
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

checkjs <- 'function checkFileName(fieldObj, shinyId) {
    var fileName  = fieldObj.value;
    var fileBase = fileName.split(/[\\\\/]/).pop();

    if (!fileBase.endsWith(".fastq.gz")) {
        fieldObj.value = "";
        Shiny.setInputValue(shinyId, "", {priority: "event"})
        return false;
    }
    return true;
}'





#' UI for navbar
#'
#' @param tabs Character vector of tab names to display
#' @param active the active tab name
#' @param logout_url url to logout from. NULL results in no Logout link
#' @return shiny.tag with html for navbar
#'
#' @export
navbarUI <- function(tabs, active, logout_url = NULL) {

  logout_li <- NULL
  if (!is.null(logout_url)) {
    # target _top break out of iframe for ShinyProxy
    logout_li <- tags$li(class = 'navbar-right', `data-toggle`="collapse",
                         `data-target`=".navbar-collapse.in", a(target="_top", href = logout_url, 'Sign out')
    )

  }

  tags$nav(class = 'navbar navbar-default navbar-static-top',
           div(class = 'container-fluid',
               div(class = 'navbar-header',
                   tags$button(type='button', class='navbar-toggle collapsed', `data-toggle`='collapse', `data-target`='#bs-navbar', `aria-expanded`='false',
                               span(class = 'sr-only', 'Toggle navigation'),
                               span(class = 'icon-bar'),
                               span(class = 'icon-bar'),
                               span(class = 'icon-bar')
                   ),
                   span(class = 'navbar-brand', title = 'dseqr',
                        tags$a(class = 'brand-icons',
                               href="/",
                               tags$img(src="favicon.png"),
                               span('seqr')
                        )
                   )
               ),
               div(id = 'bs-navbar', class = 'collapse navbar-collapse',
                   tags$ul(class = 'nav navbar-nav shiny-tab-input shiny-bound-input', `data-tabsetid` = 'tabset', id = 'tab',
                           lapply(seq_along(tabs), function(i) {
                             tab <- tabs[i]
                             is.active <- tab == active
                             tags$li(class = ifelse(is.active, 'active', ''), `data-toggle`="collapse", `data-target`=".navbar-collapse.in",
                                     a(href = paste0('#', id_from_tab(tab)), `data-toggle` = 'tab', `data-value` = tab, `aria-expanded` = ifelse(is.active, 'true', 'false'), tab)
                             )
                           }),
                           # github linkout section
                           tags$li(class = 'navbar-right', `data-toggle`="collapse", `data-target`=".navbar-collapse.in",
                                   a(href = 'https://github.com/hms-dbmi/dseqr', target="_blank", icon('github'), style = 'padding-bottom: 0px; font-size:17px;',)
                           ),
                           # docs section
                           tags$li(class = 'navbar-right', `data-toggle`="collapse", `data-target`=".navbar-collapse.in", id = 'docs-link',
                                   a(href = "https://docs.dseqr.com/docs/single-cell/add-dataset/", `target` = '_blank', 'Docs')
                           ),
                           # logout section
                           logout_li
                   )
               )
           )
  )
}

#' UI for secondary navbar
#'
#' @param hide should delete and add be disabled (for demo)
#' @return shiny.tag with html for secondary navbar
#'
#' @export
navbar2UI <- function(hide) {
  class <- 'action-button shiny-bound-input btn-intro btn navbar-btn'
  add <- ifelse(hide, 'disabled', '')

  ui <- tags$nav(
    class = 'navbar navbar-default secondary-navbar navbar-expand',
    tags$div(
      class = 'container-fluid',
      tags$div(
        tags$ul(
          class = 'nav navbar-non-responsive',
          tags$div(class = "secondary-navbar-btn-group",
                   tags$li(tags$button(class = 'btn', id = 'start_tour', class=class, icon('info', 'fa-fw'), 'Tour')),
                   tags$li(tags$button(class = 'btn', id = 'feedback', class=class, tags$i(class= 'far fa-comment-dots fa-fw'), 'Report Issue'),),
          ),
          tags$div(class='btn-group',
                   tags$button(
                     id = 'datasets_dropdown',
                     class = paste0(class, ' dropdown-toggle'),
                     `data-toggle` = 'dropdown',
                     type = 'button',
                     haspopup = 'true',
                     `aria-expanded`='false',
                     tags$i(class= 'far fa-folder-open fa-fw'),
                     tags$span(class = 'hidden-xs', 'Datasets'),
                     tags$span(class='caret')),
                   tags$ul(class="dropdown-menu",
                           tags$li(tags$a(id = 'add_dataset', role='button', class = 'action-button shiny-bound-input', icon('plus', 'fa-fw'), 'Import')),
                           tags$li(tags$a(id = 'integrate_dataset', role='button', class = 'action-button shiny-bound-input', tags$i(class= 'far fa-object-group fa-fw'), 'Integrate')),
                           tags$li(tags$a(id = 'export_dataset', role='button', class = 'action-button shiny-bound-input', tags$i(class = 'fab fa-r-project'), 'Export')),
                           tags$li(role = 'separator', class='divider'),
                           tags$li(tags$a(id = 'remove_dataset', role='button', class = 'action-button shiny-bound-input', tags$i(class= 'far fa-trash-alt fa-fw'), 'Delete'))
                   )
          )
        )
      )
    )
  )


  return(ui)
}


#' UI for a tab pane
#'
#' @param tab The name of the tab
#' @param active The name of the active tab
#' @param ... The UI elements to place in the tab
#' @return shiny div tag with UI for tab
#'
#' @keywords internal
tabPane <- function(tab, active, ...) {
  active_class <- ifelse(tab == active, 'active', '')
  tags$div(class = paste('tab-pane', active_class), `data-value` = tab, id = id_from_tab(tab), ...)
}

#' Convert tab name to formated id
#'
#' used by navbarUI and *PageUI for dseqr app
#'
#' @param tab The name of the tab (e.g. \code{'Single Cell'})
#'
#' @keywords internal
id_from_tab <- function(tab) {
  id <- tolower(tab)
  gsub(' ', '-', id)
}


#' Full width button group with validation
#'
#' @param ... actionButtons
#' @param label Character vector. Label above button group.
#' @param container_id id of container. Used to toggle 'has-error' class with shinyjs::toggleClass.
#' @param help_id id of help block. Used to show help message in change shinyjs::html
#' @param class classes to add to btn-group divs.
#'
#' @keywords internal
justifiedButtonGroup <- function(..., label, container_id = NULL, help_id = NULL, class = '') {
  btns <- list(...)

  if (length(class) == 1) class <- rep(class, length(btns))

  tags$div(class = 'form-group selectize-fh', id = container_id,
           tags$label(class = 'control-label',  label),
           tags$div(class = 'btn-group btn-group-justified', role = 'group',
                    lapply(seq_along(btns), function(i) {
                      if (is.null(btns[[i]])) return(NULL)

                      tags$div(
                        class = paste('btn-group', class[i]),
                        role = 'group',
                        id = paste0(container_id, '-', i),
                        btns[[i]]
                      )
                    })
           ),
           span(id = help_id, class = 'help-block')
  )
}


# deprecated but keep because nice design
dsLabelRowsUI <- function(id) {
  ns <- NS(id)
  div(class = 'btn-group btn-group-justified',
      div(class = 'btn-group',
          tags$button(type = 'button',
                      class = 'btn btn-default dropdown-toggle',
                      `data-toggle`='dropdown',
                      `aria-haspopup`='true',
                      `aria-expanded`='false',
                      span('Label Selected Rows')
          ),
          tags$ul(class = 'dropdown-menu', style = 'width: 100%;',
                  tags$li(class="dropdown-header", 'Files'),
                  dropdownMenuButton(ns('pair'), 'Paired'),
                  dropdownMenuButton(ns('rep'), 'Replicates'),
                  tags$li(role = 'separator', class='divider'),
                  tags$li(class="dropdown-header", 'Group'),
                  dropdownMenuButton(ns('test'), 'Test'),
                  dropdownMenuButton(ns('ctrl'), 'Control'),
                  tags$li(role = 'separator', class='divider'),
                  dropdownMenuButton(ns('reset'), 'Reset Labels')
          )
      )
  )
}


#' Dropdown menu button for dsLabelRowsUI
#'
#' @keywords internal
dropdownMenuButton <- function(id, label) {
  tags$li(
    tags$a(
      id = id, role = 'button', class = 'action-button shiny-bound-input', label
    )
  )
}


bootstrapPage(
  remoteDeps,
  if (!is_local) includeHTML("www/gtm.html"),
  if (!is_local) tags$head(includeHTML("www/gtag.html")),
  useShinyjs(),
  rintrojs::introjsUI(),
  includeScript(path = 'www/renderSelectize.js'),
  includeScript(path = 'www/isMobile.js'),
  includeScript(path = 'www/contextMenu.js'),
  tags$script(checkjs),
  includeCSS(path = 'www/custom.css'),
  includeCSS(path = 'www/drugs.css'),
  includeCSS(path = 'www/pathways.css'),
  tags$head(HTML("<title>Dseqr</title>"),
            tags$link(rel = "icon", type = "image/png", href = "https://raw.githubusercontent.com/hms-dbmi/dseqr.sp/master/favicon.png")),


  navbarUI(tabs, active, logout_url),
  navbar2UI(is_example),

  fluidPage(
    tags$div(
      class = "tab-content shiny-bound-input", `data-tabsetid` = "tabset", id = "tab",

      # tabs
      scPageUI("sc", tab = 'Single Cell', active),
      bulkPageUI('bulk', tab = 'Bulk Data', active),
      drugsPageUI("drug", tab = 'Drugs', active),
      docsPageUI('docs', tab = 'Docs', active)
    )
  )
)



