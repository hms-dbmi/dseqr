#' UI for Bulk Data page
#'
#' @inheritParams scPageUI
#'
#' @return shiny.tag with html for bulk data tab
#'
#' @export
bulkPageUI <- function(id, tab, active) {
  ns <- NS(id)
  withTags({
    tabPane(tab, active,
            div(class = 'row',
                div(class = 'col-sm-12 col-lg-4',
                    bulkFormInput(ns('form'))
                ),
                div(id = ns('mds_plotly_container'),
                    div(class = 'col-sm-12 col-lg-4 mobile-margin', id = 'bulk-intro-mds-left',
                        bulkPlotlyUI(ns('mds_plotly_unadjusted'))
                    ),
                    div(class = 'col-sm-12 col-lg-4 mobile-margin', id = 'bulk-intro-mds-right',
                        bulkPlotlyUI(ns('mds_plotly_adjusted'))
                    )
                ),
                div(id = ns('gene_plotly_container'), class = 'col-sm-12 col-lg-7 mobile-margin', style = 'display: none;',
                    bulkPlotlyUI(ns('gene_plotly'))
                ),
                div(id = ns('cells_plotly_container'), class = 'col-sm-12 col-lg-7 mobile-margin', style = 'display: none;overflow-x: auto;',
                    bulkPlotlyUI(ns('cells_plotly'))
                )
            ),
            hr(),
            div(id = ns('quant_table_container'), style = 'display: none;',
                bulkTable(ns('quant'))
            ),
            div(id = ns('anal_table_container'), style = 'display: none;',
                bulkAnnotInput(ns('anal')),
                bulkTable(ns('explore'))
            )
    )
  })
}

#' UI for Bulk Data annotation upload/download
#'
#' @keywords internal
#' @noRd
bulkAnnotInput <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(id=ns('validate-up'), class='validate-wrapper',
             shinyWidgets::actionGroupButtons(
               inputIds = c(ns('click_dl'), ns('click_up')),
               labels = list(icon('download', 'fa-fw'), icon('upload', 'fa-fw'))
             ),
             tags$div(class='help-block', id = ns('error_msg'))

    ),
    # hidden dl/upload buttons
    div(style = 'display: none',
        fileInput(ns('up_annot'), '', accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
    ),
    downloadLink(ns('dl_annot'), ''),

    shinyBS::bsTooltip(id = ns('click_dl'), title = 'Download metadata to fill: <b>Group name</b> and <b>Pair</b> (optional)', options = list(container = 'body')),
    shinyBS::bsTooltip(id = ns('click_up'), title = 'Upload filled metadata', options = list(container = 'body'))
  )
}


#' Plotly MDS output
#'
#' @keywords internal
#' @noRd
bulkPlotlyUI <- function(id) {
  ns <- NS(id)
  shinydlplot::downloadablePlotlyUI(ns('plotly'))
}

#' Input form for Bulk Data page
#'
#' @keywords internal
#' @noRd
bulkFormInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(class = "well-form well-bg",
        bulkDatasetInput(ns('selected_dataset')),
        div(id = ns('quant_dataset_panel'), style = 'display: none;',
            bulkFormQuantInput(ns('quant_form'))
        ),
        div(id = ns('anal_dataset_panel'), style = 'display: none;',
            bulkFormAnalInput(ns('anal_form'))
        )
    )
  })
}


#' Dataset selection input for bulkFormInput
#'
#' @keywords internal
#' @noRd
bulkDatasetInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(id = 'bulk-intro-dataset',
        shinypanel::selectizeInputWithButtons(
          ns('dataset_name'), 'Select a bulk dataset:',
          container_id = 'dataset_name_container',
          options = list(optgroupField = 'type'),
          svaButton(inputId = ns('show_nsv'), sliderId = ns('selected_nsv')),
          actionButton(ns('show_dtangle'), '', icon = icon('object-ungroup', 'far fa-fw'), title = 'Toggle cell-type deconvolution'),
          hide_btns = TRUE
        ),
        shinyFiles::shinyDirLink(ns('new_dataset_dir'), '', 'Select folder with fastq.gz files'),
        div(class = 'hidden-forms',
            dtangleFormInput(ns('dtangle'))
        )
    )
  })
}


#' Dataset quantification inputs for bulkFormInput
#'
#' @keywords internal
#' @noRd
bulkFormQuantInput <- function(id) {
  ns <- NS(id)

  tagList(
    tags$div(id = ns('bulk_controls'),
             justifiedButtonGroup(
               container_id = ns('quant_labels'),
               label = 'Label selected rows:',
               help_id = ns('error_msg'),
               actionButton(ns('pair'), 'Paired'),
               actionButton(ns('rep'), 'Replicate'),
               actionButton(ns('reset'), 'Reset')
             )
    ),
  )
}


#' Differential expression analysis inputs for bulkFormInput
#'
#' @keywords internal
#' @noRd
bulkFormAnalInput <- function(id) {
  ns <- NS(id)

  tagList(
    bulkAnalInput(ns('ds'), label = 'Download two-group comparison:'),
    div(id = 'bulk-intro-genes',
        selectizeInput(ns('explore_genes'), choices = NULL, width = "100%",
                       'Show expression for genes:',
                       options = list(maxItems = 6, multiple = TRUE))
    )
  )
}


#' Button with sliders for adjusting number of surrogate variables
#'
#' @keywords internal
#' @noRd
svaButton <- function(inputId, sliderId, max_svs = 0, prev_svs = 0) {

  dropdownButtonMod(
    br(),
    inputId = inputId,
    sliderInput(sliderId, 'Surrogate variables:',
                width = '100%', ticks = FALSE,
                min = 0, max = max_svs, value = prev_svs, step = 1),
    circle = FALSE, right = TRUE, label = tags$span('0', class='fa fa-fw'),
    title = 'Number of surrogate variables to model'
  )
}


#' Tables for datasets page
#'
#' @keywords internal
#' @noRd
bulkTable <- function(id) {
  ns <- NS(id)
  withTags({
    div(class = 'dt-container',
        DT::dataTableOutput(ns("pdata"))
    )
  })
}


#' Input form for single-cell deconvolution
#'
#' @keywords internal
#' @noRd
dtangleFormInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(id = ns('dtangle_form'), class = 'hidden-form', style = 'display: none;',
        selectizeInput(ns('dtangle_dataset'), 'Reference single-cell dataset:', choices = '', width = '100%'),
        shinypanel::selectizeInputWithButtons(
          ns('include_clusters'),
          label = 'Clusters to include:',
          label_title = 'Select cell types that are expected in the bulk dataset',
          options = list(multiple = TRUE, placeholder = 'Select none to include all'),
          actionButton(ns('submit_dtangle'), '', icon = icon('chevron-right', 'fa-fw'), title = 'Submit cell-type deconvolution')
        )
    )
  })
}


#' Bulk Differential expression analysis input
#'
#' @keywords internal
#' @noRd
bulkAnalInput <- function(id, with_dl = TRUE, label = 'Select groups to compare:') {
  ns <- NS(id)

  options <- list(maxItems = 2, placeholder = 'Select test then control group')
  if (with_dl) {
    input <- tags$div(
      id = 'bulk-intro-comparison',
      shinypanel::selectizeInputWithButtons(
        ns('contrast_groups'),
        label,
        actionButton(ns('click_dl'), label = NULL, icon = icon('download', 'fa-fw'), title = 'Download differential expression analysis'),
        options = options,
        container_id = ns('run_anal_container')
      ),
      downloadLink(ns('download'), '')

    )

  } else {
    input <- tags$div(
      id = ns('run_anal_container'),
      shinypanel::selectizeInputWithValidation(
        ns('contrast_groups'),
        label,
        options = options
      )

    )
  }

  return(input)
}



tabs <- getShinyOption('tabs', c('Single Cell', 'Bulk Data', 'Drugs'))
data_dir <- getShinyOption('data_dir')
logout_url <- getShinyOption('logout_url')
is_local <- getShinyOption('is_local')
is_example <- getShinyOption('is_example')
active <- tabs[1]



remoteDeps <- list()
if (!is_local) {

  dtCoreDeps <- htmltools::htmlDependency(
    'dt-core', '1.10.20',
    src = c(href = 'https://cdn.datatables.net/1.10.20/'),
    stylesheet = c('css/jquery.dataTables.min.css', 'css/jquery.dataTables.extra.css'),
    script = c('js/jquery.dataTables.min.js')
  )

  dtScrollerDeps <- htmltools::htmlDependency(
    'dt-ext-scroller', '1.10.20',
    src = c(href = 'https://cdn.datatables.net/scroller/2.0.1'),
    stylesheet = c('css/scroller.dataTables.min.css'),
    script = c('js/dataTables.scroller.min.js')
  )

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
    src = c(href = "https://d174upwcmdw9dj.cloudfront.net/"),
    script = c("picker.min.js")
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
    pickerDep,
    dtCoreDeps,
    dtScrollerDeps,
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



