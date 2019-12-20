#' UI for Bulk Data page
#' @export
#' @keywords internal
bulkPageUI <- function(id, tab, active) {
  ns <- NS(id)
  withTags({
    tabPane(tab, active,
            div(class = 'row',
                div(class = 'col-sm-5',
                    bulkFormInput(ns('form'))
                ),
                div(id = ns('mds_plotly_container'), class = 'col-sm-7 mobile-margin', style = '',
                    bulkPlotlyUI(ns('mds_plotly'))
                ),
                div(id = ns('gene_plotly_container'), class = 'col-sm-7 mobile-margin', style = 'display: none;',
                    bulkPlotlyUI(ns('gene_plotly'))
                ),
                div(id = ns('cells_plotly_container'), class = 'col-sm-7 mobile-margin', style = 'display: none;',
                    bulkPlotlyUI(ns('cells_plotly'))
                )
            ),
            hr(),
            div(id = ns('quant_table_container'), style = 'display: none;',
                bulkTable(ns('quant'))
            ),
            div(id = ns('anal_table_container'), style = 'display: none;',
                bulkTable(ns('explore'))
            )
    )
  })
}


#' Plotly MDS output
#' @export
#' @keywords internal
bulkPlotlyUI <- function(id) {
  ns <- NS(id)
  downloadablePlotlyUI(ns('plotly'))
}

#' Input form for Bulk Data page
#' @export
#' @keywords internal
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
#' @export
#' @keywords internal
bulkDatasetInput <- function(id) {
  ns <- NS(id)



  withTags({
    div(
      selectizeInputWithButtons(ns('dataset_name'), 'Select a dataset:',
                                options = list(create = TRUE, placeholder = 'Type name to add new dataset', optgroupField = 'type'),
                                button(id = ns('dataset_dir'), type = 'button', class="shinyDirectories btn btn-default action-button shiny-bound-input disabled",
                                       `data-title` = 'Folder with fastq.gz or cellranger files',
                                       title = 'Select folder with fastq.gz or cellranger files',
                                       `parent-style` = 'min-width: 0px;width: 0px;',
                                       i(class = 'far fa-folder fa-fw')

                                ),
                                svaButton(inputId = ns('show_nsv'), sliderId = ns('selected_nsv')),
                                actionButton(ns('show_dtangle'), '', icon = icon('object-ungroup', 'far fa-fw'), title = 'Toggle cell-type deconvolution')


      ),
      div(class = 'hidden-forms',
          dtangleFormInput(ns('dtangle'))
      )
    )
  })
}


#' Dataset quantification inputs for bulkFormInput
#' @export
#' @keywords internal
bulkFormQuantInput <- function(id) {
  ns <- NS(id)

  tagList(
    tags$div(id = ns('bulk_controls'),
             bulkEndTypeInput(ns('end_type')),
             justifiedButtonGroup(
               container_id = ns('quant_labels'),
               label = 'Label selected rows:',
               help_id = ns('error_msg'),
               actionButton(ns('pair'), 'Paired'),
               actionButton(ns('rep'), 'Replicate'),
               actionButton(ns('reset'), 'Reset')
             )
    ),
    actionButton(ns('run_quant'), 'Run Quantification', width = '100%', class = 'btn-warning')
  )
}

#' Dataset end-type input for bulkFormQuantInput
#' @export
#' @keywords internal
bulkEndTypeInput <- function(id) {
  ns <- NS(id)

  selectizeInput(ns('end_type'),
                 'Confirm end-type:',
                 choices = NULL, width = '100%')
}

#' Differential expression analysis inputs for bulkFormInput
#' @export
#' @keywords internal
bulkFormAnalInput <- function(id) {
  ns <- NS(id)

  tagList(
    bulkAnalInput(ns('ds')),
    selectizeInput(ns('explore_genes'), choices = NULL, width = "100%",
                   'Show expression for genes:',
                   options = list(maxItems = 6, multiple = TRUE)),


    textInputWithButtons(ns('explore_group_name'),
                         container_id = ns('validate'),
                         'Group name for selected rows:',
                         actionButton(ns('grouped'), '', icon = icon('plus', 'fa-fw'), title = 'Add group'),
                         actionButton(ns('reset_explore'), '', icon = tags$i(class='fa-trash-alt far fa-fw'), title = 'Clear all groups'),
                         actionButton(ns('run_sva'), '', icon = icon('redo-alt', 'fa-fw'), title = 'Recalculate surrogate variables'),
                         help_id = ns('error_msg')
    )
  )
}


#' Button with sliders for adjusting number of surrogate variables
#' @export
#' @keywords internal
svaButton <- function(inputId, sliderId, max_svs = 0, prev_svs = 0) {

  drugseqr::dropdownButton(
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
#' @export
#' @keywords internal
bulkTable <- function(id) {
  ns <- NS(id)
  withTags({
    div(class = 'dt-container',
        DT::dataTableOutput(ns("pdata"))
    )
  })
}


#' Input form for single-cell deconvolution
#' @export
#' @keywords internal
dtangleFormInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(id = ns('dtangle_form'), class = 'hidden-form', style = 'display: none;',
        selectizeInput(ns('dtangle_anal'), 'Reference single-cell dataset:', choices = '', width = '100%'),
        selectizeInputWithButtons(ns('include_clusters'),
                                  label = 'Clusters to include:',
                                  label_title = 'Select cell types that are expected in the bulk dataset',
                                  options = list(multiple = TRUE, placeholder = 'Select none to include all'),
                                  actionButton(ns('submit_dtangle'), '', icon = icon('chevron-right', 'fa-fw'), title = 'Submit cell-type deconvolution')
        )
    )
  })
}




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



