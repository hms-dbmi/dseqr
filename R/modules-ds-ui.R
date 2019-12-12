#' UI for datasets page
#' @export
#' @keywords internal
dsPageUI <- function(id, tab, active) {
  ns <- NS(id)
  withTags({
    tabPane(tab, active,
            div(class = 'row',
                div(class = 'col-sm-5',
                    dsFormInput(ns('form'))
                ),
                div(id = ns('mds_plotly_container'), class = 'col-sm-7 mobile-margin', style = '',
                    dsPlotlyUI(ns('mds_plotly'))
                ),
                div(id = ns('gene_plotly_container'), class = 'col-sm-7 mobile-margin', style = 'display: none;',
                    dsPlotlyUI(ns('gene_plotly'))
                ),
                div(id = ns('cells_plotly_container'), class = 'col-sm-7 mobile-margin', style = 'display: none;',
                    dsPlotlyUI(ns('cells_plotly'))
                )
            ),
            hr(),
            div(id = ns('quant_table_container'), style = 'display: none;',
                dsTable(ns('quant'))
            ),
            div(id = ns('anal_table_container'), style = 'display: none;',
                dsTable(ns('explore'))
            )
    )
  })
}


#' Plotly MDS output
#' @export
#' @keywords internal
dsPlotlyUI <- function(id) {
  ns <- NS(id)
  downloadablePlotlyUI(ns('plotly'))
}

#' Input form for datasets page
#' @export
#' @keywords internal
dsFormInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(class = "well-form well-bg",
        dsDatasetInput(ns('selected_dataset')),
        div(id = ns('quant_dataset_panel'), style = 'display: none;',
            dsFormQuantInput(ns('quant_form'))
        ),
        div(id = ns('anal_dataset_panel'), style = 'display: none;',
            dsFormAnalInput(ns('anal_form'))
        )
    )
  })
}


#' Dataset selection input for dsFormInput
#' @export
#' @keywords internal
dsDatasetInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(
      selectizeInputWithButtons(ns('dataset_name'), 'Dataset name:',
                                options = list(create = TRUE, placeholder = 'Type name to add new dataset', optgroupField = 'type'),
                                button(id = ns('dataset_dir'), type = 'button', class="shinyDirectories btn btn-default action-button shiny-bound-input disabled",
                                       `data-title` = 'Folder with fastq.gz or cellranger files',
                                       title = 'Select folder with fastq.gz or cellranger files',
                                       i(class = 'far fa-folder fa-fw')

                                )

      )
    )
  })
}



#' Dataset quantification inputs for dsFormInput
#' @export
#' @keywords internal
dsFormQuantInput <- function(id) {
  ns <- NS(id)

  tagList(
    tags$div(id = ns('bulk_controls'),
             dsEndTypeInput(ns('end_type')),
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

#' Dataset end-type input for dsQuantInput
#' @export
#' @keywords internal
dsEndTypeInput <- function(id) {
  ns <- NS(id)

  selectizeInput(ns('end_type'),
                 'Confirm end-type:',
                 choices = NULL, width = '100%')
}

#' Differential expression analysis inputs for dsFormInput
#' @export
#' @keywords internal
dsFormAnalInput <- function(id) {
  ns <- NS(id)

  tagList(
    analysisTypeToggle(ns),
    div(id = ns('diff_panel'),
        selectizeInputWithButtons(
          ns('contrast_groups'),
          'Select test then control group:',
          downloadButton(ns('download'), label = NULL, icon = icon('download', 'fa-fw'), title = 'Download differential expression results'),
          actionButton(ns('run_anal'), label = NULL, icon = icon('chevron-right', 'fa-fw'), title = 'Run differential expression analysis'),
          options = list(maxItems = 2, placeholder = 'Select test group'),
          container_id = ns('run_anal_container'),
          help_id = ns('run_anal_help')
        )
    ),
    div(id = ns('explore_panel'), style = 'display: none;',

        selectizeInputWithButtons(ns('explore_genes'),
                                  'Show expression for genes:',
                                  actionButton(ns('show_dtangle'), '', icon = icon('object-ungroup', 'far fa-fw'), title = 'Toggle cell-type deconvolution'),
                                  svaButton(inputId = ns('show_nsv'), sliderId = ns('num_svs')),
                                  options = list(maxItems = 6, multiple = TRUE)),
        div(class = 'hidden-forms',
            dtangleFormInput(ns('dtangle'))
        ),

    ),
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
    title = 'Number of surrogate variables adjusting for'
  )
}


#' Input form to control/test/all groups for integrated datasets
#' @export
#' @keywords internal
analysisTypeToggle <- function(ns) {

  shinyWidgets::radioGroupButtons(ns('anal_type'), "Select analysis type:",
                                  choices = c('exploratory', 'differential expression'),
                                  selected = 'exploratory', justified = TRUE)
}

#' Tables for datasets page
#' @export
#' @keywords internal
dsTable <- function(id) {
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

