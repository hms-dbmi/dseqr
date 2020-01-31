#' UI for Bulk Data page
#' @export
#' @keywords internal
bulkPageUI <- function(id, tab, active) {
  ns <- NS(id)
  withTags({
    tabPane(tab, active,
            div(class = 'row',
                div(class = 'col-lg-4',
                    bulkFormInput(ns('form'))
                ),
                div(id = ns('mds_plotly_container'),
                    div(class = 'col-lg-4 mobile-margin',
                        bulkPlotlyUI(ns('mds_plotly_unadjusted'))
                    ),
                    div(class = 'col-lg-4 mobile-margin',
                        bulkPlotlyUI(ns('mds_plotly_adjusted'))
                    )
                ),
                div(id = ns('gene_plotly_container'), class = 'col-lg-7 mobile-margin', style = 'display: none;',
                    bulkPlotlyUI(ns('gene_plotly'))
                ),
                div(id = ns('cells_plotly_container'), class = 'col-lg-7 mobile-margin', style = 'display: none;',
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
bulkPlotlyUI <- function(id, height = 'auto') {
  ns <- NS(id)
  downloadablePlotlyUI(ns('plotly'), height)
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
      selectizeInputWithButtons(
        ns('dataset_name'), 'Select a bulk dataset:',
        container_id = 'dataset_name_container',
        options = list(create = TRUE, placeholder = 'Type name to add new bulk dataset', optgroupField = 'type'),
        svaButton(inputId = ns('show_nsv'), sliderId = ns('selected_nsv')),
        actionButton(ns('show_dtangle'), '', icon = icon('object-ungroup', 'far fa-fw'), title = 'Toggle cell-type deconvolution'),
        hide_btns = TRUE
      ),
      shinyFiles::shinyDirLink(ns('new_dataset_dir'), '', 'Select folder fastq.gz files'),
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
             selectizeInput(ns('end_type'),
                            'Confirm end-type:',
                            choices = NULL, width = '100%'),
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


#' Differential expression analysis inputs for bulkFormInput
#' @export
#' @keywords internal
bulkFormAnalInput <- function(id) {
  ns <- NS(id)

  tagList(
    bulkAnalInput(ns('ds'), label = 'Download two-group comparison:'),
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
        selectizeInput(ns('dtangle_dataset'), 'Reference single-cell dataset:', choices = '', width = '100%'),
        selectizeInputWithButtons(ns('include_clusters'),
                                  label = 'Clusters to include:',
                                  label_title = 'Select cell types that are expected in the bulk dataset',
                                  options = list(multiple = TRUE, placeholder = 'Select none to include all'),
                                  actionButton(ns('submit_dtangle'), '', icon = icon('chevron-right', 'fa-fw'), title = 'Submit cell-type deconvolution')
        )
    )
  })
}


#' Bulk Differential expression analysis input
#' @export
#' @keywords internal
bulkAnalInput <- function(id, with_dl = TRUE, label = 'Select groups to compare:') {
  ns <- NS(id)

  options <- list(maxItems = 2, placeholder = 'Select test then control group')
  if (with_dl) {
    input <- selectizeInputWithButtons(
      ns('contrast_groups'),
      label,
      downloadButton(ns('download'), label = NULL, icon = icon('download', 'fa-fw'), title = 'Download differential expression analysis'),
      options = options,
      container_id = ns('run_anal_container')
    )

  } else {
    input <- tags$div(
      id = ns('run_anal_container'),
      selectizeInputWithValidation(
        ns('contrast_groups'),
        label,
        options = options
      )

    )
  }

  return(input)
}


