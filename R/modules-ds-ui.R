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
                div(id = ns('mds_plotly_container'), class = 'col-sm-7 mobile-margin',  style = 'display: none;',
                    dsPlotlyUI(ns('mds_plotly'))
                ),
                div(id = ns('gene_plotly_container'), class = 'col-sm-7 mobile-margin', style = 'display: none;',
                    dsPlotlyUI(ns('gene_plotly'))
                ),
                div(id = ns('cells_plotly_container'), class = 'col-sm-7 mobile-margin', style = '',
                    dsPlotlyUI(ns('cells_plotly'))
                )
            ),
            hr(),
            div(id = ns('quant_table_container'), style = 'display: none;',
                dsTable(ns('quant'))
            ),
            div(id = ns('anal_table_container'), style = 'display: none;',
                dsTable(ns('anal'))
            ),
            div(id = ns('explore_table_container'), style = 'display: none;',
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
          ns('anal_name'),
          'Analysis name:',
          downloadButton(ns('download'), label = NULL, icon = icon('download', 'fa-fw'), title = 'Download differential expression results'),
          options = list(create = TRUE, placeholder = 'Type name to create new analysis'),
          container_id = ns('anal_name_container'),
          help_id = ns('anal_name_help')
        ),
        div(id = ns('anal_buttons_panel'), style = 'display: none;',
            justifiedButtonGroup(
              container_id = ns('anal_labels'),
              label = 'Label selected rows:',
              help_id = ns('labels_help'),
              actionButton(ns('test'), 'test'),
              actionButton(ns('ctrl'), 'control'),
              actionButton(ns('reset'), 'reset')
            ),
            actionButton(ns('run_anal'), 'Run Analysis', width = '100%', class = 'btn-primary')
        )

    ),
    div(id = ns('explore_panel'), style = 'display: none;',

        selectizeInputWithButtons(ns('explore_genes'),
                                  'Show expression for:',
                                  actionButton(ns('show_dtangle'), '', icon = icon('object-ungroup', 'far fa-fw'), title = 'Toggle cell-type deconvolution'),
                                  options = list(maxItems = 6, multiple = TRUE)),
        div(class = 'hidden-forms',
            dtangleFormInput(ns('dtangle'))
        ),
        textInputWithButtons(ns('explore_group_name'),
                             container_id = ns('validate'),
                             'Group name for selected rows:',
                             actionButton(ns('grouped'), '', icon = icon('plus', 'fa-fw'), title = 'Add group'),
                             actionButton(ns('reset_explore'), '', icon = icon('redo-alt', 'fa-fw'), title = 'Clear all groups'),
                             help_id = ns('error_msg')
        )

    )
  )
}


#' Input form to control/test/all groups for integrated datasets
#' @export
#' @keywords internal
analysisTypeToggle <- function(ns) {

  shinyWidgets::radioGroupButtons(ns('anal_type'), "Select analysis type:",
                                  choices = c('differential expression', 'exploratory'),
                                  selected = 'differential expression', justified = TRUE)
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
