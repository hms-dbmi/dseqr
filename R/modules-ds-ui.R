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
                div(class = 'col-sm-7 mobile-margin',
                    dsMDSplotlyUI(ns('mds_plotly'))
                )
            ),
            hr(),
            div(id = ns('quant_table_container'), style = 'display: none;',
                dsTable(ns('quant'))
            ),
            div(id = ns('anal_table_container'), style = 'display: none;',
                dsTable(ns('anal'))
            )
    )
  })
}

#' Plotly MDS output
#' @export
#' @keywords internal
dsMDSplotlyUI <- function(id) {
  ns <- NS(id)
  plotly::plotlyOutput(ns('plotly'))
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
    selectizeInputWithButtons(ns('dataset_name'), 'Dataset name:',
                              options = list(create = TRUE, placeholder = 'Type name to add new dataset', optgroupField = 'type'),
                              button(id = ns('dataset_dir'), type = 'button', class="shinyDirectories btn btn-default action-button shiny-bound-input disabled",
                                     `data-title` = 'Folder with fastq.gz or cellranger files',
                                     title = 'Select folder with fastq.gz or cellranger files',
                                     i(class = 'far fa-folder fa-fw')

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
    selectizeInputWithButtons(
      ns('anal_name'),
      'Analysis name:',
      downloadButton(ns('download'), label = NULL, icon = icon('download', 'fa-fw'), title = 'Download results'),
      options = list(create = TRUE, placeholder = 'Type name to create new analysis'),
      container_id = ns('anal_name_container'),
      help_id = ns('anal_name_help')
    ),
    div(id = ns('anal_buttons_panel'),
        justifiedButtonGroup(
          container_id = ns('anal_labels'),
          label = 'Label selected rows:',
          help_id = ns('labels_help'),
          actionButton(ns('test'), 'Test'),
          actionButton(ns('ctrl'), 'Control'),
          actionButton(ns('reset'), 'Reset')
        ),
        actionButton(ns('run_anal'), 'Run Analysis', width = '100%', class = 'btn-primary')
    )
  )
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
