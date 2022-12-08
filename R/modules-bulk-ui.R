#' UI for Bulk Data page
#'
#' @inheritParams scPageUI
#'
#' @return shiny.tag with html for bulk data tab
#'
#' @export
#' @examples
#'
#' bulkPageUI("bulk", tab = 'Bulk Data', active = 'Single Cell')
#'
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
             actionGroupButtonsWithClickHandlers(
               inputIds = c(ns('click_dl'), ns('click_up')),
               onclicks = c(sprintf('$("#%s").get(0).click()', ns('dl_annot')),
                            sprintf('$("#%s").click()', ns('up_annot'))),
               labels = list(icon('download', 'fa-fw'), icon('upload', 'fa-fw'))
             ),
             tags$div(class='help-block', id = ns('error_msg'))

    ),
    # hidden dl/upload buttons
    div(style = 'display: none;',
        fileInput(ns('up_annot'), '', accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
    ),
    downloadLink(ns('dl_annot'), ''),

    shinyBS::bsTooltip(id = ns('click_dl'), title = 'Download metadata to fill: <b>Group name</b> and <b>Pair</b> (optional)', options = list(container = 'body')),
    shinyBS::bsTooltip(id = ns('click_up'), title = 'Upload filled metadata', options = list(container = 'body'))
  )
}

actionGroupButtonsWithClickHandlers <- function(inputIds,
                                                labels,
                                                onclicks,
                                                status = "default",
                                                size = "normal",
                                                direction = "horizontal",
                                                fullwidth = FALSE) {
  stopifnot(length(inputIds) == length(labels))
  size <- match.arg(
    arg = size,
    choices = c("xs", "sm", "normal", "lg")
  )
  direction <- match.arg(
    arg = direction,
    choices = c("horizontal", "vertical")
  )
  status <- rep(status, length.out = length(labels))
  tags$div(
    class = "btn-group",
    class = if (direction == "vertical") "btn-group-vertical",
    class = paste0("btn-group-", size),
    class = if(fullwidth) "btn-group-justified",
    role = "group",
    lapply(
      X = seq_along(labels),
      FUN = function(i) {
        value <- shiny::restoreInput(id = inputIds[i], default = NULL)
        if (fullwidth) {
          tags$div(
            class = "btn-group", role = "group",
            class = paste0("btn-group-", size),
            tags$button(
              id = inputIds[i],
              type = "button",`data-val` = value,
              class = paste0("btn action-button btn-", status[i]),
              labels[i]
            )
          )
        } else {
          tags$button(
            id = inputIds[i],
            type = "button",`data-val` = value,
            class = paste0("btn action-button btn-", status[i]),
            onclick = onclicks[i],
            labels[i]
          )
        }
      }
    )
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
        div(class = 'hidden-forms',
            dtangleFormInput(ns('dtangle'))
        )
    )
  })
}


#' Differential expression analysis inputs for bulkFormInput
#'
#' @keywords internal
#' @noRd
bulkFormAnalInput <- function(id) {
  ns <- NS(id)

  tagList(
    bulkAnalInput(ns('ds'), label = 'Download two-group comparison:'),
    div(id = ns('explore_genes_container'),
        div(id = 'bulk-intro-feature',
            tags$label(class='control-label', `for`=ns('gene_table'), 'Select features to plot:'),
            div(class = 'normal-header',
                DT::dataTableOutput(ns('gene_table'), width='100%', height='330px'),
            ),
            div(id=ns('gene_search_input'),
                style = 'height: 60px',
                textInput(
                  inputId = ns('gene_search'),
                  width = '100%',
                  label = '',
                  placeholder = 'type regex to search gene'
                )
            )
        )
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

