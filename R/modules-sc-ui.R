#' UI for Single Cell Exploration page
#' @export
#' @keywords internal
scPageUI <- function(id, tab, active) {
  ns <- NS(id)
  active_class <- ifelse(tab == active, 'active', '')
  withTags({
    div(class = paste('tab-pane', active_class), `data-value` = tab, id = id_from_tab(tab),
        div(class = 'row',
            div(class = 'col-sm-12 col-lg-6',
                scFormInput(ns('form'))
            ),
            div(class = 'col-sm-12 col-lg-6 mobile-margin',
                scClusterPlotOutput(ns('cluster_plot'))
            )
        ),
        hr(),
        div(id = ns('comparison_row'), style = 'display: none;',
            # row for cluster comparison
            div(class = 'row', id = ns('cluster_comparison_row'), style = 'display: none;',
                div(class = "col-sm-12 col-lg-6 col-lg-push-6",
                    scMarkerPlotOutput(ns('marker_plot_cluster'))
                ),
                div(class = "col-sm-12 col-lg-6 col-lg-pull-6 mobile-margin",
                    div(id = ns('biogps_container'),
                        scBioGpsPlotOutput(ns('biogps_plot'))
                    ),
                    div( style = 'display: none;', id = ns('ridge_container'),
                         scRidgePlotOutput(ns('ridge_plot'))
                    )
                )
            ),
            # row for samples comparison (integrated test vs ctrl)
            div(class = 'row', id = ns('sample_comparison_row'), style = 'display: none;',
                div(class = "col-sm-12 col-lg-6 mobile-margin",
                    scSampleMarkerPlotOutput(ns('left'))
                ),
                div(class = "col-sm-12 col-lg-6 mobile-margin",
                    scSampleMarkerPlotOutput(ns('right'))
                )
            )
        ),
        # row for labels comparison (integration before and after)
        div(class = 'row', id = ns('label_comparison_row'), style = 'display: none;')
    )
  })
}



#' Input form for Single Cell Exploration page
#' @export
#' @keywords internal
scFormInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(class = "well-form well-bg",
        scSelectedDatasetInput(ns('dataset')),
        div(class = 'hidden-forms',
            labelTransferFormInput(ns('transfer')),
            integrationFormInput(ns('integration'))
        ),
        div(id = ns('form_container'), style = 'display: none;',
            div(style = 'display: none;', id = ns('comparison_toggle_container'), class = 'selectize-fh form-group',
                comparisonTypeToggle(ns('comparison'))
            ),
            # inputs for comparing clusters
            div(id = ns('cluster_comparison_inputs'),
                clusterComparisonInput(ns('cluster')),
                selectedGeneInput(ns('gene_clusters'))
            ),
            # inputs for comparing samples
            div(id = ns('sample_comparison_inputs'), style = 'display: none',
                scSampleComparisonInput(ns('sample'), with_dl = TRUE),
                selectedGeneInput(ns('gene_samples'), sample_comparison = TRUE)
            ),
            div(id = ns('label_comparison_inputs'), style = 'display: none;',
                selectedAnnotDatasetInput(ns('annot'))
            )
        )
    )
  })
}


#' Input form to control/test/all groups for integrated datasets
#' @export
#' @keywords internal
comparisonTypeToggle <- function(id) {
  ns <- NS(id)

  shinyWidgets::radioGroupButtons(ns('comparison_type'), "Perform comparisons between:",
                                  choices = c('clusters', 'samples'),
                                  selected = 'clusters', justified = TRUE)
}

#' Input for selecting datasets to show original labels for
#' @export
#' @keywords internal
selectedAnnotDatasetInput <- function(id) {
  ns <- NS(id)
  selectizeInput(ns('integration_anals'), 'Show original labels for:', multiple = TRUE, choices = '', width = '100%', options = list(maxItems = 2))
}

#' Input form/associated buttons for selecting single cell dataset
#' @export
#' @keywords internal
scSelectedDatasetInput <- function(id) {
  ns <- NS(id)

  tagList(
    selectizeInputWithButtons(ns('selected_dataset'), 'Select a single-cell dataset:',
                              actionButton(ns('show_label_transfer'), '', icon = icon('tag', 'fa-fw'), title = 'Toggle label transfer', class = 'squashed-btn'),
                              actionButton(ns('show_integration'), '',icon = icon('object-group', 'far fa-fw'), title = 'Toggle dataset integration'),
                              options = list(placeholder = 'Type name to add new single-cell dataset', create = TRUE)),
    shinyFiles::shinyDirLink(ns('new_dataset_dir'), '', 'Select folder with single cell fastq or cell ranger files')

  )

}


#' Input form for transfering labels between single cell datasets
#' @export
#' @keywords internal
labelTransferFormInput <- function(id) {
  ns <- NS(id)
  withTags({
    div(id = ns('label-transfer-form'), class = 'hidden-form', style = 'display: none;',
        selectizeInputWithButtons(ns('ref_name'), 'Transfer labels from:',
                                  actionButton(ns('overwrite_annot'), '', icon = icon('plus', 'fa-fw'), title = 'Overwrite previous labels'),
                                  options = list(optgroupField = 'type', placeholder = 'Select none to reset labels',
                                                 render = I('{option: transferLabelOption}'))
        )
    )
  })
}


#' Input form for integrating single cell datasets
#' @export
#' @keywords internal
integrationFormInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(id = ns('integration-form'), class = 'hidden-form', style = 'display: none;',
        selectizeInput(ns('test_integration'), 'Integration test datasets:', multiple = TRUE, choices = '', width = '100%'),
        selectizeInput(ns('ctrl_integration'), 'Integration control datasets:', multiple = TRUE, choices = '', width = '100%'),
        selectizeInput(ns('exclude_clusters'), 'Integration excluded clusters:', multiple = TRUE, choices = '', width = '100%', options = list(optgroupField = 'anal')),
        textInputWithButtons(ns('integration_name'),
                             container_id = ns('validate'),
                             'Name for new integrated analysis:',
                             actionButton(ns('submit_integration'), '', icon = icon('plus', 'fa-fw'), title = 'Integrate datasets'),
                             help_id = ns('error_msg'))
    )
  })
}



#' Input form and buttons to select a cluster or contrast and rename a cluster
#' @export
#' @keywords internal
clusterComparisonInput <- function(id) {
  ns <- NS(id)

  withTags({
    tagList(
      div(id = ns('selected_cluster_panel'),
          div(id = ns('select_panel'),

              selectizeInputWithButtons(ns('selected_cluster'),
                                        label = 'Show marker genes for:',
                                        label_title = 'Cluster (n cells :: % of total)',
                                        actionButton(ns('show_contrasts'), '',
                                                     icon = icon('chevron-right', 'fa-fw'),
                                                     title = 'Toggle single group comparisons'),
                                        actionButton(ns('show_rename'), '',
                                                     icon = icon('tag', 'fa-fw'),
                                                     title = 'Toggle rename cluster'))
          )
      ),
      div(id = ns('rename_panel'), style = 'display: none',
          textInputWithButtons(ns('new_cluster_name'),
                               'New cluster name:',
                               actionButton(ns('rename_cluster'), '',
                                            icon = icon('plus', 'fa-fw'),
                                            title = 'Rename cluster'))
      )
    )
  })
}


#' Input form to select gene for scBioGpsPlotOutput and scMarkerPlotOutput
#' @export
#' @keywords internal
selectedGeneInput <- function(id, sample_comparison = FALSE) {
  ns <- NS(id)

  exclude_ambient_button <- ridge_plot_button <- NULL
  if (sample_comparison)
    exclude_ambient_button <- actionButton(ns('exclude_ambient'), '',
                                           icon = icon('ban', 'fa-fw'),
                                           title = 'Toggle excluding ambient genes', class = 'squashed-btn')

  if (!sample_comparison)
    ridge_plot_button <- actionButton(ns('show_ridge'), label = NULL, icon = icon('chart-line', 'fa-fw'), title = 'Toggle ridgeline plot')


  selectizeInputWithButtons(id = ns('selected_gene'),
                            label = 'Show expression for:',
                            exclude_ambient_button,
                            ridge_plot_button,
                            actionButton(ns('genecards'), label = NULL, icon = icon('external-link-alt', 'fa-fw'), title = 'Go to GeneCards'),
                            label_title = 'Gene (% exp test :: % exp ctrl)'
  )
}


#' Output plot of single cell clusters
#' @export
#' @keywords internal
scClusterPlotOutput <- function(id) {
  ns <- NS(id)
  downloadablePlotUI(ns('cluster_plot'))
}

#' UI for plot with downloadable data
#' @export
#' @keywords internal
downloadablePlotUI <- function(id) {
  ns <- NS(id)
  withTags({
    div(class = 'downloadable-plot', style = 'display: none;', id = ns('plot_container'),
        div(class = 'clearfix',
            span(
              id = ns('download_container'), class = 'pull-right downloadable-plot-btn',
              downloadButton(ns('download'), label = NULL, icon = icon('download', 'fa-fw')),
              shinyBS::bsTooltip(
                ns('download'),
                'Download plot data',
                placement = 'left',
                options = list(
                  container = 'body',
                  template = '<div class="tooltip ggplot" role="tooltip"><div class="tooltip-arrow"></div><div class="tooltip-inner"></div></div>'
                ))
            )
        ),
        plotOutput(ns('dl_plot'))

    )
  })
}


#' Output plot of single cell markers
#' @export
#' @keywords internal
scMarkerPlotOutput <- function(id) {
  ns <- NS(id)
  downloadablePlotUI(ns('marker_plot'))
}

#' Output plot of biogps data for a gene
#' @export
#' @keywords internal
scBioGpsPlotOutput <- function(id) {
  ns <- NS(id)
  plotOutput(ns('biogps_plot'), height = '423px')
}

scSampleMarkerPlotOutput <- function(id) {
  ns <- NS(id)
  plotOutput(ns('plot'), height = 'auto')
}

scRidgePlotOutput <- function(id) {
  ns <- NS(id)
  plotOutput(ns('ridge_plot'), height = 'auto')
}

#' Input for Single Cell analysis
#'
#' Used in Single Cell, Drugs and Pathways tab
#'
#' @export
#' @keywords internal
scSampleComparisonInput <- function(id, with_dl = FALSE) {
  ns <- NS(id)

  dl_btn <- NULL
  if (with_dl)
    dl_btn <- downloadButton(ns('download'), label = NULL, icon = icon('download', 'fa-fw'), title = 'Download results')

  selectizeInputWithButtons(
    ns('selected_clusters'),
    label = 'Compare samples for:',
    actionButton(ns('run_comparison'), '', icon = icon('chevron-right', 'far fa-fw'), title = 'Compare test to control cells'),
    dl_btn,
    #TODO: implement logic for multi-cluster differential expression
    options = list(multiple = FALSE),
    label_title = 'Cluster (n test :: n ctrl)')

}

