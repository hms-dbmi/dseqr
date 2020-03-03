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
            ),
            div(class = 'row', id = ns('labels_comparison_row'), style = 'display: none;',
                div(class = "col-sm-12 col-lg-7",
                    scLabelsPlotOutput(ns('labels_plot_cluster'))
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
            integrationFormInput(ns('integration')),
            subsetFormInput(ns('subset'))
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
            div(id = ns('labels_comparison_inputs'), style = 'display: none;',
                scSampleComparisonInput(ns('labels'))
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
                                  choices = c('clusters', 'samples', 'labels'),
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
                              actionButton(ns('show_label_transfer'), '', icon = icon('tag', 'fa-fw'), title = 'Toggle label transfer', class = 'squashed-btn', `parent-style`='display: none;'),
                              actionButton(ns('show_integration'), '',icon = icon('object-group', 'far fa-fw'), title = 'Toggle <b>once</b> to integrate or <b>twice</b> to subset dataset(s)'),
                              options = list(placeholder = 'Type name to add new single-cell dataset', optgroupField = 'type', create = TRUE)),
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
                                  options = list(optgroupField = 'type',
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
        selectizeInputWithButtons(ns('exclude_clusters'),
                                  container_id = ns('exclude-container'),
                                  label = tags$span('Clusters to', tags$span(class="text-warning", 'exclude'), 'or', tags$span(class='text-success', 'include', .noWS = 'after'), ':'),
                                  actionButton(ns('toggle_exclude'), '', icon = tags$i(id =ns('toggle_icon'), class = 'fa fa-minus fa-fw text-warning'), title = 'Toggle exclude or include'),
                                  options = list(multiple = TRUE, optgroupField = 'anal')),
        shinyWidgets::radioGroupButtons(ns('integration_type'), 'Integration type:', choices = c('harmony', 'liger', 'fastMNN'), justified = TRUE, selected = 'harmony'),
        textInputWithButtons(ns('integration_name'),
                             container_id = ns('name-container'),
                             label = 'Name for new dataset:',
                             actionButton(ns('click_dl'), '', icon = icon('download', 'fa-fw'), title = 'Download sample pairs csv to fill out'),
                             actionButton(ns('click_up'), '', icon = icon('upload', 'fa-fw'), title = 'Upload filled in sample pairs csv'),
                             actionButton(ns('submit_integration'), '', icon = icon('plus', 'fa-fw'), title = 'Integrate datasets'),
                             help_id = ns('error_msg')),

        div(style = 'display: none',
            fileInput(ns('up_pairs'), '', accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
        ),
        downloadLink(ns('dl_samples'), '')
    )
  })
}

subsetFormInput <- function(id) {

  ns <- NS(id)

  withTags({
    div(id = ns('subset-form'), class = 'hidden-form', style = 'display: none;',
        selectizeInputWithButtons(ns('exclude_clusters'),
                                  container_id = ns('exclude-container'),
                                  label = tags$span('Clusters to', tags$span(class="text-warning", 'exclude'), 'or', tags$span(class='text-success', 'include', .noWS = 'after'), ':'),
                                  actionButton(ns('toggle_exclude'), '', icon = tags$i(id =ns('toggle_icon'), class = 'fa fa-minus fa-fw text-warning'), title = 'Toggle exclude or include'),
                                  options = list(multiple = TRUE, optgroupField = 'anal')),
        textInputWithButtons(ns('subset_name'),
                             container_id = ns('name-container'),
                             label = 'Name for new dataset:',
                             actionButton(ns('submit_subset'), '', icon = icon('plus', 'fa-fw'), title = 'Subset dataset'),
                             help_id = ns('error_msg'),
                             placeholder = 'eg: QC2 (appended to founder dataset name)'),

        div(style = 'display: none',
            fileInput(ns('up_pairs'), '', accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
        ),
        downloadLink(ns('dl_samples'), '')
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
                                        actionButton(ns('show_rename'), '',
                                                     icon = icon('tag', 'fa-fw'),
                                                     title = 'Toggle rename cluster'),
                                        actionButton(ns('show_contrasts'), '',
                                                     icon = icon('chevron-right', 'fa-fw'),
                                                     title = 'Toggle single group comparisons')
              )
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

  exclude_ambient_button <- ridge_plot_button <- custom_button <- NULL
  if (sample_comparison)
    exclude_ambient_button <- actionButton(ns('exclude_ambient'), '',
                                           icon = icon('ban', 'fa-fw'),
                                           title = 'Toggle excluding ambient genes', class = 'squashed-btn')

  if (!sample_comparison) {
    ridge_plot_button <- actionButton(ns('show_ridge'), label = NULL, icon = icon('chart-line', 'fa-fw'), title = 'Toggle BioGPS plot',  `parent-style`='display: none;')
    custom_button <- actionButton(ns('show_custom_metric'), label = NULL, icon = tags$i(class ='far fa-fw fa-edit'), title = 'Toggle custom metric')
  }



  tagList(
    selectizeInputWithButtons(id = ns('selected_gene'),
                              label = 'Feature to plot:',
                              exclude_ambient_button,
                              ridge_plot_button,
                              options = list(optgroupField = 'type'),
                              actionButton(ns('genecards'), label = NULL, icon = icon('external-link-alt', 'fa-fw'), title = 'Go to GeneCards', `parent-style`='display: none;'),
                              custom_button
    ),
    div(id = ns('custom_metric_panel'), class = 'hidden-form', style = 'display: none',
        textInputWithButtons(ns('custom_metric'),
                             'Custom metric:',
                             placeholder = 'e.g: PF4 >= 2.2',
                             actionButton(ns('update_custom_metric'), '',
                                          icon = icon('redo', 'fa-fw'),
                                          title = 'Reload custom metric'),
                             actionButton(ns('save_custom_metric'), '',
                                          icon = icon('plus', 'fa-fw'),
                                          title = 'Save custom metric'))
    )

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
  tagList(
    tags$div(plotOutput(ns('plot'), height = 'auto'), style = 'line-height: 0px;'),
    plotly::plotlyOutput(ns('plotly'), height = 'auto')
  )
}


scRidgePlotOutput <- function(id) {
  ns <- NS(id)
  plotOutput(ns('ridge_plot'), height = 'auto')
}

scLabelsPlotOutput <- function(id) {
  ns <- NS(id)
  plotly::plotlyOutput(ns('labels_plot'))
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
    ns('selected_cluster'),
    label = 'Compare samples for:',
    dl_btn,
    #TODO: implement logic for multi-cluster differential expression
    options = list(multiple = FALSE),
    label_title = '(ntest :: nctrl **<b>hover for samples</b>**) [<b>if reps:</b> #p<0.05 <b>else:</b> #logFC>1]')

}


