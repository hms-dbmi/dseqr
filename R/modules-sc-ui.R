#' UI for Single Cell Exploration page
#'
#' @param id Identification string that is names-paced using \link[shiny]{NS}.
#' @param tab Name to appear on tab
#' @param active Name of current active \code{tab}
#'
#' @return shiny.tag with html for single-cell tab
#'
#' @export
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
        div(id = ns('comparison_row'), style = '',
            # row for cluster comparison
            div(class = 'row', id = ns('cluster_comparison_row'), style = 'display: none;',
                div(class = "col-sm-12 col-lg-6 col-lg-push-6",
                    scMarkerPlotOutput(ns('marker_plot_cluster'))
                ),
                div(class = "col-sm-12 col-lg-6 col-lg-pull-6 mobile-margin beside-downloadable",
                    div(id = ns('biogps_container'),
                        scBioGpsPlotOutput(ns('biogps_plot'))
                    ),
                    div(style = 'display: none;', id = ns('ridge_container'),
                        scRidgePlotOutput(ns('ridge_plot'))
                    )
                )
            ),
            # row for samples comparison (integrated test vs ctrl)
            div(class = 'row', id = ns('sample_comparison_row'), style = 'display: none;',
                div(class = "col-sm-12 col-lg-6 mobile-margin", id = ns('col_left'),
                    scSampleMarkerPlotOutput(ns('left'))
                ),
                div(class = "col-sm-12 col-lg-6 mobile-margin",
                    scSampleMarkerPlotOutput(ns('right')),
                    br(),
                    scSampleMarkerPlotOutput(ns('right_bottom'))
                )
            ),
            div(class = 'row', id = ns('labels_comparison_row'), style = 'display: none;',
                div(class = "col-sm-12 col-lg-7",
                    scLabelsPlotOutput(ns('labels_plot_cluster'))
                )
            )
        )
    )
  })
}



#' Input form for Single Cell Exploration page
#'
#' @keywords internal
#' @noRd
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
#'
#' @keywords internal
#' @noRd
comparisonTypeToggle <- function(id) {
  ns <- NS(id)

  shinyWidgets::radioGroupButtons(ns('comparison_type'), "Perform comparisons between:",
                                  choices = c('clusters', 'samples', 'labels'),
                                  selected = 'clusters', justified = TRUE)
}

#' Input for selecting datasets to show original labels for
#'
#' @keywords internal
#' @noRd
selectedAnnotDatasetInput <- function(id) {
  ns <- NS(id)
  selectizeInput(ns('integration_anals'), 'Show original labels for:', multiple = TRUE, choices = '', width = '100%', options = list(maxItems = 2))
}

#' Input form/associated buttons for selecting single cell dataset
#'
#' @keywords internal
#' @noRd
scSelectedDatasetInput <- function(id) {
  ns <- NS(id)

  tagList(
    div(id = 'sc-intro-dataset',
        shinypanel::selectizeInputWithButtons(
          inputId = ns('selected_dataset'),
          label = 'Select a single-cell dataset:',
          actionButton(
            ns('show_label_transfer'), '',
            icon = icon('tag', 'fa-fw'),
            title = 'Toggle label transfer',
            class = 'squashed-btn',
            `parent-style`='display: none;'
          ),
          actionButton(
            ns('show_integration'), '',
            icon = icon('object-ungroup', 'far fa-fw'),
            title = 'Toggle <b>once</b> to subset or <b>twice</b> to integrate dataset(s)'
          )
        )
    ),

    shinyFiles::shinyDirLink(ns('new_dataset_dir'), '', 'Select folder with single cell fastq or cell ranger files')
  )

}


#' Input form for transferring labels between single cell datasets
#'
#' @keywords internal
#' @noRd
labelTransferFormInput <- function(id) {
  ns <- NS(id)
  withTags({
    div(id = ns('label-transfer-form'), class = 'hidden-form', style = 'display: none;',
        shinypanel::selectizeInputWithButtons(
          ns('ref_name'), 'Transfer labels from:',
          actionButton(ns('overwrite_annot'), '', icon = icon('plus', 'fa-fw'), title = 'Overwrite previous labels'),
          options = list(optgroupField = 'type',
                         render = I('{option: transferLabelOption}'))
        )
    )
  })
}


#' Input form for integrating single cell datasets
#'
#' @keywords internal
#' @noRd
integrationFormInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(id = ns('integration-form'), class = 'hidden-form', style = 'display: none;',
        selectizeInput(ns('test_integration'), 'Integration test datasets:', multiple = TRUE, choices = '', width = '100%'),
        selectizeInput(ns('ctrl_integration'), 'Integration control datasets:', multiple = TRUE, choices = '', width = '100%'),
        shinypanel::selectizeInputWithButtons(
          ns('subset_clusters'),
          container_id = ns('exclude-container'),
          label = 'Clusters to subset on:',
          actionButton(ns('toggle_exclude'), '', icon = tags$i(id =ns('toggle_icon'), class = 'fa fa-minus fa-fw text-warning'), title = 'Toggle to <span class="text-warning">exclude</span> or <span class="text-success">include</span> selected clusters'),
          options = list(multiple = TRUE, optgroupField = 'anal', placeholder = 'select none to keep all clusters')),
        shinyWidgets::checkboxGroupButtons(ns('integration_types'), 'Integration types:', choices = c('harmony', 'fastMNN'), justified = TRUE, selected = 'harmony'),
        shinypanel::textInputWithButtons(
          ns('integration_name'),
          container_id = ns('name-container'),
          label = 'Name for integrated dataset:',
          actionButton(ns('click_dl'), '', icon = icon('download', 'fa-fw'), title = 'Download sample pairs csv to fill out'),
          actionButton(ns('click_up'), '', icon = icon('upload', 'fa-fw'), title = 'Upload filled in sample pairs csv'),
          actionButton(ns('submit_integration'), '', icon = icon('plus', 'fa-fw'), title = 'Integrate datasets'),
          help_id = ns('error_msg'),
          placeholder = 'Prepended to integration type(s)'),

        div(style = 'display: none',
            fileInput(ns('up_pairs'), '', accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
        ),
        downloadLink(ns('dl_samples'), '')
    )
  })
}

#' Input form for subsetting single cell datasets
#'
#' @keywords internal
#' @noRd
subsetFormInput <- function(id) {

  ns <- NS(id)

  withTags({
    div(id = ns('subset-form'), class = 'hidden-form', style = 'display: none;',
        shinypanel::selectizeInputWithButtons(
          ns('subset_clusters'),
          container_id = ns('exclude-container'),
          label = 'Features to subset on:',
          actionButton(ns('toggle_exclude'), '', icon = tags$i(id =ns('toggle_icon'), class = 'fa fa-minus fa-fw text-warning'), title = 'Toggle to <span class="text-warning">exclude</span> or <span class="text-success">include</span> selected features'),
          options = list(multiple = TRUE, optgroupField = 'type')),
        shinypanel::textInputWithButtons(
          ns('subset_name'),
          container_id = ns('name-container'),
          label = 'Name for subset dataset:',
          actionButton(ns('click_up'), '', icon = icon('upload', 'fa-fw'), title = 'Upload custom genes for clustering'),
          actionButton(ns('submit_subset'), '', icon = icon('plus', 'fa-fw'), title = 'Subset dataset'),
          help_id = ns('error_msg'),
          placeholder = 'eg: QC2 (appended to founder name)'),

        div(style = 'display: none',
            fileInput(ns('up_hvgs'), '', accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
        )
    )
  })

}


#' Input form and buttons to select a cluster or contrast and rename a cluster
#'
#' @keywords internal
#' @noRd
clusterComparisonInput <- function(id) {
  ns <- NS(id)

  withTags({
    div(id = 'sc-intro-cluster',
        div(id = ns('selected_cluster_panel'),
            div(id = ns('select_panel'),
                shinypanel::selectizeInputWithButtons(
                  ns('selected_cluster'),
                  label = 'Sort marker genes for:',
                  label_title = 'Cluster (n cells :: % of total)',
                  actionButton(ns('show_rename'), '',
                               icon = icon('tag', 'fa-fw'),
                               title = 'Toggle rename cluster'
                  ),
                  actionButton(ns('show_contrasts'), '',
                               icon = icon('chevron-right', 'fa-fw'),
                               title = 'Toggle single group comparisons'
                  )
                )
            )
        ),
        div(id = ns('rename_panel'), style = 'display: none',
            shinypanel::textInputWithButtons(
              ns('new_cluster_name'),
              'New cluster name:',
              actionButton(ns('rename_cluster'), '',
                           icon = icon('plus', 'fa-fw'),
                           title = 'Rename cluster'))
        )

    )
  })
}


#' Input form to select gene for scBioGpsPlotOutput and scMarkerPlotOutput
#'
#' @keywords internal
#' @noRd
selectedGeneInput <- function(id, sample_comparison = FALSE) {
  ns <- NS(id)
  btn1 <- btn2 <- btn3 <- NULL
  gene_btn   <- actionButton(ns('genecards'), label = NULL, icon = icon('external-link-alt', 'fa-fw'), title = 'Go to GeneCards', `parent-style`='display: none;')
  custom_btn <- actionButton(ns('show_custom_metric'), label = NULL, icon = tags$i(class ='far fa-fw fa-edit'), title = 'Toggle custom metric')
  amb_btn    <- actionButton(ns('exclude_ambient'), label = NULL, icon = icon('ban', 'fa-fw'), title = 'Toggle excluding ambient genes')
  ridge_btn  <- actionButton(ns('show_ridge'), label = NULL, icon = icon('chart-line', 'fa-fw'), title = 'Toggle BioGPS plot',  `parent-style`='display: none;')

  if (sample_comparison) {
    btn1 <- gene_btn
    btn2 <- amb_btn
  } else {
    btn1 <- ridge_btn
    btn2 <- gene_btn
    btn3 <-custom_btn
  }

  div(id = 'sc-intro-feature',
      shinypanel::selectizeInputWithButtons(
        inputId = ns('selected_gene'),
        label = 'Feature to plot:',
        btn1, btn2, btn3,
        options = list(optgroupField = 'type')
      ),
      div(id = ns('custom_metric_panel'), class = 'hidden-form', style = 'display: none',
          shinypanel::textAreaInputWithButtons(
            inputId = ns('custom_metric'),
            label = 'Custom metric:',
            placeholder = 'e.g: PF4>2.2',
            actionButton(ns('save_custom_metric'), '',
                         icon = icon('plus', 'fa-fw'),
                         title = 'Save custom metric'))
      )

  )
}


#' Output plot of single cell clusters
#'
#' @keywords internal
#' @noRd
scClusterPlotOutput <- function(id) {
  ns <- NS(id)
  shinydlplot::downloadablePlotUI(ns('cluster_plot'))
}


#' Output plot of single cell markers
#'
#' @keywords internal
#' @noRd
scMarkerPlotOutput <- function(id) {
  ns <- NS(id)
  shinydlplot::downloadablePlotUI(ns('marker_plot'))
}

#' Output plot of biogps data for a gene
#'
#' @keywords internal
#' @noRd
scBioGpsPlotOutput <- function(id) {
  ns <- NS(id)
  plotOutput(ns('biogps_plot'), height = '423px')
}

#' Output plot/plotly for samples comparison with integrated datasets
#'
#' @keywords internal
#' @noRd
scSampleMarkerPlotOutput <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(shinydlplot::downloadablePlotUI(ns('plot'), height = 'auto'), style = 'line-height: 0px;'),
    plotly::plotlyOutput(ns('plotly'), height = 'auto')
  )
}

#' Output Ridgeline plot
#'
#' @keywords internal
#' @noRd
scRidgePlotOutput <- function(id) {
  ns <- NS(id)
  shinydlplot::downloadablePlotUI(ns('ridge_plot'), height = 'auto')
  # plotOutput(ns('ridge_plot'), height = 'auto')
}

#' Output plotly for labels comparison with integrated datasets
#'
#' @keywords internal
#' @noRd
scLabelsPlotOutput <- function(id) {
  ns <- NS(id)
  plotly::plotlyOutput(ns('labels_plot'))
}

#' Input for Single Cell analysis
#'
#' Used in Single Cell and Drugs tab
#'
#' @keywords internal
#' @noRd
scSampleComparisonInput <- function(id, with_dl = FALSE) {
  ns <- NS(id)

  dl_btn <- NULL
  if (with_dl)
    dl_btn <- downloadButton(ns('download'), label = NULL, icon = icon('download', 'fa-fw'), title = 'Download results')

  shinypanel::selectizeInputWithButtons(
    inputId = ns('selected_cluster'),
    label = 'Compare samples for:',
    dl_btn,
    #TODO: implement logic for multi-cluster differential expression
    options = list(multiple = FALSE),
    label_title = '(ntest :: nctrl **<b>hover for samples</b>**) [<b>if reps:</b> #p<0.05 <b>else:</b> #logFC>1]')

}

