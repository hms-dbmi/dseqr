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
        div(class = 'row row-fluid',
            div(class = 'col-sm-12 col-lg-6',
                style = 'align-self: start',
                scFormInput(ns('form'))
            ),
            div(class = 'col-sm-12 col-lg-6 mobile-margin',
                div(style='margin-bottom:35px;',
                    scAbundancePlotOutput(ns('abundance_plot'))
                ),
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
                div(class = "col-sm-12 col-lg-6 col-lg-pull-6 mobile-margin",
                    div(id = ns('biogps_container'), class="beside-downloadable",
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
            div(id = ns('label-resolution-form'), class = 'hidden-form', style = 'display: none;',
                labelTransferFormInput(ns('transfer')),
                hr(),
                resolutionFormInput(ns('resolution'))
            ),
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
                scSampleGroupsInput(ns('sample_groups')),

                div(id = ns('sample_cluster_input'),
                    scSampleClustersInput(ns('sample_clusters'), with_dl = TRUE),
                ),
                div(id = ns('sample_gene_input'),
                    selectedGeneInput(ns('gene_samples'), sample_comparison = TRUE)
                )
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
                                  choices = c('clusters', 'samples'),
                                  selected = 'clusters', justified = TRUE)
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
          label = 'Select or upload a single-cell dataset:',
          actionButton(
            ns('show_label_resoln'), '',
            icon = icon('cog', 'fa-fw'),
            title = 'Toggle for label transfer and cluster resolution',
            class = 'squashed-btn',
            `parent-style`='display: none;'
          ),
          actionButton(
            ns('show_integration'), '',
            icon = icon('object-ungroup', 'far fa-fw'),
            title = 'Toggle <b>once</b> to subset or <b>twice</b> to integrate datasets'
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
    shinypanel::selectizeInputWithButtons(
      ns('ref_name'), 'Transfer labels from:',
      actionButton(ns('overwrite_annot'), '', icon = icon('plus', 'fa-fw'), title = 'Overwrite previous labels'),
      options = list(optgroupField = 'type',
                     render = I('{option: transferLabelOption, item: scDatasetItemDF}'))
    )

  })
}


#' Input form for specifying leiden resolution parameter
#'
#' @keywords internal
#' @noRd
resolutionFormInput <- function(id) {
  ns <- NS(id)
  withTags({
    div(
      div(id=ns('resoln_container'),
          numericInput(
            ns('resoln'),
            label = HTML(paste0('Cluster resolution [n=<span id="', ns('nclus'),'">0</span>]:')),
            min=0.1, value=1, max=5.1, step = 0.1, width = '100%')
      ),
      div(id = ns('resoln_azi_container'), style='display: none;',
          numericInput(
            ns('resoln_azi'),
            label = HTML(paste0('Cluster resolution [n=<span id="', ns('nclus_azi'),'">0</span>]:')),
            min=1, value=2, max=3, step = 1, width = '100%')
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
        shinyWidgets::pickerInput(
          inputId = ns('integration_datasets'),
          label = 'Select datasets to integrate:',
          choices = '',
          width = '100%',
          options = shinyWidgets::pickerOptions(
            `selected-text-format` = "count > 0",
            actionsBox = TRUE,
            liveSearch = TRUE,
            size=14),
          multiple = TRUE),
        shinyWidgets::checkboxGroupButtons(ns('integration_types'), 'Integration types:', choices = c('harmony', 'fastMNN', 'Azimuth'), justified = TRUE, selected = 'harmony'),
        div(id=ns('azimuth_ref_container'), style='display: none;',
            selectizeInput(
              ns('azimuth_ref'),
              HTML('Select Azimuth reference:'),
              choices = c('', 'human_pbmc'), width = '100%')
        ),
        shinypanel::textInputWithButtons(
          ns('integration_name'),
          container_id = ns('name-container'),
          label = 'Name for integrated dataset:',
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
  gene_btn   <- actionButton(ns('genecards'), label = NULL, icon = icon('external-link-alt', 'fa-fw'), title = 'Go to GeneCards')
  custom_btn <- actionButton(ns('show_custom_metric'), label = NULL, icon = tags$i(class ='far fa-fw fa-edit'), title = 'Toggle custom metric')
  ridge_btn  <- actionButton(ns('show_ridge'), label = NULL, icon = icon('chart-line', 'fa-fw'), title = 'Toggle BioGPS plot')
  dprimes_btn  <- actionButton(ns('show_dprimes'), label = NULL, icon = icon('chart-bar', 'fa-fw'), title = 'Toggle Pseudobulk plot')

  if (sample_comparison) {
    btn1 <- NULL
    btn2 <- dprimes_btn
  } else {
    btn1 <- ridge_btn
    btn2 <- custom_btn
  }

  div(id = 'sc-intro-feature',
      tags$label(class='control-label', `for`=ns('gene_table'), 'Select feature to plot:'),
      DT::dataTableOutput(ns('gene_table'), width='100%', height='291px'),
      div(id=ns('gene_search_input'),
          style = 'height: 60px',
          shinypanel::textInputWithButtons(
            inputId = ns('gene_search'),
            label = '',
            btn1, btn2,
            placeholder = 'type regex to search'
          )
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


#' Output plot of single cell abundances
#'
#' @keywords internal
#' @noRd
scAbundancePlotOutput <- function(id) {
  ns <- NS(id)
  shiny::plotOutput(ns('abundance_plot'))
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


#' Input for Single Cell analysis
#'
#' Used in Single Cell and Drugs tab
#'
#' @keywords internal
#' @noRd
scSampleClustersInput <- function(id, with_dl = FALSE) {
  ns <- NS(id)

  dl_btn <- NULL
  if (with_dl)
    dl_btn <- actionButton(ns('click_dl_anal'), label = NULL, icon = icon('download', 'fa-fw'), title = 'Download results')

  tags$div(
    shinypanel::selectizeInputWithButtons(
      inputId = ns('selected_cluster'),
      label = 'Cluster for group comparison:',
      dl_btn,
      #TODO: implement logic for multi-cluster differential expression
      options = list(multiple = FALSE),
      label_title = '(ntest :: nctrl **<b>hover for samples</b>**) [<b>if N>2:</b> #p<0.05 <b>else:</b> #logFC>1]'
    ),

    # hidden dl/upload buttons
    downloadLink(ns('dl_anal'), ''),
    downloadLink(ns('dl_meta'), ''),
    div(style = 'display: none',
        fileInput(ns('up_meta'), '', accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
    )
  )

}


scSampleGroupsInput <- function(id) {
  ns <- NS(id)

  shinypanel::selectizeInputWithButtons(
    inputId = ns('compare_groups'),
    label = 'Groups to compare:',
    actionButton(ns('click_dl_meta'), label = NULL, icon = icon('download', 'fa-fw'), title = 'Download metadata to fill: <b>Group name</b> and <b>Pair</b> (optional)'),
    actionButton(ns('click_up_meta'), label = NULL, icon = icon('upload', 'fa-fw'), title = 'Upload filled metadata'),
    options = list(maxItems = 2, placeholder = 'Select test then control group'),
    container_id = ns('validate-up'),
    help_id = ns('error_msg')
  )
}
