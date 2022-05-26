#' UI for Single Cell Exploration page
#'
#' @param id Identification string that is names-paced using \link[shiny]{NS}.
#' @param tab Name to appear on tab
#' @param active Name of current active \code{tab}
#'
#' @return shiny.tag with html for single-cell tab
#'
#' @export
#' @examples
#'
#' scPageUI('sc', tab = 'Single Cell', active = 'Single Cell')
#'
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
                scClusterPlotOutput(ns('cluster_plot'))
            )
        ),

        hr(),
        div(id = ns('comparison_row'),
            # row for cluster comparison
            div(class = 'row', id = ns('cluster_comparison_row'),
                div(class = "col-sm-12 col-lg-6 col-lg-push-6 mobile-margin",
                    scMarkerPlotOutput(ns('marker_plot_cluster'))
                ),
                div(class = "col-sm-12 col-lg-6 col-lg-pull-6 mobile-margin",
                    div(id = ns('biogps_container'), class="beside-downloadable",
                        scBioGpsPlotOutput(ns('biogps_plot'))
                    ),
                    div(style = 'display: none;', id = ns('violin_container'),
                        scViolinPlotOutput(ns('violin_plot'))
                    )
                )
            ),
            # row for samples comparison (integrated test vs ctrl)
            div(id = ns('sample_comparison_row'), class = 'invisible',
                div(class = 'row',
                    div(class = "col-sm-12 col-lg-6 mobile-margin", id = ns('col_left'),
                        scMarkerPlotOutput(ns('expr_test')),
                        br()
                    ),
                    div(class = "col-sm-12 col-lg-6 mobile-margin",
                        scMarkerPlotOutput(ns('expr_ctrl')),
                        br()
                    )
                ),
                div(class = 'row',
                    div(class = "col-sm-12 col-lg-6 mobile-margin",
                        scSamplePlotOutput(ns('expr_sample_violin'))
                    )
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
          label = 'Select a single-cell dataset:',
          actionButton(
            ns('show_label_resoln'), '',
            icon = icon('cog', 'fa-fw'),
            title = 'Toggle label transfer and cluster resolution',
            class = 'squashed-btn',
            `parent-style`='display: none;'
          ),
          actionButton(
            ns('show_subset'), '',
            icon = icon('project-diagram', 'fa-fw'),
            title = 'Toggle to subset current dataset or align to reference'
          )
        ),
        downloadLink(ns("download_dataset"), '')
    )
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
      actionButton(
        ns('overwrite_annot'),
        '',
        icon = tags$i(class = 'fa fa-plus fa-fw', tags$i(class='fa fa-ban fa-fw fa-hide')),
        title = 'Overwrite previous labels'
      ),
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
          sliderInput(
            ns('resoln'),
            label = HTML(paste0('Cluster resolution [n=<span id="', ns('nclus'),'">0</span>]:')),
            min=0.1, value=1, max=5.1, step = 0.1, width = '100%'),
          span(
            id = ns('provided_clusters_warning'),
            style='color: grey; font-style: italic; display: none;',
            tags$i(class = 'fas fa-exclamation-triangle text-warning'),
            ' using provided clusters.'
          )
      ),
      div(id = ns('resoln_ref_container'), style='display: none;',
          selectizeInput(
            ns('resoln_ref'),
            label = HTML(paste0('Cluster resolution [n=<span id="', ns('nclus_ref'),'">0</span>]:')), choices = '', width = '100%')
      )
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
          ns('subset_features'),
          container_id = ns('exclude-container'),
          label = 'Features to subset on:',
          actionButton(ns('toggle_exclude'), '', icon = tags$i(id =ns('toggle_icon'), class = 'fa fa-minus fa-fw text-warning'), title = 'Toggle to <span class="text-warning">exclude</span> or <span class="text-success">include</span> selected features'),
          options = list(multiple = TRUE,
                         optgroupField = 'type',
                         placeholder = 'Select none to recluster')),
        shiny::selectizeInput(ns('ref_name'), 'Reference:', choices = NULL, width = '100%', options = list(placeholder = 'optional')),
        shinypanel::textInputWithButtons(
          ns('subset_name'),
          container_id = ns('name-container'),
          label = 'Suffix for new dataset:',
          actionButton(ns('click_up'), '', icon = icon('upload', 'fa-fw'), title = 'Upload custom genes for clustering (optional)'),
          actionButton(ns('submit_subset'), '',
                       icon = tags$i(class = 'fa fa-plus fa-fw', tags$i(class='fa fa-ban fa-fw fa-hide')), title = 'Submit'),
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
                               icon = tags$i(class = 'fa fa-tag fa-fw', tags$i(class='fa fa-ban fa-fw fa-hide')),
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
                           icon = tags$i(class ='far fa-fw fa-window-close'))
            )
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
  biogps_btn  <- actionButton(ns('show_biogps'), label = NULL, icon = icon('chart-line', 'fa-fw'), title = 'Toggle BioGPS plot')
  pbulk_btn  <- actionButton(ns('show_pbulk'), label = NULL, icon = icon('chart-bar', 'fa-fw'), title = 'Toggle between \U0394 EXPRESSION and \U0394 CELLS layer')

  if (sample_comparison) {
    btn1 <- NULL
    btn2 <- pbulk_btn
  } else {
    btn1 <- biogps_btn
    btn2 <- custom_btn
  }

  div(id = 'sc-intro-feature',
      tags$label(class='control-label', `for`=ns('gene_table'), 'Select feature to plot:'),
      div(class = 'normal-header',
          DT::dataTableOutput(ns('gene_table'), width='100%', height='330px'),
      ),
      div(id=ns('gene_search_input'),
          style = 'height: 60px',
          shinypanel::textInputWithButtons(
            inputId = ns('gene_search'),
            label = '',
            btn1, btn2,
            placeholder = 'type regex to search gene'
          )
      ),
      div(id = ns('custom_metric_panel'), class = 'hidden-form', style = 'display: none',
          shinypanel::textAreaInputWithButtons(
            inputId = ns('custom_metric'),
            label = 'Custom metric:',
            placeholder = 'e.g: PF4\u{003E}2.2',
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
  div(class = 'invisible',
      id = ns('cluster_plot_container'),
      picker::pickerOutput(ns('cluster_plot'))
  )
}


#' Output plot of single cell abundances
#'
#' @keywords internal
#' @noRd
scAbundancePlotOutput <- function(id) {
  ns <- NS(id)

  div(style='margin-bottom:35px; display: none;',
      id = ns('abundance_plot_container'),
      shiny::plotOutput(ns('abundance_plot'))
  )
}


#' Output plot of single cell markers
#'
#' @keywords internal
#' @noRd
scMarkerPlotOutput <- function(id) {
  ns <- NS(id)
  div(class = 'invisible', id = ns('marker_plot_container'),
      picker::pickerOutput(ns('marker_plot'))
  )
}


#' Output plot of biogps data for a gene
#'
#' @keywords internal
#' @noRd
scBioGpsPlotOutput <- function(id) {
  ns <- NS(id)
  plotOutput(ns('biogps_plot'), height = '423px')
}


#' Output plots for samples comparison with integrated datasets
#'
#' @keywords internal
#' @noRd
scSamplePlotOutput <- function(id) {
  ns <- NS(id)
  div(style='margin-bottom: 20px',
      shiny::plotOutput(ns('plot'), height = 'auto')
  )
}


#' Output violin plot
#'
#' @keywords internal
#' @noRd
scViolinPlotOutput <- function(id) {
  ns <- NS(id)
  shiny::plotOutput(ns('violin_plot'), height = 'auto')
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
      label_title = 'ntest :: nctrl [<b>if N>2:</b> #p<0.05 <b>else:</b> #logFC>1]'
    ),

    # hidden dl/upload buttons
    downloadLink(ns('dl_anal'), '')
  )
}


scSampleGroupsInput <- function(id) {
  ns <- NS(id)

  tags$div(
    class = 'selectize-with-ht',
    shinypanel::selectizeInputWithButtons(
      inputId = ns('compare_groups'),
      label = 'Groups to compare:',
      actionButton(ns('edit_groups'), label = NULL, icon = tags$i(class ='far fa-fw fa-edit'), title = 'edit sample groups'),
      options = list(maxItems = 2, placeholder = 'Select test then control group'),
      container_id = ns('validate-up'),
      help_id = ns('error_msg')
    ),
    div(id=ns('groups_table_container'), class='handsontable-container',
        rhandsontable::rHandsontableOutput(ns('groups_table'), width='100%'),
        br(),
        span(class='pull-left', tags$i(class = 'fas fa-exclamation-triangle', style='color: orange;'), ' collapse to save changes.', style='color: grey; font-style: italic;'),
        hr()
    )
  )
}

