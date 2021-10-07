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
          label = 'Select a single-cell dataset:',
          actionButton(
            ns('show_label_resoln'), '',
            icon = icon('cog', 'fa-fw'),
            title = 'Toggle label transfer and cluster resolution',
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
          selectizeInput(
            ns('resoln_azi'),
            label = HTML(paste0('Cluster resolution [n=<span id="', ns('nclus_azi'),'">0</span>]:')), choices = '', width = '100%')
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
              choices = c('', unname(azimuth_refs)), width = '100%')
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
          options = list(multiple = TRUE,
                         optgroupField = 'type',
                         placeholder = 'Select none to recluster')),
        shiny::selectizeInput(ns('azimuth_ref'), 'Azimuth reference:', choices = c('', unname(azimuth_refs)), width = '100%', options = list(placeholder = 'optional')),
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
  pbulk_btn  <- actionButton(ns('show_pbulk'), label = NULL, icon = icon('chart-bar', 'fa-fw'), title = 'Toggle between Δ EXPRESSION and Δ CELLS layer')

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
      label_title = '(ntest :: nctrl **<b>hover for samples</b>**) [<b>if N>2:</b> #p<0.05 <b>else:</b> #logFC>1]'
    ),

    # hidden dl/upload buttons
    downloadLink(ns('dl_anal'), '')
  )
}




scSampleGroupsInput <- function(id) {
  ns <- NS(id)

  shinypanel::selectizeInputWithButtons(
    inputId = ns('compare_groups'),
    label = 'Groups to compare:',
    actionButton(ns('edit_groups'), label = NULL, icon = tags$i(class ='far fa-fw fa-edit'), title = 'edit sample groups and pairs (optional)'),
    options = list(maxItems = 2, placeholder = 'Select test then control group'),
    container_id = ns('validate-up'),
    help_id = ns('error_msg')

  )
}



tabs <- getShinyOption('tabs', c('Single Cell', 'Bulk Data', 'Drugs'))
data_dir <- getShinyOption('data_dir')
logout_url <- getShinyOption('logout_url')
is_example <- getShinyOption('is_example')
active <- tabs[1]

remoteDeps <- list()
if (!is.null(logout_url)) {

  dtCoreDeps <- htmltools::htmlDependency(
    'dt-core', '1.10.20',
    src = c(href = 'https://cdn.datatables.net/1.10.20/'),
    stylesheet = c('css/jquery.dataTables.min.css', 'css/jquery.dataTables.extra.css'),
    script = c('js/jquery.dataTables.min.js')
  )

  dtScrollerDeps <- htmltools::htmlDependency(
    'dt-ext-scroller', '1.10.20',
    src = c(href = 'https://cdn.datatables.net/scroller/2.0.1'),
    stylesheet = c('css/scroller.dataTables.min.css'),
    script = c('js/dataTables.scroller.min.js')
  )

  selectizeDep <- htmltools::htmlDependency(
    'selectize', '0.12.4',
    src = c(href = 'https://d174upwcmdw9dj.cloudfront.net/'),
    stylesheet = 'selectize.bootstrap3.min.css',
    script = c('selectize.min.js', 'selectize-plugin-a11y.js')
  )

  shinypanelDep <- htmltools::htmlDependency(
    'shinypanel', '0.1.5',
    src = c(href = 'https://d174upwcmdw9dj.cloudfront.net/'),
    stylesheet = c('shinypanel.css')
  )

  dseqrDep <- htmltools::htmlDependency(
    'dseqr', '0.20.27',
    src = c(href = 'https://d174upwcmdw9dj.cloudfront.net/'),
    stylesheet = c('custom.css')
  )

  htmlwidgetsDep <- htmltools::htmlDependency(
    'htmlwidgets', '1.5.3',
    src = c(href = "https://d174upwcmdw9dj.cloudfront.net/"),
    script = c("htmlwidgets.js")
  )
  pickerDep <- htmltools::htmlDependency(
    'picker', '0.2.5',
    src = c(href = "https://d174upwcmdw9dj.cloudfront.net/"),
    script = c("picker.min.js")
  )

  bootstrapDep <- htmltools::htmlDependency(
    "bootstrap", "3.4.1",
    src = c(href = "https://d174upwcmdw9dj.cloudfront.net/"),
    script = "bootstrap.min.js", stylesheet = "bootstrap.min.css"
  )

  fontawesomeDep <- htmltools::htmlDependency(
    'font-awesome', '5.13.0',
    src = c(href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.13.0/css/"),
    stylesheet = c("all.min.css", "v4-shims.min.css")
  )

  bootstrapSelectDep <- htmltools::htmlDependency(
    'bootstrap-select', '1.13.18',
    src = c(href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap-select/1.13.18/"),
    stylesheet = c("css/bootstrap-select.min.css"),
    script = c("js/bootstrap-select.min.js")
  )

  remoteDeps <- list(
    selectizeDep,
    bootstrapDep,
    fontawesomeDep,
    shinypanelDep,
    dseqrDep,
    htmlwidgetsDep,
    pickerDep,
    dtCoreDeps,
    dtScrollerDeps,
    bootstrapSelectDep
  )
}

bootstrapPage(
  remoteDeps,
  if (!is.null(logout_url)) includeHTML("www/gtm.html"),
  if (!is.null(logout_url)) tags$head(includeHTML("www/gtag.html")),
  useShinyjs(),
  rintrojs::introjsUI(),
  includeScript(path = 'www/renderSelectize.js'),
  includeScript(path = 'www/isMobile.js'),
  includeScript(path = 'www/contextMenu.js'),
  includeCSS(path = 'www/custom.css'),
  includeCSS(path = 'www/drugs.css'),
  includeCSS(path = 'www/pathways.css'),
  tags$head(HTML("<title>Dseqr</title>"),
            tags$link(rel = "icon", type = "image/png", href = "https://raw.githubusercontent.com/hms-dbmi/dseqr.sp/master/favicon.png")),


  navbarUI(tabs, active, logout_url),
  navbar2UI(is_example),

  fluidPage(
    tags$div(
      class = "tab-content shiny-bound-input", `data-tabsetid` = "tabset", id = "tab",

      # tabs
      scPageUI("sc", tab = 'Single Cell', active),
      bulkPageUI('bulk', tab = 'Bulk Data', active),
      drugsPageUI("drug", tab = 'Drugs', active),
      docsPageUI('docs', tab = 'Docs', active)
    )
  )

)




