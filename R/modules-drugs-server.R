#' Logic for Drugs page
#' @export
#' @keywords internal
drugsPage <- function(input, output, session) {

  # the form area inputs/results
  form <- callModule(drugsForm, 'form')

  # the output table
  callModule(drugsTable, 'table',
             cmap_res = form$cmap_res,
             l1000_res = form$l1000_res,
             drug_study = form$drug_study,
             cells = form$cells,
             show_clinical = form$show_clinical)




}

# Logic for form on drugs page
#' @export
#' @keywords internal
drugsForm <- function(input, output, session) {

  drugStudy <- callModule(selectedDrugStudy, 'drug_study')

  advancedOptions <- callModule(advancedOptions, 'advanced',
                                drugStudy$cmap_res,
                                drugStudy$l1000_res,
                                drugStudy$drug_study,
                                drugStudy$show_advanced)




  return(list(
    cmap_res = drugStudy$cmap_res,
    l1000_res = drugStudy$l1000_res,
    drug_study = drugStudy$drug_study,
    cells = advancedOptions$cells,
    show_clinical = drugStudy$show_clinical
  ))


}

#' Logic for selected drug study
#' @export
#' @keywords internal
selectedDrugStudy <- function(input, output, session) {

  data_dir <- '/home/alex/Documents/Batcave/zaklab/drugseqr/data-raw/bulk/example-data/IBD'

  cmap_res <- readRDS(file.path(data_dir, 'cmap_res.rds'))
  l1000_res <- readRDS(file.path(data_dir, 'l1000_res.rds'))

  drug_study <- reactive(input$study)

  # boolean for advanced options
  show_advanced <- reactive({
    input$advanced %% 2 != 0
  })


  study_choices <- c('CMAP02', 'L1000')
  updateSelectizeInput(session, 'study', choices = study_choices, selected = 'CMAP02')


  # toggle for clinical status
  show_clinical <- reactive({
    input$clinical %% 2 != 0
  })
  observe({
    toggleClass('advanced', 'btn-primary', condition = show_advanced())
    toggleClass('clinical', 'btn-primary', condition = show_clinical())
  })



  return(list(
    cmap_res = cmap_res,
    l1000_res = l1000_res,
    drug_study = reactive(input$study),
    show_clinical = show_clinical,
    show_advanced = show_advanced
  ))

}

#' Logic for advanced options for selectedDrugStudy
#' @export
#' @keywords internal
advancedOptions <- function(input, output, session, cmap_res, l1000_res, drug_study, show_advanced) {


  null_cmap <- is.null(cmap_res)
  null_l1000 <- is.null(l1000_res)


  # available cell lines
  if (!null_cmap)
    cmap_cells  <- unique(gsub('^[^_]+_([^_]+)_.+?$', '\\1', names(cmap_res)))

  if (!null_l1000)
    l1000_cells <- unique(gsub('^[^_]+_([^_]+)_.+?$', '\\1', names(l1000_res)))

  #  toggle button styling and showing advanced options
  shiny::observe({
    toggle('advanced-panel', condition = show_advanced(), anim = TRUE)
  })


  # update choices for cell lines
  shiny::observe({
    req(drug_study())
    if (drug_study() == 'L1000') {
      cell_choices <- l1000_cells

    } else if (drug_study() == 'CMAP02') {
      cell_choices <- cmap_cells
    }
    shiny::updateSelectizeInput(session, 'cells', choices = cell_choices, selected = NULL)
  })

  return(list(
    cells = reactive(input$cells)
  ))

}

#' Logic for drug table
#' @export
#' @keywords internal
drugsTable <- function(input, output, session, cmap_res, l1000_res, drug_study, cells, show_clinical) {

  # generate table to display
  query_res <- shiny::reactive({
    drug_study <- drug_study()

    req(drug_study)
    if (drug_study == 'L1000') {
      query_res <- study_table(l1000_res, 'L1000', cells())

    } else if (drug_study == 'CMAP02') {
      query_res <- study_table(cmap_res, 'CMAP02', cells())
    }

    # for removing entries without a clinical phase
    if (show_clinical()) {
      query_res <- tibble::as_tibble(query_res)
      query_res <- dplyr::filter(query_res, !is.na(.data$`Clinical Phase`))
    }

    return(query_res)
  })

  # show query data
  output$query_res <- DT::renderDataTable({

    query_res <- query_res()
    wide_cols <- c('MOA', 'Target', 'Disease Area', 'Indication', 'Vendor', 'Catalog #', 'Vendor Name')
    # -1 needed with rownames = FALSE
    elipsis_targets <- which(colnames(query_res) %in% wide_cols) - 1

    DT::datatable(
      query_res,
      class = 'cell-border',
      rownames = FALSE,
      selection = 'none',
      escape = FALSE, # to allow HTML in table
      options = list(
        columnDefs = list(list(className = 'dt-nopad sim-cell', height=38, width=120, targets = 0),
                          list(targets = elipsis_targets, render = DT::JS(
                            "function(data, type, row, meta) {",
                            "return type === 'display' && data !== null && data.length > 17 ?",
                            "'<span title=\"' + data + '\">' + data.substr(0, 17) + '...</span>' : data;",
                            "}"))),
        ordering=FALSE,
        scrollX = TRUE,
        pageLength = 20,
        paging = TRUE,
        bInfo = 0,
        dom = 'ftp'
      )
    )
  }, server = TRUE)
}

