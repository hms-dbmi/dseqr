#' Run app with drug query table only
#'
#' Used to show full query table and use programmatic filter options
#'
#' @inheritParams format_query_res
#'
#' @export
#'
run_drugseqr_table <- function(query_res,
                               drug_study = 'L1000 Drugs',
                               cells = NULL,
                               sort_by = 'avg_cor',
                               show_clinical = TRUE,
                               min_signatures = 3,
                               is_pert = FALSE,
                               direction = 'opposing',
                               ntop = 'all',
                               app_dir = system.file('dt', package = 'drugseqr', mustWork = TRUE),
                               host = '0.0.0.0',
                               port = 3838) {


  # pass arguments to app through options then run
  shiny::shinyOptions(query_res = query_res,
                      cells = cells,
                      sort_by = sort_by,
                      show_clinical = show_clinical,
                      min_signatures = min_signatures,
                      is_pert = is_pert,
                      ntop = ntop,
                      direction = direction)



  # auto-reload if update app files
  options(shiny.autoreload = TRUE)
  shiny::runApp(app_dir, launch.browser = TRUE, host = host, port = port)

}
