server <- function(input, output, session) {

  # get arguments from calling function
  anal_name <- reactiveVal(NULL)
  sorted_query <- reactiveVal(NULL)

  query_res <- reactiveVal(getShinyOption('query_res'))
  drug_study <- reactiveVal(getShinyOption('drug_study', 'L1000 Drugs'))
  cells <- reactiveVal(getShinyOption('cells'))
  sort_by <- reactiveVal(getShinyOption('sort_by', 'avg_cor'))
  show_clinical <- reactiveVal(getShinyOption('show_clinical', TRUE))
  min_signatures <- reactiveVal(getShinyOption('min_signatures', 3))
  is_pert <- reactiveVal(getShinyOption('is_pert', FALSE))
  direction <- reactiveVal(getShinyOption('direction', 'opposing'))
  ntop <- getShinyOption('ntop', 'all')


  # the output table
  callModule(drugsTable, 'table',
             query_res = query_res,
             sorted_query = sorted_query,
             drug_study = drug_study,
             anal_name = anal_name,
             cells = cells,
             sort_by = sort_by,
             show_clinical = show_clinical,
             min_signatures = min_signatures,
             is_pert = is_pert,
             direction = direction,
             ntop = ntop)



}
