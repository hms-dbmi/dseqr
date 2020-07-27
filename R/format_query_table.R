#' Format perturbation query result like web app.
#'
#' Used to get unfiltered query table from downloaded perturbation query result. A condensed top-hits is shown in the
#' web app for performance reasons.
#'
#' @param query_res Numeric vector of correlation values with names as HGNC symbols
#' @param drug_study Queried database giving rise to \code{query_res}.
#'   One of \code{'CMAP02'},\code{'L1000 Genetic'}, or \code{'L1000 Drugs'}.
#' @param sort_by Metric to sort by. Either \code{'avg_cor'} (default) or \code{'min_cor'}.
#' @param direction  Direction of correlation to sort results by. One of \code{'both'}, \code{'similar'}, or \code{'opposing'} (default).
#' @param sort_abs Should results be sorted based on absolute correlation? Default is \code{FALSE}.
#' @param show_clinical Should result only contain drugs with clinical phase annotation? Default is \code{FALSE}.
#' @param min_signatures Number of independent perturbagen signatures below which a perturbation is excluded. Default is \code{3}.
#' @param cells Character vector of cell types to include. Default (\code{NULL}) includes all cell types.
#'
#' @return \code{data.frame} with annotated, filtered, and sorted query result.
#' @export
#'
#' @examples
#' ## Not run:
#'
#' dl_res <- read.csv('path/to/downloaded_full_query_result.csv', row.names = 2)
#' query_res <- dl_res$correlation
#' names(query_res) <- row.names(dl_res)
#'
#' query_table <- format_query_res(query_res)
#'
#'
format_query_res <- function(query_res,
                               drug_study,
                               sort_by = c('avg_cor', 'min_cor'),
                               direction = c('opposing', 'similar', 'both'),
                               sort_abs = FALSE,
                               show_clinical = FALSE,
                               min_signatures = 3,
                               cells = NULL
) {

  # if didn't make selection, choose defaults
  sort_by <- sort_by[1]
  direction <- direction[1]

  is_genetic <- drug_study == 'L1000 Genetic'

  # get annotation for selected study
  annot_study <- switch(drug_study,
                        'CMAP02' = 'CMAP02',
                        'L1000 Genetic' = 'L1000_genes',
                        'L1000 Drugs' = 'L1000_drugs')

  drug_annot <- get_drugs_table(annot_study)
  query_table_annot <- annot_query_res(query_res, drug_annot)

  # subset to selected cells, summarize by compound, and add html
  query_table_summarised <- summarise_query_table(query_table_annot,
                                                  is_genetic = is_genetic,
                                                  cells = cells,
                                                  sort_abs = sort_abs,
                                                  remove_html = TRUE,
                                                  ntop = 'all')

  # subset by min signatures
  query_table_nsig <- filter_nsig(query_table_summarised, min_signatures)

  # subset by clinical phase
  query_table_clin <- filter_clinical(query_table_nsig, show_clinical)

  # final sorting/filtering
  query_table_final <- sort_query_table_clin(query_table_clin,
                                             sort_by = sort_by,
                                             sort_abs = sort_abs,
                                             direction = direction,
                                             drug_study = drug_study,
                                             remove_html = TRUE)

  return(query_table_final)

}


#' Annotate query result
#'
#' @inheritParams format_query_res
#' @param drug_annot result of \link{get_drugs_table}
#'
#' @return Annotated query table
#' @export
#' @keywords internal
annot_query_res <- function(query_res, drug_annot) {

  drug_annot <- drug_annot[drug_annot$title %in% names(query_res), ]
  query_res <- query_res[drug_annot$title]

  if (!all.equal(drug_annot$title, names(query_res))) stop('query_res not complete')

  tibble::add_column(drug_annot,
                     Rank = NA,
                     Correlation = query_res,
                     .before=0)

}

#' Summarise annotated query table
#'
#' @param query_table_annot result of \link{query_table_annot}
#' @param remove_html Should html columns be removed? For non-webapp usage.
#' @inheritParams format_query_res
#' @inheritParams get_top_cors
#'
#' @return summarised query table
#' @export
#' @keywords internal
summarise_query_table <- function(query_table_annot, is_genetic, cells, sort_abs, remove_html = FALSE, ntop = 1500) {

  if (ntop == 'all') ntop <- length(unique(query_table_annot$Compound))

  cols <- get_query_cols(is_genetic)

  res <- query_table_annot %>%
    limit_cells(cells) %>%
    summarize_compound(is_genetic = sort_abs, ntop = ntop) %>%
    add_table_html() %>%
    dplyr::select(cols, dplyr::everything())

  if (remove_html)
    res <- res %>% dplyr::select(-c('External Links', 'Correlation'))

  return(res)
}

#' Get query columns
#'
#' Gets appropriate column names for query study.
#'
#' @param is_genetic Should gene columns be returned? If \code{FALSE} drug columns are returned.
#'
#' @return Character vector of column names.
#' @export
#' @keywords internal
get_query_cols <- function(is_genetic) {
  drug_cols <- c('Rank', 'Correlation', 'Compound', 'Clinical Phase', 'External Links', 'MOA', 'Target', 'Disease Area', 'Indication', 'Vendor', 'Catalog #', 'Vendor Name')
  gene_cols <- c('Rank', 'Correlation', 'Compound', 'External Links', 'Description')

  cols <- if (is_genetic) gene_cols else drug_cols
  return(cols)
}


#' Subset by min signatures
#'
#' @param query_table_summarised  result of \link{summarise_query_table}
#' @inheritParams format_query_res
#'
#' @return \code{query_table_summarised} with drugs filtered that contain fewer than \code{min_signatures}.
#' @export
#' @keywords internal
filter_nsig <- function(query_table_summarised, min_signatures) {
  query_table_summarised %>% dplyr::filter(n >= min_signatures)
}

#' Subset by clinical phase
#'
#' @param query_table_nsig result of \link{filter_nsig}
#' @param show_clinical if \code{TRUE}, only drugs with a clinical phase annotation are kept.
#'
#' @return \code{query_table_nsig}.
#'   Drugs without a clinical phase annotation are removed if \code{show_clinical} is \code{FALSE}
#' @export
#' @keywords internal
filter_clinical <- function(query_table_nsig, show_clinical) {
  query_table_nsig %>% {
    if (show_clinical && 'Clinical Phase' %in% colnames(.))
      dplyr::filter(., !is.na(`Clinical Phase`))
    else .
  }
}

#' final sorting/filtering of query table
#'
#' @param query_table_clin result of \code{filter_clinical}
#' @inheritParams format_query_res
#' @inheritParams summarise_query_table
#'
#' @return Final query table
#' @export
#' @keywords internal
sort_query_table_clin <- function(query_table_clin, sort_by, sort_abs, direction, drug_study, remove_html = FALSE) {

  # final sorting/filtering
  q <- query_table_clin

  if (sort_by == 'avg_cor' & !remove_html)
    q$Correlation <- gsub('simplot', 'simplot show-meanline', q$Correlation)

  # show largest absolute correlations first for genetic and pert queries
  # as both directions are informative
  if (sort_abs) {

    # filter none, opposing, or similar signatures based on direction toggle
    if (sort_by == 'avg_cor') {
      is.sim <- q$avg_cor > 0

    } else {
      mm <- q[, c('min_cor', 'max_cor')]
      mcol <- max.col(abs(mm), ties.method = 'last')
      is.sim <- mcol == 2 & mm[, 2] > 0
    }

    q <- switch(direction, 'both' = q, 'similar' = q[is.sim, ], 'opposing' = q[!is.sim, ]) %>%
      dplyr::mutate(min_cor = -pmax(abs(min_cor), abs(max_cor))) %>%
      dplyr::mutate(avg_cor = -abs(avg_cor))
  }

  # indicate total number of unique perts in title for rank
  # skip if not for web app

  rank_title <- switch(drug_study,
                       'CMAP02' = 'out of 1,309',
                       'L1000 Genetic' = 'out of 6,943',
                       'L1000 Drugs' = 'out of 19,360')

  # sort as desired then add rank
  q <- q %>%
    dplyr::arrange(!!sym(sort_by)) %>%
    dplyr::select(-min_cor, -avg_cor, -max_cor, -n)

  if (!remove_html) q <- q %>%
    dplyr::mutate(Rank = paste0('<span class="rank-label label label-default" title="', rank_title, '">', 1:nrow(q), '</span>'))

  return(q)
}


