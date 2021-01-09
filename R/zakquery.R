#' Run Drug Query with Zak Specifications
#'
#' Only clinical drugs, cells types: all, haematopoietic, monocytic, and lung
#'
#' @inheritParams run_drugseqr_table
#' @param interactive if \code{TRUE} (default) then starts shiny app with first \code{cell_types}. Otherwise returns
#'  data.frames with results for all \code{cell_types}.
#'
#' @return List of data.frames or result of run_drugseqr_table
#' @keywords internal
#'
zakquery <- function(query_res,
                     drug_study = c('CMAP02', 'L1000 Drugs', 'L1000 Genes'),
                     cell_types = c('all', 'haematopoietic', 'monocytic', 'lung'),
                     min_signatures = 1,
                     interactive = TRUE) {

  # HL60: predominantly neutrophilic promyelocytic morphology
  # PL21: myeloid leukemia
  # SKM1: myeloid leukemia
  # WSUDLCL2: Diffuse large B-cell lymphoma

  # NOMO1: monocytic leukemia
  # THP1: monocytic leukemia
  # U937: monocytic leukemia
  cells <- list(
    haematopoietic = c('HL60', 'JURKAT', 'NOMO1', 'PL21', 'SKM1', 'THP1', 'U937', 'WSUDLCL2'),
    monocytic = c('NOMO1', 'THP1', 'U937'),
    lung = c('HCC515', 'A549', 'CORL23', 'DV90', 'H1299', 'HCC15', 'NCIH1694', 'NCIH1836', 'NCIH2073', 'NCIH596', 'SKLU1', 'T3M10')
  )

  if (interactive) {
    cells <- cells[[cell_types[1]]]
    run_drugseqr_table(
      query_res,
      drug_study,
      cells = cells,
      min_signatures = min_signatures)

  } else {
    res <- lapply(c('all', names(cells)),
                  function(x) format_query_res(
                    query_res,
                    drug_study,
                    cells = cells[[x]],
                    min_signatures = min_signatures))

    names(res) <- c('all', names(cells))
    return(res)
  }

}
