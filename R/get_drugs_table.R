#' Get Drugs metadata table for CMAP02 or L1000
#'
#' @param study either \code{'CMAP02'} or \code{'L1000_genes'} or \code{'L1000_drugs'}
#'
#' @return \code{data.frame} of meta data for CMAP02 or L1000
#' @export
#' @keywords internal
#' @examples
#'
#' cmap_meta <- get_drugs_table('CMAP02')
#' l1000_meta <- get_drugs_table('L1000_genes')
#'
get_drugs_table <- function(study) {

  # load pdata for study
  pdata_fname <- paste0(study, '_pdata.rds')
  pdata_path <- system.file('extdata', pdata_fname,
                            package = 'drugseqr.data', mustWork = TRUE)
  pdata <- readRDS(pdata_path)

  # append BRH annotation
  drugs_table <- append_annot(pdata, study)

  # cell line used for sorting
  drugs_table$cell_line <- gsub('^[^_]+_([^_]+)_.+?$', '\\1', drugs_table$title)

  # compound used for grouping
  drugs_table <- tibble::add_column(drugs_table,
                                    Compound = stringr::str_extract(drugs_table$title, '^[^_]+'), .before = 1)

  # titles for correlation points
  drugs_table <- drugs_table %>%
    dplyr::mutate(cor_title = paste(stringr::str_replace(title, '^[^_]+_', ''), `Samples(n)`, sep = '_')) %>%
    dplyr::select(-`Samples(n)`)

  if (study == 'L1000_genes') {
    # remove nonsense for L1000 genetic and reorder
    drugs_table$cor_title <- gsub('_-700-666.0', '', drugs_table$cor_title)
    drugs_table <- drugs_table[, c('Compound', 'title', 'Genecards', 'Description', 'cell_line', 'cor_title')]
  }


  return(drugs_table)
}

#' Bind annotation data to pdata
#'
#' For the moment, annotation data is Broad Repurposing Hub data.
#'
#' @param pdata \code{data.frame} of pdata for either CMAP02 or L1000.
#' @param study Character identifying study. Either \code{'CMAP02'} or \code{'L1000'}.
#'
#' @return \code{pdata} with bound annotation columns.
#' @export
append_annot <- function(pdata, study) {
  annot_fname <- paste0(study, '_annot.rds')
  annot_path <- system.file('extdata', annot_fname,
                            package = 'drugseqr.data', mustWork = TRUE)

  annot <- readRDS(annot_path)
  pdata <- cbind(pdata, annot)
  return(pdata)
}
