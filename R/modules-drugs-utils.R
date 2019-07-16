
#' Limit Query Results to Specific Cell Lines
#'
#' @param query_table \code{data.frame} with \code{cell_line} column.
#' @param cells Character vector of cell lines to limit \code{query_table} by.
#' @importFrom magrittr "%>%"
#'
#' @return \code{query_table} for specified \code{cells}
#' @export
limit_cells <- function(query_table, cells) {
  if (!is.null(cells))
    query_table <- query_table %>%
      dplyr::filter(cell_line %in% cells)

  query_table %>%
    dplyr::select(-cell_line)
}


#' Reduce entries in query results
#'
#' Gets all entries for compounds that have a correlation less than the top entries with a clinical phase.
#'
#' @param query_cors \code{data.frame} with columns \code{'Compound'} and \code{arrange_by}.
#' @param arrange_by String indicating column name in \code{query_cors} do sort by.
#' @param ntop Integer indicating the number of rows to keep after sorting by \code{arrange_by}.
#'
#' @return Character vector of \code{ntop} compounds.
#' @export
#'
#' @importFrom magrittr "%>%"
get_top <- function(query_cors, arrange_by, ntop = 300) {
  query_cors %>%
    dplyr::as_tibble() %>%
    dplyr::arrange(!!sym(arrange_by)) %>%
    head(ntop) %>%
    dplyr::pull(Compound)
}


#' Summarize query results and annotations by perturbation
#'
#' Takes a \code{data.frame} with one row per signatures and summarizes to one row per compound.
#'
#' Variables related to individual signatures (cell line, dose, duration, and sample number) are
#' pasted together and added as a list to \code{'title'} column. Query correlation values are also added as a list to
#' the \code{'Correlation'} column.
#'
#' Clinical status is summarized by keeping the most advanced phase only (e.g. Launched > Phase 3). For all other variables,
#' all unique entries are paste together seperated by \code{'|'}.
#'
#' @param query_table \code{data.frame} of perturbation correlations and annotations.
#'
#' @return \code{data.frame} of perturbation correlations and annotations summarized by perturbation.
#' @export
#' @keywords internal
#'
#' @importFrom magrittr "%>%"
#'
summarize_compound <- function(query_table) {

  # group by compound
  query_table <- query_table %>%
    dplyr::group_by(Compound)

  # put all correlations together in list
  # keep minimum and average correlation for sorting
  query_cors <- query_table %>%
    dplyr::select(Correlation, Compound) %>%
    dplyr::summarise(min_cor = min(Correlation),
                     avg_cor = mean(Correlation),
                     Correlation = I(list(Correlation)))

  # compounds in top min or avg cor
  top_min <- get_top(query_cors, 'min_cor')
  top_avg <- get_top(query_cors, 'avg_cor')
  top <- unique(c(top_min, top_avg))

  # filter based on top
  query_table <- query_table %>%
    dplyr::filter(Compound %in% top)

  query_cors <- query_cors %>%
    dplyr::filter(Compound %in% top) %>%
    dplyr::select(-Compound)

  query_title <- query_table %>%
    dplyr::select(cor_title, Compound) %>%
    dplyr::summarize(cor_title = I(list(cor_title))) %>%
    dplyr::select(cor_title)

  # keep furthest clinical phase
  summarise_phase <- function(phases) {
    if (all(is.na(phases))) return(phases[1])
    max(phases, na.rm = TRUE)
  }

  query_phase <- query_table %>%
    dplyr::select(`Clinical Phase`, Compound) %>%
    dplyr::summarise(`Clinical Phase` = summarise_phase(`Clinical Phase`)) %>%
    dplyr::pull(`Clinical Phase`)


  # summarize rest
  query_rest <- query_table %>%
    dplyr::select(-Correlation, -`Clinical Phase`, -title, -cor_title) %>%
    dplyr::summarize_all(function(x) {
      unqx <- na.omit(x)
      unqx <- as.character(unqx)
      unqx <- unlist(strsplit(unqx, '|', fixed = TRUE))
      unqx <- unique(unqx)

      # keep as NA if they all are
      if (!length(unqx)) return(NA_character_)

      # collapse distinct non-NA entries
      return(paste(unqx, collapse = ' | '))
    }) %>%
    tibble::add_column('Clinical Phase' = query_phase, .after = 'Compound')

  query_table <- dplyr::bind_cols(query_cors, query_rest, query_title)
  return(query_table)
}

#' Add linkout HTML
#'
#' Non \code{NA} values in \code{id_col} of \code{query_res} are inserted between
#' \code{pre_url} and \code{post_url} to form hyperlinks. Relevant HTML markup is also added.
#'
#'
#' @param query_res \code{data.frame} returned by \code{\link{summarize_compound}}.
#' @param id_col Character with column in \code{query_res} that contains ids to
#'   be inserted between \code{pre_url} and \code{post_url} to form the link. \code{NA}
#'   values will be ignored.
#' @param img_url Character with url to an image to display as hyperlink.
#' @param pre_url Character with url portion to paste before \code{id_col} column values.
#' @param post_url Character with url portion to paste before \code{id_col} column values.
#' @param title Character that will be added to hyperlink title attribute. Default is \code{id_col}.
#'
#' @return \code{query_res} with HTML for hyperlinks in \code{id_col}.
#' @export
add_linkout <- function(query_res, id_col, img_url, pre_url, post_url = NULL, title = id_col) {

  ids <- query_res[[id_col]]
  have_ids <- !is.na(ids)

  query_res[[id_col]][have_ids] <- paste0(get_open_a(pre_url, ids[have_ids], post_url, title),
                                          '<img src="', img_url, '" height="20px" hspace="4px"></img>',
                                          # ids[have_ids],
                                          '</a>')

  return(query_res)
}

#' Get Opening HTML a Tag
#'
#' If there are multiple-entry \code{ids} (e.g. \code{'2250, | 60823'}), the first entry
#' is added to the href attribute and subsequent entries are added to onclick javascript.
#'
#' @inheritParams add_linkout
#' @param ids non \code{NA} ids to be inserted between \code{pre_url} and \code{post_url} to form the link.
#'
#' @return Character vector of opening HTML a tags
#' @export
get_open_a <- function(pre_url, ids, post_url, title) {

  # are some cases with e.g. multiple pubchem cids
  ids <- strsplit(ids, ' | ', fixed = TRUE)

  open_a <- sapply(ids, function(id) {
    # add first id to href
    res <- paste0('<a href="', pre_url, id[1], post_url, '" target="_blank"')

    if (length(id) > 1) {
      # open next ids with javascript
      onclick_js <- paste0("window.open('", pre_url, id[-1], post_url, "')", collapse = '; ')
      res <- paste0(res, ' onclick="', onclick_js, '" title="', paste('Go to', title, paste(id, collapse = ', '), '">'))
    } else {
      res <- paste0(res, ' title="', paste('Go to', title, '">'))
    }
    return(res)
  })
}

#' Add HTML to query results table
#'
#' @param query_res \code{data.frame} returned by \code{\link{summarize_compound}}
#'
#' @return \code{query_res} with pubchem cid links and correlation plot HTML.
#' @export
add_table_html <- function(query_res) {

  pre_urls <- c('https://pubchem.ncbi.nlm.nih.gov/compound/',
                'http://sideeffects.embl.de/drugs/',
                'https://www.drugbank.ca/drugs/',
                'https://en.wikipedia.org/wiki/')

  img_urls <- c('https://pubchem.ncbi.nlm.nih.gov/pcfe/favicon/favicon.ico',
                'http://sideeffects.embl.de/media/images/EMBL_Logo.png',
                'https://www.drugbank.ca/favicons/favicon.ico',
                'https://en.wikipedia.org/static/favicon/wikipedia.ico')


  # add linkout to Pubchem, SIDER, and DrugBank
  query_res <- add_linkout(query_res, 'Pubchem CID', img_urls[1], pre_urls[1], title = 'Pubchem')
  query_res <- add_linkout(query_res, 'SIDER', img_urls[2], pre_urls[2])
  query_res <- add_linkout(query_res, 'DrugBank', img_urls[3], pre_urls[3])
  query_res <- add_linkout(query_res, 'Wikipedia', img_urls[4], pre_urls[4])

  # merge linkouts into single column
  query_res <- merge_linkouts(query_res, c('Pubchem CID', 'Wikipedia', 'DrugBank', 'SIDER'))

  # replace correlation with svg element
  cors <- query_res$Correlation

  # move titles to plots
  titles <- query_res$cor_title
  query_res$cor_title <- NULL

  cors_range <- range(unlist(cors))
  # centerline x
  cx <- calcx(0, cors_range)

  query_res$Correlation <- paste0('<svg class="simplot" width="180" height="38">
                            <line class="centerline" x1="',cx,'" x2="',cx,'" y1="0" y2="38"></line>',
                                  get_cors_html(cors, titles, cors_range),
                                  '</svg>')

  return(query_res)
}

#' Merge columns with image links
#'
#' @param query_res \code{data.frame} after calling \code{\link{add_linkout}} to \code{cols}.
#' @param cols Character vector of columns in \code{query_res} that \code{\link{add_linkout}} has been called on.
#
#' @importFrom magrittr "%>%"
#'
#' @return \code{query_res} with column \code{'External Links'} formed from pasting \code{cols} together. \code{cols} are removed.
#' @export
merge_linkouts <- function(query_res, cols) {

  # paste cols with non-NA values
  paste.na <- function(x) paste(x[!is.na(x)], collapse = '')
  new_vals <- apply(query_res[ ,cols] , 1, paste.na)

  # remove cols that pasted
  query_res <- query_res %>%
    tibble::add_column(`External Links` = new_vals, .before = cols[1]) %>%
    dplyr::select(-cols)

  return(query_res)
}

#' Get HTML for correlation values.
#'
#' @param cors List of numeric vectors of pearson correlations.
#' @param titles List of character vectors of treatment titles for pearson correlations (e.g. MCF7_1e-05M_6h_3).
#' @param cors_range Numeric vector of length two specifying the range of correlation values.
#'
#' @return Character vector of HTML markup for the title/circle/text for a correlation plot.
#' @export
get_cors_html <- function(cors, titles, cors_range) {

  cors_html <- sapply(seq_along(cors), function(i) {
    x <- cors[[i]]
    xtitle <- titles[[i]]
    paste0('<g class="cor-point">
              <title>', xtitle, '</title>
              <g><text x="', calcx(x, cors_range), '" y="38" class="x text" dy="-2">', signif(x, 3), '</text></g>
              <g><circle cx="', calcx(x, cors_range), '" cy="19" r="5" class="cor"></circle></g>
            </g>', collapse = '\n')
  })

  return(cors_html)
}




#' Calculate x position for correlation plot
#'
#' @param cor Numeric vector of correlation values.
#' @param range Numeric vector of length 2 specifying maximum and minimum values of \code{cor}.
#' @param width Plot width to scale correlation values to.
#' @param pad Numeric value that is respectively, subtracted and added to values in \code{range}. Make it so that circles and
#' correlation text values don't get cut off.
#'
#' @return Numeric vector giving x position for correlation plot in \code{\link{explore_search}}
#' @export
calcx <- function(cor, range = c(-1, 1), width = 180, pad = 0.1) {
  range[1] <- range[1] - pad
  range[2] <- range[2] + pad
  (cor - range[1])/diff(range) * width
}
