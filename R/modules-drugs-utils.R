#' Get file paths to drug query results
#'
#' @param data_dir Folder with drug query results
#' @param suffix string that appears after e.g. \code{'cmap_res_'} in the drug query result filenames.
#'
#' @return List with file paths to drug query results
#'
#' @keywords internal
get_drug_paths <- function(data_dir, suffix) {
  suffix <- paste0(suffix, '.rds')
  list(
    cmap = file.path(data_dir, paste0('cmap_res_', suffix)),
    l1000_drugs = file.path(data_dir, paste0('l1000_drugs_res_', suffix)),
    l1000_genes = file.path(data_dir, paste0('l1000_genes_res_', suffix))
  )
}


#' Run custom query
#'
#' @param query_genes Named list of character vectors with \code{'up'} indicating genes to upregulate and \code{'dn'}
#'  containing genes to downregulate.
#' @param res_paths List of paths to save results to.
#' @param session session object from Shiny. Used to indicate progress.
#'
#' @return List with query results.
#'
#' @keywords internal
run_custom_query <- function(query_genes, res_paths, session) {
  res <- list()

  progress <- Progress$new(session, min = 0, max = 4)
  progress$set(message = "Querying drugs", value = 1)
  on.exit(progress$close())

  cmap_es  <- drugseqr.data::load_drug_es('cmap_es_ind.rds')
  progress$inc(1)

  # get correlations between query and drug signatures
  res$cmap <- query_budger(query_genes, cmap_es)
  rm(cmap_es)

  saveRDS(res$cmap, res_paths$cmap)
  saveRDS(query_genes, res_paths$query_genes)

  # don't run l1000 if no genes from it selected
  run.l1000 <- any(unlist(query_genes) %in% genes$common)
  progress$inc(ifelse(run.l1000, 1, 2))

  if (run.l1000) {

    l1000_drugs_es <- drugseqr.data::load_drug_es('l1000_drugs_es.rds')
    res$l1000_drugs <- query_budger(query_genes, l1000_drugs_es)
    rm(l1000_drugs_es)

    l1000_genes_es <- drugseqr.data::load_drug_es('l1000_genes_es.rds')
    res$l1000_genes <- query_budger(query_genes, l1000_genes_es)
    rm(l1000_genes_es)

    saveRDS(res$l1000_drugs, res_paths$l1000_drugs)
    saveRDS(res$l1000_genes, res_paths$l1000_genes)
    progress$inc(1)
  }

  return(res)
}

#' Validate custom query
#'
#' @param dn_genes Character vector of genes to upregulate
#' @param up_genes Character vector of genes to downregulate
#' @param custom_name Name for new custom analysis
#'
#' @return \code{NULL} if valid, otherwise an error message
#'
#' @keywords internal
validate_custom_query <- function(dn_genes, up_genes, custom_name) {
  msg <- NULL

  # genes in both lists
  in.both <- intersect(dn_genes, up_genes)

  # genes not in pert studies
  na.genes <- setdiff(c(dn_genes, up_genes), c(genes$common, genes$cmap_only))

  # make sure custom name provided
  if (is.null(custom_name) || custom_name == '') {
    msg <- 'Provide a name for custom query'

  } else if (length(in.both)) {
    msg <- paste0('Genes must be in only one of up/down lists: ',
                  paste(in.both, collapse = ', '))

  } else if (length(na.genes)) {
    msg <- paste0('Genes not measured in perturbation studies: ',
                  paste(na.genes, collapse = ', '))

  } else if (!length(c(dn_genes, up_genes))) {
    msg <- 'Please provide genes to up and/or down regulate'
  }

  return(msg)
}

#' Validate uploaded custom query
#'
#' @param top_table data.frame with logFC, dprime, or effect sizes as column 1
#'   and row.names as HGNC symbols
#' @param custom_name Name for new custom analysis
#'
#' @return \code{NULL} if valid, otherwise an error message
#'
#' @keywords internal
validate_up_custom <- function(top_table, custom_name) {
  msg <- NULL

  if (is.null(top_table)) return('Need two columns: unique HGNC and effect size')

  # need at least one L1000/CMAP gene
  genes <- unique(unlist(genes))
  common <- intersect(genes, row.names(top_table))

  if (!length(common)) {
    msg <- 'First column: HGNC symbols in L1000/CMAP02'
  }

  # need logFC/dprimes effect size as second column
  tt_cols <- colnames(top_table)
  dpcol <- match.arg(c('dprime', 'logFC', tt_cols[1]), tt_cols, several.ok = TRUE)[1]
  es <- top_table[[dpcol]]

  es.abs <- max(abs(es))
  if (is.character(es)) {
    msg <- 'Second column: numeric logFC or dprime values'

  } else if (es.abs > 1000) {
    msg <- 'Second column: values too large'

  } else if (es.abs <= 1) {
    msg <- 'Second column: values too small'
  }

  return(msg)

}

#' Format uploaded custom query signature
#'
#' @inheritParams validate_up_custom
#'
#' @return \code{top_table} subsetted to rows with HGNC symbols in CMAP02/L1000 and
#'  with dprime column set to either \code{top_table$dprime}, \code{top_table$logFC},
#'  or \code{top_table[,1]}.
#'
#' @keywords internal
#'
format_up_custom <- function(top_table) {

  # use common genes
  genes <- unique(unlist(genes))
  common <- intersect(genes, row.names(top_table))
  top_table <- top_table[row.names(top_table) %in% common,, drop = FALSE]

  # set logfc/dprime to dprime for drug queries
  tt_cols <- colnames(top_table)
  dpcol <- match.arg(c('dprime', 'logFC', tt_cols[1]), tt_cols, several.ok = TRUE)[1]
  top_table$dprime <- top_table[[dpcol]]
  return(top_table)
}



#' Load results of previous custom query
#'
#' Will only load \code{res_paths} that exist. L1000 results won't exist
#' if no genes selected for custom query were in the L1000 measured genes.
#'
#' @param res_paths List of named strings with paths to results.
#'
#' @return List of names results loaded from \code{res_paths}.
#'
#' @keywords internal
load_custom_results <- function(res_paths, is_pert) {
  res <- list()
  for (i in seq_along(res_paths)) {
    res_path <- res_paths[[i]]
    res_name <- names(res_paths)[i]

    # download requested pert result
    if (!file.exists(res_path) && is_pert && !grepl('diff_expr_symbol_.+?.rds$', res_path))
      dl_pert_result(res_path)

    if (file.exists(res_path))
      res[[res_name]] <- readRDS(res_path)
  }
  return(res)
}

#' Download CMAP02/L1000 pert query result from S3
#'
#' @param res_path Path to download file to.
#'
#' @return NULL
#' @keywords internal
#' @examples
#' data_dir <- tempdir()
#' res_path <- file.path(data_dir, 'cmap_res_BRD-K45319408_PC3_5um_24h.rds')
#' drugseqr:::dl_pert_result(res_path)
#'
dl_pert_result <- function(res_path) {
  # name of the file being requested
  if (file.exists(res_path)) return(NULL)
  dl_url <- paste0('https://s3.us-east-2.amazonaws.com/drugseqr/pert_query_dir/', basename(res_path))
  dl_url <- utils::URLencode(dl_url)
  dl_url <- gsub('+', '%2B', dl_url, fixed = TRUE)
  utils::download.file(dl_url, res_path)
}

#' Load data.frame of custom anals for selectizeInput
#'
#' @param data_dir Path to directory containing \code{'custom_queries'} folder.
#'
#' @return \code{data.frame} used to show available custom queries in Drugs tab
#'
#' @keywords internal
load_custom_anals <- function(data_dir) {
  custom_dir <- file.path(data_dir, 'custom_queries')

  anals <- NULL
  if (dir.exists(custom_dir)) {
    custom_names <- list.files(custom_dir, pattern = '^cmap_res_.+?.rds')
    custom_names <- gsub('^cmap_res_(.+?).rds$', '\\1', custom_names)


    anals <- data.frame(matrix(ncol = 5, nrow = 0), stringsAsFactors = FALSE)
    colnames(anals) <- c("dataset_name", "dataset_dir", "label", "value", "type")
    custom_dir <- file.path(data_dir, 'custom_queries')

    for (i in seq_along(custom_names))
      anals[i, ] <- c(NA, 'custom_queries', custom_names[i], custom_names[i], 'Custom')

  }

  return(anals)
}

#' Load data.frame of CMAP02/L1000 perturbations for right click load signature on correlation points
#'
#' @return \code{data.frame} used to show available pertubation queries in Drugs tab
#'
#' @keywords internal
load_pert_anals <- function() {
  anals <- data.frame(matrix(NA, ncol = 5, nrow = length(pert_names)), stringsAsFactors = FALSE)
  colnames(anals) <- c("dataset_name", "dataset_dir", "label", "value", "type")

  anals$label <- anals$value <- pert_names
  anals$type <- 'CMAP02/L1000 Perturbations'

  return(anals)
}


#' Get cell choices data.frame for CMAP02 or L1000
#'
#' @param drug_study either 'CMAP02', 'L1000 Drugs', or 'L1000 Genetic'
#'
#' @return data.frame
#'
#' @keywords internal
get_cell_choices <- function(drug_study) {
  if (drug_study == 'CMAP02') return(cell_info$cmap)
  else if (drug_study == 'L1000 Drugs') return(cell_info$l1000_drugs)
  else if (drug_study == 'L1000 Genetic') return(cell_info$l1000_genes)
}


#' Limit Query Results to Specific Cell Lines
#'
#' @param query_table \code{data.frame} with \code{cell_line} column.
#' @param cells Character vector of cell lines to limit \code{query_table} by. If \code{NULL} (default),
#'  no filtering occurs.
#' @importFrom magrittr "%>%"
#'
#' @return \code{query_table} for specified \code{cells}
#' @keywords internal
limit_cells <- function(query_table, cells = NULL) {
  cell_line <- NULL

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
#' @param decreasing Should compounds be ordered from most similar to most
#'   dissimilar. \code{TRUE} used for genetic perturbations.
#'
#' @return Character vector of \code{ntop} compounds.
#' @keywords internal
#'
#' @importFrom magrittr "%>%"
get_top <- function(query_cors, arrange_by, ntop, decreasing = FALSE) {
  Compound <- NULL

  pre <- ifelse(decreasing, -1, 1)
  query_cors %>%
    dplyr::as_tibble() %>%
    dplyr::arrange(pre*!!rlang::sym(arrange_by)) %>%
    utils::head(ntop) %>%
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
#' @param is_genetic is \code{query_table} from L1000 genetic perts?
#' @inheritParams get_top
#'
#' @return \code{data.frame} of perturbation correlations and annotations summarized by perturbation.
#'
#' @keywords internal
#'
#' @importFrom magrittr "%>%"
#' @importFrom rlang :=
#'
summarize_compound <- function(query_table, is_genetic = FALSE, ntop = 1500) {
  Compound <- `Clinical Phase` <- .SD <- . <- NULL

  query_table <- data.table::data.table(query_table, key = 'Compound')


  query_cors <- get_top_cors(query_table, is_genetic = is_genetic, ntop = ntop)
  query_table <- query_table[Compound %in% query_cors$Compound, ]

  # split joined cols
  joined_cols <- intersect(colnames(query_table), c('MOA', 'Target', 'Disease Area', 'Indication', 'Vendor', 'Catalog #', 'Vendor Name'))
  if (length(joined_cols))
    query_table[,
                (joined_cols) := lapply(.SD, strsplit, ' | ', fixed = TRUE),
                .SDcols = joined_cols]

  # summarize rest cols
  rest_cols <- setdiff(colnames(query_table), c('Compound', 'Correlation', 'Clinical Phase', 'cor_title', 'title'))
  query <- c('title = I(list(title))',
             'cor_title = I(list(cor_title))',
             paste0('`', rest_cols, '`', ' = paste(unique(unlist(`', rest_cols, '`)), collapse = " | ")'))
  query <- paste(query, collapse = ', ')
  query <- paste0('query_table[, .(', query, '), by = Compound]')

  query_rest <- eval(parse(text = query))
  query_rest[query_rest == 'NA'] <- NA_character_


  # not in L1000 Genetic
  if ('Clinical Phase' %in% colnames(query_table)) {

    summarise_phase <- function(phases) {
      if (all(is.na(phases))) return(phases[1])
      max(phases, na.rm = TRUE)
    }
    query_phase <- query_table[,
                               .(summarise_phase(`Clinical Phase`)),
                               by = Compound]$V1


    query_rest <- query_rest %>%
      tibble::add_column('Clinical Phase' = query_phase, .after = 'Compound')

  }

  query_rest <- query_rest[, -('Compound')]


  query_table <- dplyr::bind_cols(query_cors, query_rest) %>%
    as.data.frame()

  return(query_table)
}

#' Get top correlated compounds for drug query
#'
#' @param query_table data.frame with columns \code{'Correlation'} and \code{'Compound'}.
#' @param is_genetic Boolean indicating if the query results are from L1000 genetic perturbations.
#' @inheritParams get_top
#'
#' @return data.table summarized by Compound with top results.
#'
#' @keywords internal
get_top_cors <- function(query_table, ntop, is_genetic = FALSE) {
  Correlation <- .N <- Compound <- NULL

  query_table <- data.table::data.table(query_table, key = 'Compound')

  # put all correlations together in list
  query_cors <- query_table[, list(Correlation = I(list(Correlation)),
                                   max_cor = max(Correlation),
                                   min_cor = min(Correlation),
                                   avg_cor = mean(Correlation),
                                   n = .N),
                            by = Compound]


  # compounds in top min or avg cor
  top_max <- if (is_genetic) get_top(query_cors, 'max_cor', decreasing=TRUE, ntop=ntop) else NULL
  top_min <- get_top(query_cors, 'min_cor', ntop=ntop)

  top_avg_sim <- if (is_genetic) get_top(query_cors, 'avg_cor', decreasing=TRUE, ntop=ntop) else NULL
  top_avg_dis <- get_top(query_cors, 'avg_cor', ntop=ntop)

  top <- unique(c(top_min, top_max, top_avg_sim, top_avg_dis))

  # filter based on top
  query_cors <- query_cors[Compound %in% top, ]
  return(query_cors)
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
#' @keywords internal
add_linkout <- function(query_res, id_col, img_url, pre_url, post_url = NULL, title = id_col) {

  if (!id_col %in% colnames(query_res)) return(query_res)

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
#' @keywords internal
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
#' @keywords internal
add_table_html <- function(query_res) {

  pre_urls <- c('https://pubchem.ncbi.nlm.nih.gov/compound/',
                'http://sideeffects.embl.de/drugs/',
                'https://www.drugbank.ca/drugs/',
                'https://en.wikipedia.org/wiki/',
                'https://www.genecards.org/cgi-bin/carddisp.pl?gene=')

  img_urls <- c('pubchem_logo.ico',
                'EMBL_Logo.png',
                'drugbank_logo.ico',
                'wiki_logo.svg',
                'genecards_logo.ico')


  # add linkout to Pubchem, SIDER, and DrugBank
  query_res <- add_linkout(query_res, 'Pubchem CID', img_urls[1], pre_urls[1], title = 'Pubchem')
  query_res <- add_linkout(query_res, 'SIDER', img_urls[2], pre_urls[2])
  query_res <- add_linkout(query_res, 'DrugBank', img_urls[3], pre_urls[3])
  query_res <- add_linkout(query_res, 'Wikipedia', img_urls[4], pre_urls[4])
  query_res <- add_linkout(query_res, 'Genecards', img_urls[5], pre_urls[5])

  # merge linkouts into single column
  cols_to_merge <- intersect(colnames(query_res), c('Pubchem CID', 'Wikipedia', 'DrugBank', 'SIDER', 'Genecards'))
  query_res <- merge_linkouts(query_res, cols_to_merge)

  # replace correlation with svg element
  cors <- query_res$Correlation

  # move cor titles to plots
  titles <- query_res$title
  cor_titles <- query_res$cor_title
  query_res$title <- query_res$cor_title <- NULL


  cors_range <- range(unlist(cors))
  xcenter <- calcx(0, cors_range)
  xmean <- calcx(query_res$avg_cor, cors_range)
  compounds <- query_res$Compound

  query_res$Correlation <- paste0('<svg class="simplot" width="180" height="38">',
                                  paste0('<line class="meanline" x1="', xmean,'" x2="', xmean,'" y1="0" y2="8"></line>'),
                                  '<line class="centerline" x1="', xcenter,'" x2="', xcenter,'" y1="0" y2="38"></line>',
                                  get_cors_html(cors, titles, cor_titles, cors_range),
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
#' @keywords internal
merge_linkouts <- function(query_res, cols) {

  # paste cols with non-NA values
  paste.na <- function(x) paste(x[!is.na(x)], collapse = '')
  new_vals <- apply(query_res[ ,cols, drop = FALSE] , 1, paste.na)

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
#' @param cor_titles List of character vectors used for title element. Used
#'   to allow right click on correlation point to load query for perturbation signature.
#' @param cors_range Numeric vector of length two specifying the range of correlation values.
#'
#' @return Character vector of HTML markup for the title/circle/text for a correlation plot.
#' @keywords internal
get_cors_html <- function(cors, titles, cor_titles, cors_range) {


  cors_html <- sapply(seq_along(cors), function(i) {
    x <- cors[[i]]
    title <- titles[[i]]
    cor_title <- cor_titles[[i]]
    paste0('<g class="cor-point">
           <title>', cor_title, '</title>
           <g><text x="', calcx(x, cors_range), '" y="38" class="x text" dy="-2">', signif(x, 3), '</text></g>
           <g><circle cx="', calcx(x, cors_range), '" cy="19" r="5" class="cor" title="', title, '"></circle></g>
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
#' @return Numeric vector giving x position for correlation plot in Drugs tab.
#' @keywords internal
calcx <- function(cor, range = c(-1, 1), width = 180, pad = 0.1) {
  range[1] <- range[1] - pad
  range[2] <- range[2] + pad
  (cor - range[1])/diff(range) * width
}


#' Generate plotly of dprimes values for Drugs tab.
#'
#' Also used for pseudobulk single-cell datasets.
#'
#' @param path_df result of \link{get_path_df}.
#'
#' @return plotly
#'
#' @keywords internal
plot_dprimes <- function(path_df, drugs = TRUE) {


  # hovertemplate tooltips miss-behaves when single row
  # fix is direct substitution
  if (nrow(path_df) == 1) {
    sd <- path_df$sd
    fdr <- path_df$fdr
    pval <- path_df$pval
    text <- path_df$Gene
    dprime <- path_df$Dprime
    logfc <- path_df$logfc
    ambient <- path_df$ambient
    direction <- path_df$direction
    description <- path_df$description

  } else {
    sd <- '%{customdata.sd:.2f}'
    fdr <- '%{customdata.fdr}'
    pval <- '%{customdata.pval}'
    text <- '%{text}'
    dprime <- '%{x:.2f}'
    logfc <- '%{customdata.logfc:.2f}'
    ambient <- '%{customdata.ambient}'
    direction <- '%{customdata.direction}'
    description <- '%{customdata.description}'
  }

  if (drugs) {
    fontsize = 12
    pergene = 25
    title <- 'Standardized Effect Size for Top Query Genes'
    hovertemplate = paste0(
      '<span style="color: crimson; font-weight: bold; text-align: left;">Gene</span>: ', text, '<br>',
      '<span style="color: crimson; font-weight: bold; text-align: left;">Description</span>: ', description, '<br>',
      '<span style="color: crimson; font-weight: bold; text-align: left;">Dprime</span>: ', dprime, '<br>',
      '<span style="color: crimson; font-weight: bold; text-align: left;">SD</span>: ', sd, '<br>',
      '<span style="color: crimson; font-weight: bold; text-align: left;">FDR</span>: ', fdr, '<br>',
      '<span style="color: crimson; font-weight: bold; text-align: left;">Pvalue</span>: ', pval,
      '<extra></extra>')

  } else {
    fontsize = 14
    pergene = 35
    title <- 'Standardized Pseudobulk Effect Size for Each Cluster'
    hovertemplate = paste0(
      '<span style="color: crimson; font-weight: bold; text-align: left;">Cluster</span>: ', text, '<br>',
      '<span style="color: crimson; font-weight: bold; text-align: left;">Dprime</span>: ', dprime, '<br>',
      '<span style="color: crimson; font-weight: bold; text-align: left;">SD</span>: ', sd, '<br>',
      '<span style="color: crimson; font-weight: bold; text-align: left;">logFC</span>: ', logfc, '<br>',
      '<span style="color: crimson; font-weight: bold; text-align: left;">FDR</span>: ', fdr, '<br>',
      '<span style="color: crimson; font-weight: bold; text-align: left;">Pvalue</span>: ', pval, '<br>',
      '<span style="color: crimson; font-weight: bold; text-align: left;">Ambient</span>: ', ambient,
      '<extra></extra>')
  }

  xrange <- c(min(floor(c(path_df$Dprime - path_df$sd, path_df$dprime_sum)), -1),
              max(ceiling(c(path_df$Dprime + path_df$sd, path_df$dprime_sum)),  1))


  # 30 pixels width per gene in pathway
  genes <- unique(path_df$Gene)
  genes <- as.character(genes)
  ngenes <- length(genes)
  plot_height <- max(400, ngenes*pergene + 125)
  customdata <- apply(path_df, 1, as.list)

  left <- max(nchar(genes)+5)*7

  (pl <- plotly::plot_ly(data = path_df,
                         y = ~Gene,
                         x = ~Dprime,
                         text = ~Gene,
                         customdata = customdata,
                         type = 'scatter',
                         mode = 'markers',
                         height = plot_height,
                         marker = list(size = 6, color = path_df$color, opacity = path_df$opacity),
                         error_x = ~list(array = sd, color = path_df$color, thickness = 0.5, width = 0, opacity = path_df$opacity),
                         hoverlabel = list(bgcolor = '#000000', align = 'left'),
                         hovertemplate = hovertemplate
  ) %>%
      plotly::config(displayModeBar = 'hover',
                     displaylogo = FALSE,
                     modeBarButtonsToRemove = c('lasso2d',
                                                'select2d',
                                                'toggleSpikelines',
                                                'hoverClosestCartesian',
                                                'hoverCompareCartesian'),
                     toImageButtonOptions = list(format = "png", filename = 'blah')) %>%
      plotly::layout(hoverdistance = -1,
                     hovermode = 'y',
                     margin = list(t = 80, r = 80, l = left, pad = 0, autoexpand = FALSE),
                     title = list(text = title, y = .95, x = 0),
                     xaxis = list(fixedrange = TRUE,
                                  range = xrange,
                                  rangemode = "tozero",
                                  side = 'top',
                                  title = '',
                                  tickfont = list(size = fontsize)),
                     yaxis = list(fixedrange = TRUE,
                                  title = '',
                                  range = c(ngenes, -1),
                                  tickmode = 'array',
                                  tickvals = 0:ngenes,
                                  ticktext = ~Link,
                                  tickfont = list(size = fontsize)),
                     autosize = TRUE))


  # add arrow to show drug effect
  if ('dprime_sum' %in% colnames(path_df))
    pl <- pl %>%
    plotly::add_annotations(x = ~dprime_sum,
                            y = ~Gene,
                            xref = "x", yref = "y",
                            axref = "x", ayref = "y",
                            text = "",
                            showarrow = TRUE,
                            arrowcolor = ~arrow_color,
                            arrowwidth = 1,
                            ay = ~Gene,
                            ax = ~Dprime)

  return(pl)
}

