#' Validate pdata before differential expression analysis
#'
#' @param pdata data.frame with column \code{'Group'}
#'
#' @return NULL if valid, otherwise a character vector indicating what's wrong
#' @export
#' @keywords internal
validate_pdata <- function(pdata) {
  group <- pdata$Group
  group <- group[!is.na(group)]

  if (length(unique(group)) != 2) {
    msg <- 'Analysis requires test and control groups'

  } else if (length(group) < 3) {
    msg <- 'At least three samples are required for analysis'

  } else {
    msg <- NULL
  }
  return(msg)
}

#' Generate boxplotly for vsd normalized gene and cell-type deconvolution plots.
#'
#' @param df \code{data.frame} with columns: \itemize{
#'  \item x Factor for x-labels.
#'  \item y Numeric column used for y-values.
#'  \item text Character column used for hoverinfo.
#'  \item name Character column used for legend names of \code{x} values.
#'  \item color Factor column used to generate ordered colors for boxplots for each \code{'x'} value.
#' }
#' @param boxgap Used for plotly layout.
#' @param boxgroupgap Used for plotly layout.
#' @param plot_fname Name to save plot as.
#' @param ytitle Y axis title.
#' @param xtitle X axis title.
#'
#' @return plotly
#' @export
#' @keywords internal
boxPlotly <- function(df, boxgap, boxgroupgap, plot_fname, ytitle, xtitle) {

  # legend styling
  l <- list(
    font = list(
      family = "sans-serif",
      size = 12,
      color = "#000"),
    bgcolor = "#f8f8f8",
    bordercolor = "#e7e7e7",
    borderwidth = 1)


  df %>%
    plotly::plot_ly() %>%
    plotly::add_trace(x = ~x,
                      y = ~y,
                      color = ~color,
                      text = ~text,
                      name = ~name,
                      type = 'violin',
                      points = 'all',
                      jitter = 1,
                      pointpos = 0,
                      fillcolor = 'transparent',
                      hoverinfo = 'text',
                      hoveron = 'points',
                      marker = list(color = "rgba(0, 0, 0, 0.6)")) %>%
    plotly::layout(violinmode = 'group', violingroupgap = boxgroupgap, violingap = boxgap,
                   xaxis = list(fixedrange=TRUE, title = xtitle),
                   yaxis = list(fixedrange=TRUE, title = ytitle),
                   legend = l) %>%
    plotly::config(displaylogo = FALSE,
                   displayModeBar = 'hover',
                   modeBarButtonsToRemove = c('lasso2d',
                                              'select2d',
                                              'toggleSpikelines',
                                              'hoverClosestCartesian',
                                              'hoverCompareCartesian'),
                   toImageButtonOptions = list(format = "png", filename = plot_fname))

}

#' Get arguments for gene boxplotly
#'
#' @param eset ExpressionSet object with \code{'adjusted'} assayDataElement
#' @param explore_genes Character vector of genes to plot
#' @param dataset_name Name of bulk dataset.
#'
#' @return List with items \code{'df'}, \code{'boxgap'}, \code{'boxgroupgap'}, and \code{'plot_fname'}.
#' @export
get_boxplotly_gene_args <- function(eset, explore_genes, dataset_name) {
  dat <- Biobase::assayDataElement(eset, 'adjusted')
  pdata <- Biobase::pData(eset)

  dfs <- list()
  for (gene in explore_genes) {
    dfs[[gene]] <- data.frame(text = row.names(pdata),
                              x = gene,
                              y = dat[gene, row.names(pdata)],
                              name = as.character(pdata$`Group name`),
                              color = pdata$Group,
                              stringsAsFactors = FALSE)
  }

  df <- do.call(rbind, dfs)

  group_levels <- as.character(sort(unique(df$color)))

  # adjust gaps within/between group based on number of boxs (chosen by trial and error)
  nbox <- length(group_levels) * length(explore_genes)
  boxgap <- ifelse(nbox > 5, 0.4, 0.6)
  boxgroupgap <- ifelse(nbox > 6, 0.3, 0.6)

  # plotly bug when two groups uses first and third color in RColorBrewer Set2 pallette
  if (length(group_levels) <= 2)
    group_levels <- c(group_levels, 'NA')

  df$color <- factor(df$color, levels = group_levels)
  df$x <- factor(df$x, levels = explore_genes)

  # name for saving plot
  fname <- paste(explore_genes, collapse = '_')
  fname <- paste('bulk', dataset_name, fname, Sys.Date(), sep='_')

  return(list(df = df,
              boxgap = boxgap,
              boxgroupgap = boxgroupgap,
              plot_fname = fname))
}


#' Get arguments for cells boxplotly
#'
#' @param pdata data.frame of phenotype data.
#' @param dtangle_est data.frame of proportion estimates from dtangle.
#' @param dataset_name Name of bulk dataset.
#'
#' @return List with items \code{'df'}, \code{'boxgap'}, \code{'boxgroupgap'}, and \code{'plot_fname'}.
#' @export
get_boxplotly_cell_args <- function(pdata, dtangle_est, dataset_name) {

  common <- intersect(row.names(dtangle_est), row.names(pdata))
  dtangle_est <- as.data.frame(dtangle_est)[common, ]
  pdata <- pdata[common, c('Title', 'Group name', 'Group')]
  colnames(pdata) <- c('text', 'name', 'color')

  df <- stack(dtangle_est)
  colnames(df) <- c('y', 'x')

  df <- cbind(pdata, df)
  group_levels <- as.character(sort(unique(df$color)))
  clus_levels <- colnames(dtangle_est)

  # adjust gaps within/between group based on number of boxs (chosen by trial and error)
  nbox <- length(group_levels) * length(clus_levels)
  boxgap <- ifelse(nbox > 5, 0.4, 0.6)
  boxgroupgap <- ifelse(nbox > 6, 0.3, 0.6)

  # plotly bug when two groups uses first and third color in RColorBrewer Set2 pallette
  if (length(group_levels) <= 2)
    group_levels <- c(group_levels, 'NA')


  df$color <- factor(df$color, levels = group_levels)
  df$x <- factor(df$x, levels = clus_levels)

  # name for saving plot
  fname <- paste(clus_levels, collapse = '_')
  fname <- gsub(' ', '', fname)
  fname <- paste(dataset_name, fname, Sys.Date(), sep='_')

  return(list(df = df,
              boxgap = boxgap,
              boxgroupgap = boxgroupgap,
              plot_fname = fname))
}



#' Load previous bulk anals dataframe
#'
#' used for determining available anals for download
#'
#' @param data_dir Path to folder container \code{'bulk'} and \code{'single-cell'} directories
#' @return data.frame with columns "dataset_name" "dataset_dir" and "anal_name".
#' @export
#' @keywords internal
load_bulk_anals <- function(data_dir) {
  anals_path <- file.path(data_dir, 'bulk', 'anals.rds')

  if (file.exists(anals_path)) {
    anals <- readRDS(anals_path)

  } else {
    anals <- data.frame(matrix(ncol = 3, nrow = 0), stringsAsFactors = FALSE)
    colnames(anals) <- c("dataset_name", "dataset_dir", "anal_name")
    saveRDS(anals, anals_path)
  }

  anals$label <- anals$anal_name
  anals$value <- seq_len(nrow(anals))


  return(anals)
}

#' Save new analysis info to anals dataframe
#'
#' @param dataset_name Name of dataset
#' @param dataset_dir Folder name for dataset
#' @param anal_name Name of new analysis
#' @param data_dir Path to folder container \code{'bulk'} and \code{'single-cell'} directories
#' @return NULL
#' @export
#' @keywords internal
save_bulk_anals <- function(dataset_name, dataset_dir, anal_name, data_dir) {
  anals_path <- file.path(data_dir, 'bulk', 'anals.rds')
  anals <- readRDS(anals_path)

  anals[nrow(anals)+1, ] <- c(dataset_name, dataset_dir, anal_name)
  anals <- anals[!duplicated(anals), ]
  saveRDS(anals, anals_path)
}

#' Remove analysis info from anals dataframe
#'
#' Used to clear previous analysis when groupings change
#'
#' @inheritParams save_bulk_anals
#'
#' @return NULL
#' @export
#' @keywords internal
remove_bulk_anals <- function(dataset_name, data_dir) {
  anals_path <- file.path(data_dir, 'bulk', 'anals.rds')
  anals <- readRDS(anals_path)
  is.ds <- anals$dataset_name %in% dataset_name
  anals <- anals[!is.ds, ]
  saveRDS(anals, anals_path)
}

#' Load previous bulk datasets dataframe
#'
#' @param data_dir Path to folder container \code{'bulk'} and \code{'single-cell'} directories
#' @return data.frame with columns "dataset_name" and "dataset_dir"
#' @export
#' @keywords internal
load_bulk_datasets <-function(data_dir) {
  datasets_path <- file.path(data_dir, 'bulk', 'datasets.rds')

  datasets <- data.frame(matrix(ncol = 2, nrow = 0), stringsAsFactors = FALSE)
  colnames(datasets) <- c("dataset_name", "dataset_dir")

  dataset_names <- list.dirs(file.path(data_dir, 'bulk'), full.names = FALSE, recursive = FALSE)
  has.eset <- file.exists(file.path(data_dir, 'bulk', dataset_names, 'eset.rds'))
  dataset_names <- dataset_names[has.eset]

  datasets <- data.frame(dataset_name = dataset_names,
                         dataset_dir = file.path('bulk', dataset_names), stringsAsFactors = FALSE)

  datasets$value <-  datasets$label <- datasets$dataset_name
  if (nrow(datasets)) datasets$type <- 'Bulk Data'

  return(datasets)
}


#' Delete stale dataset files after changing sample groups or running sva
#'
#' @param data_dir Path to folder with dataset files
#' @param patterns patterns to remove. Default is all.
#'
#' @export
#' @keywords internal
remove_dataset_files <- function(data_dir, patterns = c('^adjusted_\\d+svs.rds$',
                                                        '^iqr_keep_\\d+svs.rds$',
                                                        '^vsd.rds$',
                                                        '^svobj.rds$',
                                                        '^numsv.rds$',
                                                        '^lm_fit_\\d+svs.rds$',
                                                        'cmap_res_.+_\\d+svs.rds$',
                                                        'go_.+_\\d+svs.rds$',
                                                        'kegg_.+_\\d+svs.rds$',
                                                        'l1000_drugs_res_.+_\\d+svs.rds$',
                                                        'l1000_genes_res_.+_\\d+svs.rds$',
                                                        'diff_expr_symbol_.+_\\d+svs.rds$'), exclude = NULL) {
  for (pattern in patterns) {
    fpaths <- list.files(data_dir, pattern)
    if (!is.null(exclude)) {
      keep <- grepl(exclude, fpaths)
      fpaths <- fpaths[!keep]
    }
    unlink(file.path(data_dir, fpaths))
  }
}


