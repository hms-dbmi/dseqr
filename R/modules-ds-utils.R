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

#' Generate boxplotly for vsd normalized gene data by group for Datasets tab
#'
#' @param eset ExpressionSet object with \code{'vsd'} assayDataElement
#' @param explore_genes Character vector of genes to plot
#' @param pdata data.frame with columns \code{'Group name'} (character) and \code{'Group'} (integer).
#'
#' @return plotly
#' @export
#'
plotlyGene <- function(eset, explore_genes, pdata, dataset_name) {

  dat <- Biobase::assayDataElement(eset, 'vsd')

  dfs <- list()
  for (gene in explore_genes) {
    dfs[[gene]] <- data.frame(Sample = row.names(pdata),
                              Gene = gene,
                              Expression = dat[gene, row.names(pdata)],
                              Name = as.character(pdata$`Group name`),
                              Num = pdata$Group,
                              stringsAsFactors = FALSE)

  }

  df <- do.call(rbind, dfs)
  group_levels <- as.character(sort(unique(df$Num)))

  # adjust gaps within/between group based on number of boxs (chosen by trial and error)
  nbox <- length(group_levels) * length(explore_genes)
  boxgap <- ifelse(nbox > 5, 0.4, 0.6)
  boxgroupgap <- ifelse(nbox > 6, 0.3, 0.6)

  # plotly bug when two groups uses first and third color in RColorBrewer Set2 pallette
  if (length(group_levels) <= 2)
    group_levels <- c(group_levels, 'NA')

  df$Num <- factor(df$Num, levels = group_levels)
  df$Gene <- factor(df$Gene, levels = explore_genes)

  l <- list(
    font = list(
      family = "sans-serif",
      size = 12,
      color = "#000"),
    bgcolor = "#f8f8f8",
    bordercolor = "#e7e7e7",
    borderwidth = 1)

  # name for saving plot
  fname <- paste(explore_genes, collapse = '_')
  fname <- paste('bulk', dataset_name, fname, Sys.Date(), sep='_')

  df %>%
    plotly::plot_ly() %>%
    plotly::add_trace(x = ~ Gene,
                      y = ~Expression,
                      color = ~Num,
                      text = ~Sample,
                      name = ~Name,
                      type = 'box',
                      boxpoints = 'all',
                      jitter = 0.8,
                      pointpos = 0,
                      fillcolor = 'transparent',
                      hoverinfo = 'text',
                      whiskerwidth = 0.1,
                      hoveron = 'points',
                      marker = list(color = "rgba(0, 0, 0, 0.6)")) %>%
    plotly::layout(boxmode = 'group', boxgroupgap = boxgroupgap, boxgap = boxgap,
                   xaxis = list(fixedrange=TRUE),
                   yaxis = list(fixedrange=TRUE, title = 'Normalized Expression'),
                   legend = l) %>%
    plotly::config(displaylogo = FALSE,
                   displayModeBar = 'hover',
                   modeBarButtonsToRemove = c('lasso2d',
                                              'select2d',
                                              'toggleSpikelines',
                                              'hoverClosestCartesian',
                                              'hoverCompareCartesian'),
                   toImageButtonOptions = list(format = "svg", filename = fname))

}



#' Load previous bulk anals dataframe
#'
#' @param data_dir Path to folder container \code{'bulk'} and \code{'single-cell'} directories
#' @return data.frame with columns "dataset_name" "dataset_dir" and "anal_name".
#' @export
#' @keywords internal
load_bulk_anals <- function(data_dir, with_type = FALSE) {
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

  if (with_type) {
    anals$type <- anals$dataset_name
    if (nrow(anals)) anals$type <- paste0('Bulk - ', anals$type)

  }

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

  if (file.exists(datasets_path)) {
    datasets <- readRDS(datasets_path)

  } else {
    datasets <- data.frame(matrix(ncol = 2, nrow = 0), stringsAsFactors = FALSE)
    colnames(datasets) <- c("dataset_name", "dataset_dir")
    saveRDS(datasets, datasets_path)
  }

  datasets$value <-  datasets$label <- datasets$dataset_name
  if (nrow(datasets)) datasets$type <- 'Bulk'

  return(datasets)
}


#' Save new dataset info to datasets dataframe
#'
#' @param dataset_name Name of dataset
#' @param dataset_dir Folder name for dataset
#' @param data_dir Path to folder container \code{'bulk'} and \code{'single-cell'} directories
#' @return NULL
#' @export
#' @keywords internal
save_bulk_dataset <- function(dataset_name, dataset_dir, data_dir) {
  datasets_path <- file.path(data_dir, 'bulk', 'datasets.rds')
  datasets <- readRDS(datasets_path)

  datasets[nrow(datasets)+1, ] <- c(dataset_name, dataset_dir)
  saveRDS(datasets, datasets_path)
}

