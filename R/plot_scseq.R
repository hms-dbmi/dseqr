#' Plot TSNE coloured by cluster
#'
#'
#' @return \code{ggplot}
#' @export
plot_tsne_cluster <- function(scseq, selected_clusters = levels(scseq$cluster), legend_title = 'Cluster', cols = NULL) {

  SingleCellExperiment::logcounts(scseq) <- as(SingleCellExperiment::logcounts(scseq), 'dgCMatrix')
  scseq <- Seurat::as.Seurat(scseq, counts = NULL)
  Idents(scseq) <- scseq$cluster

  if (is.null(cols)) cols <- get_palette(levels(scseq$cluster))

  # make selected cluster and groups stand out
  cols <- ggplot2::alpha(cols, alpha = ifelse(levels(scseq$cluster) %in% selected_clusters, 1, 0.1))
  pt.size <- min(6000/ncol(scseq), 2)

  cluster_plot <- Seurat::DimPlot(scseq, reduction = 'TSNE', cols = cols, pt.size = pt.size, label = TRUE, label.size = 6, repel = TRUE) +
    theme_no_axis_vals() +
    ggplot2::xlab('TSNE1') +
    ggplot2::ylab('TSNE2') +
    theme_dimgray(with_nums = FALSE) +
    ggplot2::theme(legend.position = 'none', text = ggplot2::element_text(color = 'dimgray'))

  return(cluster_plot)
}

#' Plot UMAP coloured by HGNC symbol
#'
#' @param scseq \code{SingleCellExperiment} object.
#' @param gene Character vector specifying gene to colour cells by.
#' @param selected_idents The groups in \code{scseq$orig.ident} to show cell for. The default \code{NULL} shows all cells.
#' @param pt.size Numeric scalar, specifying the size of the points. Defaults to 3.
#'
#' @return \code{ggplot}
#' @export
plot_tsne_gene <- function(scseq, gene, selected_idents = unique(scseq$project)) {

  SingleCellExperiment::logcounts(scseq) <- as(SingleCellExperiment::logcounts(scseq), 'dgCMatrix')
  scseq <- Seurat::as.Seurat(scseq, counts = NULL)
  Idents(scseq) <- scseq$cluster

  # make selected groups stand out
  cells <- colnames(scseq)
  cells <- cells[scseq$project %in% selected_idents]

  pt.size <- min(6000/ncol(scseq), 2)

  gene_plot <- Seurat::FeaturePlot(scseq, gene, reduction = 'TSNE', pt.size = pt.size, order = TRUE) +
    theme_no_axis_vals() +
    ggplot2::xlab('TSNE1') +
    ggplot2::ylab('TSNE2') +
    ggplot2::theme(plot.title = ggplot2::element_blank()) +
    theme_dimgray(with_nums = FALSE) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour = 'black'))

  # need until bug fixed in Seurat, then use cells argument in FeaturePlot
  gene_plot$data <- gene_plot$data[row.names(gene_plot$data) %in% cells, ]
  gene_plot$labels$colour <- gene

  return(gene_plot)
}

#' Format gene plots for sample comparison for drugseqr app
#'
#' @param plot Returned by
#' @param group Level in \code{scseq$orig.ident} to show cells for. Either \code{'ctrl'} or \code{'test'}
#' @param scseq \code{Seurat} object.
#'
#' @return \code{plot} formatted for drugseqr app
#' @export
#' @keywords internal
format_sample_gene_plot <- function(plot, group, selected_gene, scseq) {

  selected_gene <- make.names(selected_gene)

  # the min and max gene expression value
  lims <- range(plot$data[[selected_gene]])

  # show selected group only
  sel.cells <- colnames(scseq)[scseq$orig.ident == group]
  plot$data <- plot$data[row.names(plot$data) %in% sel.cells, ]

  # add selected group as title
  plot <- plot + ggplot2::ggtitle(toupper(group)) +
    ggplot2::theme(plot.title = ggplot2::element_text(color = 'black'))

  # use the same scale for each sample so that comparable
  suppressMessages(plot <- plot + ggplot2::scale_color_continuous(low ="lightgray", high = "blue", limits = lims))

  # remove control plot labels and legend
  if (group == 'ctrl')
    plot <- plot + ggplot2::xlab('') + ggplot2::ylab('') +
    ggplot2::theme(legend.position = 'none')

  return(plot)
}

theme_no_axis_vals <- function() {
  ggplot2::theme(axis.text = ggplot2::element_blank(),
                 axis.ticks = ggplot2::element_blank(),
                 legend.text = ggplot2::element_text(size = 14, margin = ggplot2::margin(b=7)),
                 legend.title = ggplot2::element_text(size = 16, margin = ggplot2::margin(b=9)))
}


#' Plot BioGPS data for a HGNC symbol
#'
#' @param gene Character vector of gene name.
#' @keywords internal
#'
#' @return ggplot
#' @export
plot_biogps <- function(gene) {
  if (!length(gene) || !gene %in% biogps[, SYMBOL]) return(NULL)

  gene_dat <- unlist(biogps[gene, -c('ENTREZID', 'SYMBOL')])
  gene_dat <- sort(gene_dat, decreasing = TRUE)[1:15]
  gene_dat <- tibble::tibble(mean = gene_dat,
                             source = factor(names(gene_dat), levels = rev(names(gene_dat))))

  ggplot2::ggplot(gene_dat, ggplot2::aes(x = source, y = mean, fill = mean)) +
    ggplot2::theme_minimal() +
    ggplot2::geom_bar(stat = "identity", color = 'black', size = 0.1, width = 0.7) +
    ggplot2::scale_fill_distiller(palette = 'Reds', name = '', direction = 1) +
    ggplot2::xlab('') +
    ggplot2::ylab('') +
    ggplot2::ggtitle('BioGPS Human Gene Atlas Expression') +
    ggplot2::scale_y_continuous(expand = c(0, 0)) + # Set the axes to cross at 0
    ggplot2::coord_flip() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 16),
                   axis.text = ggplot2::element_text(size = 14),
                   axis.text.x = ggplot2::element_blank(),
                   legend.position = "none")
}


#' Get a pallete for cluster plots
#'
#' @param levs Character vector of levels to get colour pallete for.
#'
#' @return Character vector with colour codes of \code{length(levs)}.
#' @export
#' @keywords internal
get_palette <- function(levs) {

  # palettes from scater
  tableau20 <- c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
                 "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
                 "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
                 "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

  tableau10medium <- c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                       "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                       "#CDCC5D", "#6DCCDA")


  nlevs <- length(levs)
  if (nlevs <= 10) {
    values <- head(tableau10medium, nlevs)

  } else if (nlevs <= 20) {
    values <- head(tableau20, nlevs)

  } else {
    set.seed(0)
    values <- randomcoloR::distinctColorPalette(nlevs)

  }
  return(values)
}


