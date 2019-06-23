#' Plot UMAP coloured by cluster
#'
#' @param ... Additional arguments passed to \code{Seurat::DimPlot}.
#' @inheritParams plot_umap_gene
#'
#' @return \code{ggplot}
#' @export
#'
#' @examples
plot_umap_cluster <- function(scseq, selected_groups = NULL, pt.size = 3) {

  cols <- get_colour_values(levels(scseq$seurat_clusters))

  # make selected cluster and groups stand out
  if (!is.null(selected_groups))
    cols <- ggplot2::alpha(cols, alpha = ifelse(levels(scseq$orig.ident) %in% selected_groups, 1, 0.1))

  Seurat::DimPlot(scseq, reduction = 'umap', cols = cols, pt.size = pt.size) +
    ggplot2::xlab('UMAP1') +
    ggplot2::ylab('UMAP2') +
    ggplot2::guides(colour=ggplot2::guide_legend(title='Cluster')) +
    ggplot2::theme(plot.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())

}


#' Plot UMAP coloured by HGNC symbol
#'
#' @param scseq \code{Seurat}.
#' @param gene Character vector specifying gene to colour cells by.
#' @param selected_groups The groups in \code{scseq$orig.ident} to show cell for. The default \code{NULL} shows all cells.
#' @param pt.size Numeric scalar, specifying the size of the points. Defaults to 3.
#' @inheritParams explore_scseq_clusters
#'
#' @return \code{ggplot}
#' @export
#'
#' @examples
plot_umap_gene <- function(scseq, gene, selected_groups = NULL, pt.size = 3) {

  # make selected groups stand out
  cells <- NULL
  if (!is.null(cells))
    cells <- colnames(scseq)[scseq$orig.ident %in% selected_groups]

  gene_plot <- Seurat::FeaturePlot(scseq,
                                   gene,
                                   cells = cells,
                                   cols = c('lightgray', 'red'),
                                   reduction = 'umap',
                                   order = TRUE,
                                   pt.size = pt.size) +
    ggplot2::xlab('UMAP1') +
    ggplot2::ylab('UMAP2') +
    ggplot2::theme(plot.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())

  gene_plot$labels$colour <- gene
  return(gene_plot)

}


#' Plot BioGPS data for a HGNC symbol
#'
#' @param gene Character vector of gene name.
#' @keywords internal
#'
#' @return ggplot
#' @export
#'
#' @examples
plot_biogps <- function(gene) {
  if (!length(gene) || !gene %in% biogps[, SYMBOL]) return(NULL)

  gene_dat <- unlist(biogps[gene, -c('ENTREZID', 'SYMBOL')])
  gene_dat <- sort(gene_dat, decreasing = TRUE)[1:20]
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
                   plot.title = ggplot2::element_text(),
                   axis.text = ggplot2::element_text(size = 12),
                   axis.text.x = ggplot2::element_blank(),
                   legend.position = "none")
}


#' Remove ggplot xaxis title, text, and ticks
#'
#' @return \code{theme}
#' @export
#' @keywords internal
#'
#' @examples
theme_no_xaxis <- function() {
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank())
}

#' Remove ggplot yaxis title, text, and ticks
#'
#' @return \code{theme}
#' @export
#' @keywords internal
#'
#' @examples
theme_no_yaxis <- function() {
  ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank())
}

#' Make ggplot axes and text dimgray
#'
#' @return \code{theme}
#' @export
#' @keywords internal
#'
#' @examples
theme_dimgray <- function() {
  ggplot2::theme(axis.line.y = ggplot2::element_line(size = 0.1, color = 'dimgray'),
                 axis.line.x = ggplot2::element_line(size = 0.1, color = 'dimgray'),
                 axis.ticks.x = ggplot2::element_line(size = 0.1, color = 'dimgray'),
                 axis.ticks.y = ggplot2::element_line(size = 0.1, color = 'dimgray'),
                 axis.text = ggplot2::element_text(colour = 'dimgray'),
                 axis.title = ggplot2::element_text(colour = 'dimgray'),
                 text = ggplot2::element_text(colour = 'dimgray'))
}
