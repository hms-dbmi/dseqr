#' Plot TSNE coloured by cluster
#'
#' The groups are either \code{scseq$cluster} (default) or \code{scseq@metadata$colour_by}.
#'
#' @param legend_title Title for the legend. The default is \code{scseq@metadata$colour_by} with title case.
#' @inheritParams plot_tsne_gene
#'
#' @return \code{ggplot}
#' @export
#'
#' @examples
plot_tsne_cluster <- function(scseq, legend_title = NULL, selected_groups = NULL, point_size = 3, assay.type = 'logcounts') {
  if (class(scseq) == 'Seurat') scseq <- srt_to_sce(scseq, 'SCT')

  # make selected cluster and groups stand out
  point_alpha <- rep(1, ncol(scseq))
  if (!is.null(selected_groups)) point_alpha[!scseq$orig.ident %in% selected_groups] <- 0.1

  colour_by <- ifelse(is.null(scseq@metadata$colour_by), 'cluster', scseq@metadata$colour_by)
  legend_title <- ifelse(is.null(legend_title), tools::toTitleCase(colour_by), legend_title)

  cluster_plot <- scater::plotTSNE(scseq,
                                   by_exprs_values = assay.type,
                                   colour_by = colour_by,
                                   point_size = point_size,
                                   point_alpha = point_alpha,
                                   theme_size = 14) +
    ggplot2::guides(fill = ggplot2::guide_legend(legend_title))

  return(cluster_plot)
}

#' Plot TSNE coloured by HGNC symbol
#'
#' @param scseq \code{SingleCellExperiment} or \code{Seurat}.
#' @param gene Character vector specifying gene to colour cells by.
#' @param selected_groups The groups in \code{scseq$orig.ident} to show cell for. The default \code{NULL} shows all cells.
#' @param hide_mask A boolean vector indicating samples to hide. This is useful if you need to compare e.g. one sample to
#' another but want to maintain the same colour range for expression values. Passing in subsets of \code{scseq} would not
#' achieve this.
#' @param point_size Numeric scalar, specifying the size of the points. Defaults to 3.
#' @inheritParams explore_scseq_clusters
#'
#' @return \code{ggplot}
#' @export
#'
#' @examples
plot_tsne_gene <- function(scseq, gene, selected_groups = NULL, hide_mask = NULL, assay.type = 'logcounts', point_size = 3) {
  if (class(scseq) == 'Seurat') scseq <- srt_to_sce(scseq, 'SCT')

  # make selected groups stand out
  point_alpha <- rep(1, ncol(scseq))
  if (!is.null(selected_groups)) point_alpha[!scseq$orig.ident %in% selected_groups] <- 0
  if (!is.null(hide_mask)) point_alpha[hide_mask] <- 0

  # make sure not all zero for gene (otherwise all black)
  scseq@assays$data[[assay.type]][gene, ] <- scseq@assays$data$logcounts[gene, ] +
    runif(ncol(scseq), min=0.00001, max=0.0001)

  suppressMessages(scater::plotTSNE(scseq,
                                    by_exprs_values = assay.type,
                                    colour_by = gene,
                                    point_size = point_size,
                                    point_alpha = point_alpha,
                                    theme_size = 14) +
                     ggplot2::scale_fill_distiller(palette = 'Reds', name = gene, direction = 1))

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
