#' Generate single histogram for whitelist factors
#'
#' @param sce_df \code{tibble} with columns \code{'whitelist'} and \code{get(x)}.
#' @param x Character vector indicating column in \code{sce_df} to use as x aesthetic.
#' @param xlab Label for x-axis
#' @param ylab Label for y-axis
#'
#' @return
#' @export
#'
#' @examples
geom_hist_whitelist <- function(sce_df, x, xlab = '', ylab = '') {
  ggplot2::ggplot(sce_df, ggplot2::aes(x=get(x), fill=whitelist)) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
    ggplot2::geom_histogram(position="identity", colour="black", alpha = 0.5, size = 0.05) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::theme(legend.position = "none")
}

#' Plot histograms for whitelist factors
#'
#' Currently plots histograms for percent ribosomal reads, percent mitochondrial reads, log of total counts,
#' and log of total expressed features.
#'
#' @param scseq \code{SingleCellExperiment} or \code{Seurat} returned by \code{load_scseq}
#'
#' @return plot
#' @export
#'
#' @examples
hist_scseq_whitelist <- function(scseq) {

  if (class(scseq) == 'Seurat') sce <- srt_to_sce(scseq)

  # make sure have qc metrics
  sce <- qc_scseq(sce)

  # construct tibble for ggploting
  sce_df <- tibble::tibble(log10_total_counts = sce$log10_total_counts,
                           log10_total_features_by_counts = sce$log10_total_features_by_counts,
                           pct_counts_ribo = sce$pct_counts_ribo,
                           pct_counts_mito = sce$pct_counts_mito,
                           whitelist = sce$whitelist)

  ribo <- geom_hist_whitelist(sce_df, 'pct_counts_ribo', 'Percent ribosomal', 'count')
  mito <- geom_hist_whitelist(sce_df, 'pct_counts_mito', 'Percent mitochondrial')
  ncnt <- geom_hist_whitelist(sce_df, 'log10_total_counts', 'Log total counts', 'count')
  nexp <- geom_hist_whitelist(sce_df, 'log10_total_features_by_counts', 'Log expressed genes')

  suppressMessages(legend <- cowplot::get_legend(nexp + ggplot2::theme(legend.position="top")))
  suppressMessages(cowplot::plot_grid(legend, NULL, ncnt, nexp, ribo, mito, nrow = 3, rel_heights = c(0.2, 1, 1)))
}


#' Generates individual tsne plots for tsne_scseq_whitelist
#'
#' @param sce \code{SingleCellExperiment}
#' @param colour_by List within \code{sce} to color by.
#' @param name The name of the legend.
#' @param xlab Character, x-axis label. Default is \code{''}.
#' @param ylab Character, y-axis label. Default is \code{''}.
#' @param scale One if either \code{'distiller'} (default), \code{'diverge'} (for mito/ribo), or \code{'manual'} (for whitelist).
#'
#' @return ggplot2 plot
#' @export
#'
#' @examples
geom_tsne <- function(sce, colour_by, name = colour_by, xlab = '', ylab = '', scale = 'distiller') {

  # make really low mito/ribo content more prominent as well
  if (scale == 'diverge')
    colorscale <- get_diverge(sce[[colour_by]], name)

  if (scale == 'distiller')
    colorscale <- ggplot2::scale_fill_distiller(palette = 'YlOrRd', name = name)

  if (scale == 'manual')
    colorscale <- ggplot2::scale_fill_manual(values = c("#E31A1C", "#FFFFCC"), name = name)

  suppressMessages(scater::plotTSNE(sce, colour_by=colour_by, point_alpha = 1) +
                     colorscale +
                     ggplot2::xlab(xlab) +
                     ggplot2::ylab(ylab) +
                     ggplot2::theme(legend.position = 'top',
                                    axis.text=ggplot2::element_blank(),
                                    axis.ticks=ggplot2::element_blank()))
}

#' tSNE plots of whitelisted cells
#'
#' @param scseq \code{SingleCellExperiment} or \code{Seurat} returned by \code{load_scseq}
#'
#' @return
#' @export
#'
#' @examples
tsne_scseq_whitelist <- function(scseq) {

  # SCTransform from Seurat: stabilized counts end up in "logcounts" (default for runTSNE)
  # SingleCellExperiment: uses denoised PCA (aka variance stabilized) expression matrix

  if (class(scseq) == 'Seurat') {
    sce <- srt_to_sce(scseq)
    use_dimred <- NULL

  } else {
    sce <- scseq
    use_dimred <- 'PCA'
  }

  # make sure have qc metrics
  sce <- qc_scseq(sce)

  set.seed(1000)
  sce <- scater::runTSNE(sce, use_dimred=use_dimred)

  tsne_white <- geom_tsne(sce, 'whitelist', scale = 'manual')
  tsne_mito <- geom_tsne(sce, 'pct_counts_mito', 'Mitochondrial percent', ylab = 'Dimension 2', scale = 'diverge')
  tsne_ribo <- geom_tsne(sce, 'pct_counts_ribo', 'Ribosomal percent', scale = 'diverge')
  tsne_ncnt <- geom_tsne(sce, 'log10_total_counts', 'Log of total counts', xlab = 'Dimension 1')
  tsne_nexp <- geom_tsne(sce, 'log10_total_features_by_counts', 'Log of expressed genes')

  suppressWarnings(cowplot::plot_grid(tsne_white, NULL, tsne_mito, tsne_ribo, tsne_ncnt, tsne_nexp, nrow = 3))
}



#' Get a continuous color scale that emphasizes 0-10th and 90-100th percentiles.
#'
#' @param x Numeric vector
#' @param name Character, name of scale
#'
#' @return ggplot2 \code{ScaleContinuous}
#' @export
#'
#' @examples
get_diverge <- function(x, name) {
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  qx <- quantile(x, probs = seq(0, 1, 0.1))

  diverge <- ggplot2::scale_fill_gradientn(
    name = name,
    colors = rev(RColorBrewer::brewer.pal(9, "Spectral")),
    values = range01(c(seq(qx[1], qx[2], length.out = 25),
                       seq(qx[2], qx[10], length.out = 25),
                       seq(qx[10], qx[11], length.out = 25))))

  return(diverge)
}




#' Plot TSNE coloured by cluster
#'
#' The groups are either \code{sce$cluster} (default) or \code{sce@metadata$colour_by}.
#'
#' @inheritParams plot_tsne_gene
#' @param legend_title Title for the legend.
#'
#' @return ggplot
#' @export
#'
#' @examples
plot_tsne_cluster <- function(sce, legend_title, selected_groups = NULL, point_size = 3, assay.type = 'logcounts', remove_axes = 'xaxis') {
  if (class(sce) == 'Seurat') sce <- srt_to_sce(sce, 'SCT')

  # make selected cluster and groups stand out
  point_alpha <- rep(1, ncol(sce))
  if (!is.null(selected_groups)) point_alpha[!sce$orig.ident %in% selected_groups] <- 0.2

  colour_by <- ifelse(is.null(sce@metadata$colour_by), 'cluster', sce@metadata$colour_by)

  cluster_plot <- scater::plotTSNE(sce,
                                   by_exprs_values = assay.type,
                                   colour_by = colour_by,
                                   point_size = point_size,
                                   point_alpha = point_alpha,
                                   theme_size = 14) +
    ggplot2::guides(fill = ggplot2::guide_legend(legend_title)) +
    theme_minimal_axes()

  cluster_plot <- remove_axes(cluster_plot, type = remove_axes)
  return(cluster_plot)
}

#' Plot TSNE coloured by HGNC symbol
#'
#' @param sce \code{SingleCellExperiment}
#' @param gene Character vector to colour cells by
#' @param selected_groups The groups in \code{sce$orig.ident} to show cell for. The default \code{NULL} shows all cells.
#' @inheritParams explore_scseq_clusters
#'
#' @return ggplot
#' @export
#'
#' @examples
plot_tsne_gene <- function(sce, gene, selected_groups = NULL, assay.type = 'logcounts', point_size = 3) {
  if (class(sce) == 'Seurat') sce <- srt_to_sce(sce, 'SCT')

  # make selected groups stand out
  point_alpha <- rep(1, ncol(sce))
  if (!is.null(selected_groups)) point_alpha[!sce$orig.ident %in% selected_groups] <- 0


  suppressMessages(scater::plotTSNE(sce,
                                    by_exprs_values = assay.type,
                                    colour_by = gene,
                                    point_size = point_size,
                                    point_alpha = point_alpha,
                                    theme_size = 14) +
                     ggplot2::scale_fill_distiller(palette = 'Reds', name = gene, direction = 1))


}

#' Remove ggplot axis title, text, and ticks
#'
#' @param plot ggplot
#' @param types Character vector of axis types to remove. Either or both of \code{'xaxis'} and \code{'yaxis'}.
#'
#' @return ggplot
#' @export
#'
#' @examples
remove_axes <- function(plot, types) {

  if ('xaxis' %in% types)
    plot <- plot +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())

  if ('yaxis' %in% types)
    plot <- plot +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())

  return(plot)
}

#' Make ggplot axes and text dimgray
#'
#' @param plot ggplot
#'
#' @return ggplot
#' @export
#'
#' @examples
theme_dimgray <- function(plot) {
  ggplot2::theme(axis.line.y = ggplot2::element_line(size = 0.1, color = 'dimgray'),
                 axis.line.x = ggplot2::element_line(size = 0.1, color = 'dimgray'),
                 axis.ticks.x = ggplot2::element_line(size = 0.1, color = 'dimgray'),
                 axis.ticks.y = ggplot2::element_line(size = 0.1, color = 'dimgray'),
                 axis.text = ggplot2::element_text(colour = 'dimgray'),
                 axis.title = ggplot2::element_text(colour = 'dimgray'),
                 text = ggplot2::element_text(colour = 'dimgray'))
}


#' Plot grid with scRNA-seq cluster and gene markers together
#'
#' @param sce
#' @param markers
#' @param point_size
#'
#' @return
#' @export
#'
#' @examples
plot_scseq_report <- function(sce, markers, point_size = 3) {
  if (class(sce) == 'Seurat') sce <- srt_to_sce(sce, 'SCT')
  selected_group <- names(markers)
  genes <- markers[[1]]

  # make orig.ident the clusters so that can highlight
  sce$orig.ident <- sce$cluster

  # get grid of gene plots and a cluster plot
  gene_plots <- lapply(genes, function(gene) plot_tsne_gene(sce, gene, point_size = point_size) + theme_dimgray())
  gene_plots <- lapply(gene_plots, remove_axes, c('xaxis', 'yaxis'))

  gene_grid <- cowplot::plot_grid(plotlist = gene_plots, align = 'vh', ncol = 2)
  cluster_plot <- plot_tsne_cluster(sce, legend_title = 'Cell Type', selected_groups = selected_group, remove_axes = 'none', point_size = point_size) +
    ggplot2::theme(plot.margin=unit(c(20, 5.5, 20, 5.5), "points")) +
    theme_dimgray()

  # now add the title
  title <- cowplot::ggdraw() +
    cowplot::draw_label(selected_group, x = 0, hjust = 0)
    ggplot2::theme(plot.margin = margin(0, 0, 0, 7))

  # combine
  marker_nrows <- ceiling(length(genes)/2)
  cowplot::plot_grid(title, cluster_plot, gene_grid, ncol = 1, rel_heights = c(0.1, 1.5, marker_nrows))
}

#' Save PDF of scRNA-seq reports
#'
#' @param sce
#' @param markers
#' @param fname
#' @param point_size
#'
#' @return
#' @export
#'
#' @examples
save_scseq_reports <- function(sce, markers, fname, point_size = 3) {
  pdf(file = fname, paper = 'US', width = 8.50, height = 11.0, title = 'cluster markers')
  for (i in seq_along(markers)) {
    plot(plot_scseq_report(sce, markers[i], point_size = point_size))
  }
  dev.off()
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


