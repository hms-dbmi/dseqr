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
#' @param sce \code{SingleCellExperiment} returned by \code{load_scseq}
#'
#' @return plot
#' @export
#'
#' @examples
hist_scseq_whitelist <- function(sce) {

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
#' @param sce
#'
#' @return
#' @export
#'
#' @examples
tsne_scseq_whitelist <- function(sce) {

  # make sure have qc metrics
  sce <- qc_scseq(sce)

  # uses denoised expression matrix
  set.seed(1000)
  sce <- scater::runTSNE(sce, use_dimred="PCA")

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

