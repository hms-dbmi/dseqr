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
  ggplot2::ggplot(sce_df, aes(x=get(x), fill=whitelist)) +
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

  suppressMessages(tsne_white <- scater::plotTSNE(sce, colour_by='whitelist', point_alpha = 1) +
                     ggplot2::scale_fill_manual(values = c("#E31A1C", "#FFFFCC")) +
                     guides(fill=guide_legend(title='whitelist')) +
                     theme(legend.position = 'top',
                           axis.text=element_blank(),
                           axis.ticks=element_blank()) +
                     xlab('') +
                     ylab(''))

  # make really low mito/ribo content more prominent as well
  mito_scale <- get_diverge(sce$pct_counts_mito, 'Mitochondrial percent')
  ribo_scale <- get_diverge(sce$pct_counts_ribo, 'Ribosomal percent')

  suppressMessages(tsne_mito <- scater::plotTSNE(sce, colour_by='pct_counts_mito', point_alpha = 1) +
                     mito_scale +
                     theme(legend.position = 'top',
                           axis.text=element_blank(),
                           axis.ticks=element_blank()) +
                     xlab(''))

  suppressMessages(tsne_ribo <- scater::plotTSNE(sce, colour_by='pct_counts_ribo', point_alpha = 1) +
                     ribo_scale +
                     theme(legend.position = 'top',
                           axis.title.x=element_blank(),
                           axis.text=element_blank(),
                           axis.ticks=element_blank()) +
                     xlab('') +
                     ylab(''))

  suppressMessages(tsne_ncounts <- scater::plotTSNE(sce, colour_by='log10_total_counts', point_alpha = 1) +
                     ggplot2::scale_fill_distiller(palette = 'YlOrRd', name = 'Log of total counts') +
                     theme(legend.position = 'top',
                           axis.text=element_blank(),
                           axis.ticks=element_blank()) +
                     ylab(''))


  suppressMessages(tsne_nexp <- scater::plotTSNE(sce, colour_by='log10_total_features_by_counts', point_alpha = 1) +
                     ggplot2::scale_fill_distiller(palette = 'YlOrRd', name = 'Log of expressed genes') +
                     theme(legend.position = 'top',
                           axis.text=element_blank(),
                           axis.ticks=element_blank()) +
                     xlab('') +
                     ylab(''))

  suppressWarnings(cowplot::plot_grid(tsne_white, NULL, tsne_mito, tsne_ribo, tsne_ncounts, tsne_nexp, nrow = 3))
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
    colors = rev(brewer.pal(9, "Spectral")),
    values = range01(c(seq(qx[1], qx[2], length.out = 25),
                       seq(qx[2], qx[10], length.out = 25),
                       seq(qx[10], qx[11], length.out = 25))))

  return(diverge)
}

