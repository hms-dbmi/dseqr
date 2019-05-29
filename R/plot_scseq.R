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
hist_scseq <- function(sce) {

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
