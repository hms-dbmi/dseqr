#' Save PDF of scRNA-seq reports
#'
#' @param fname String giving the name of the PDF file to save.
#' @inheritParams explore_scseq_clusters
#' @inheritParams plot_scseq_report
#'
#' @return Saves report to \code{fname}.
#' @export
#'
#' @examples
save_scseq_reports <- function(scseq, markers, fname, point_size = 3) {
  pdf(file = fname, paper = 'US', width = 8.50, height = 11.0, title = 'cluster markers')

  for (i in seq_along(markers)) {
    plot(plot_scseq_report(scseq, markers[i], point_size = point_size))
  }
  dev.off()
}


#' Plot grid with scRNA-seq cluster and gene markers together
#'
#' @param scseq \code{SingleCellExperiment} or \code{Seurat}.
#' @param point_size Numeric scalar, specifying the size of the points. Defaults to 3.
#' @inheritParams explore_scseq_clusters
#'
#' @return \code{ggplot}
#' @export
#'
#' @examples
plot_scseq_report <- function(scseq, markers, point_size = 3) {
  if (class(scseq) == 'Seurat') scseq <- srt_to_sce(scseq, 'SCT')
  selected_group <- names(markers)
  genes <- markers[[1]]

  # make orig.ident the clusters so that can highlight
  scseq$orig.ident <- scseq$cluster

  # get grid of gene plots
  gene_plots <- lapply(genes, function(gene) plot_tsne_gene(scseq, gene, point_size = point_size) +
                         theme_dimgray() + theme_no_xaxis() + theme_no_yaxis() +
                         ggplot2::theme(legend.position = 'none', plot.title=ggplot2::element_text(size=12, hjust = 0)) +
                         ggplot2::ggtitle(gene))

  gene_grid <- cowplot::plot_grid(plotlist = gene_plots, align = 'vh', ncol = 2)

  # get cluster plot
  cluster_plot <- plot_tsne_cluster(scseq, legend_title = 'Cell Type', selected_groups = selected_group, point_size = point_size) +
    ggplot2::theme(plot.margin=ggplot2::unit(c(20, 5.5, 20, 5.5), "points")) +
    theme_dimgray() + theme_no_xaxis() + theme_no_yaxis()

  # get the title
  # include number of cells and percentage of total
  ncells <- sum(scseq$cluster == selected_group)
  pcells <- round(ncells/ncol(scseq) * 100)

  title <- cowplot::ggdraw() +
    cowplot::draw_label(selected_group, x = 0, hjust = 0)

  label <- cowplot::ggdraw() +
    cowplot::draw_label('TSNE PLOTS', x=1, hjust = 1,  colour = 'dimgray', size = 12)

  title <- cowplot::plot_grid(title, label)

  subtitle <- cowplot::ggdraw() +
    cowplot::draw_label(paste0(ncells, ' cells - ', pcells, '%'), x = 0, hjust = 0, size = 12, colour = 'dimgray')

  # combine titles and plots
  marker_nrows <- ceiling(length(genes)/2)
  cowplot::plot_grid(title, subtitle, cluster_plot, gene_grid, ncol = 1, rel_heights = c(0.1, 0.1, 1, marker_nrows))
}
