#' Save PDF of scRNA-seq reports
#'
#' @param scseq \code{Seurat}.
#' @param markers Named list of character vectors specifying the marker genes to plot for each cluster. The names are the
#'  clusters that reports will be generated for in order.
#' @param fname String giving the name of the PDF file to save.
#' @param pt.size Numeric scalar, specifying the size of the points. Defaults to 3.
#'
#' @return Saves report to \code{fname}.
#' @export
save_scseq_reports <- function(scseq, markers, fname, pt.size = 3) {

  if (!all(names(markers) %in% levels(scseq$seurat_clusters)))
    stop('All names of markers must also be in seurat cluster levels.')

  # exclude any clusters that markers are not provided for
  scseq <- scseq[, scseq$seurat_clusters %in% names(markers)]

  # put levels in order of markers
  scseq$seurat_clusters <- factor(scseq$seurat_clusters, levels = names(markers))
  Seurat::Idents(scseq) <- 'seurat_clusters'

  pdf(file = fname, paper = 'US', width = 8.50, height = 11.0, title = 'cluster markers')
  for (i in seq_along(markers)) {
    plot(plot_scseq_report(scseq, markers[i], pt.size = pt.size))
  }
  dev.off()
}


#' Plot grid with scRNA-seq cluster and gene markers together
#'
#' @param scseq \code{Seurat} object.
#' @param markers Named list with a single character vector specifying the marker genes to plot.
#'  The name is the cluster that will be highlighted.
#' @param point_size Numeric scalar, specifying the size of the points. Defaults to 3.
#'
#' @return \code{ggplot}
#' @export
#' @keywords internal
plot_scseq_report <- function(scseq, markers, pt.size = 3) {

  selected_cluster <- names(markers)
  genes <- markers[[1]]

  # make orig.ident the clusters so that can highlight
  scseq$orig.ident <- scseq$seurat_clusters

  # get grid of gene plots
  gene_plots <- lapply(genes, function(gene) plot_umap_gene(scseq, gene, pt.size = pt.size) +
                         theme_dimgray() +
                         theme_no_xaxis() +
                         theme_no_yaxis() +
                         ggplot2::theme(legend.position = 'none', plot.title=ggplot2::element_text(size=12, hjust = 0)) +
                         ggplot2::ggtitle(gene))

  gene_grid <- cowplot::plot_grid(plotlist = gene_plots, align = 'vh', ncol = 2)

  # get cluster plot
  cluster_plot <- plot_umap_cluster(scseq, selected_clusters = selected_cluster, pt.size = pt.size) +
    ggplot2::theme(plot.margin=ggplot2::unit(c(20, 5.5, 20, 5.5), "points")) +
    theme_dimgray() + theme_no_xaxis() + theme_no_yaxis()

  # get the title
  # include number of cells and percentage of total
  ncells <- sum(scseq$seurat_clusters == selected_cluster)
  pcells <- round(ncells/ncol(scseq) * 100)

  title <- cowplot::ggdraw() +
    cowplot::draw_label(selected_cluster, x = 0, hjust = 0)

  label <- cowplot::ggdraw() +
    cowplot::draw_label('UMAP PLOTS', x=1, hjust = 1,  colour = 'dimgray', size = 12)

  title <- cowplot::plot_grid(title, label)

  subtitle <- cowplot::ggdraw() +
    cowplot::draw_label(paste0(ncells, ' cells - ', pcells, '%'), x = 0, hjust = 0, size = 12, colour = 'dimgray')

  # combine titles and plots
  marker_nrows <- ceiling(length(genes)/2)
  cowplot::plot_grid(title, subtitle, cluster_plot, gene_grid, ncol = 1, rel_heights = c(0.1, 0.1, 1, marker_nrows))
}


#' Remove ggplot xaxis title, text, and ticks
#'
#' @return \code{theme}
#' @export
#' @keywords internal
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
theme_no_yaxis <- function() {
  ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank())
}

#' Make ggplot axes and text dimgray
#'
#' @param with_nums Include axis ticks/text? Default is TRUE.
#'
#' @return \code{theme}
#' @export
#' @keywords internal
theme_dimgray <- function(with_nums = TRUE) {

  axis.line <- ggplot2::element_line(size = 0.1, color = 'dimgray')

  if (with_nums) {
    axis.text <- ggplot2::element_text(colour = 'dimgray')
    axis.ticks <- ggplot2::element_line(size = 0.1, color = 'dimgray')

  } else {
    axis.text <- axis.ticks <- ggplot2::element_blank()

  }

  ggplot2::theme(axis.line.y = axis.line,
                 axis.line.x = axis.line,
                 axis.ticks.x = axis.ticks,
                 axis.ticks.y = axis.ticks,
                 axis.text = axis.text,
                 axis.title = ggplot2::element_text(colour = 'dimgray'),
                 text = ggplot2::element_text(colour = 'dimgray'))
}
