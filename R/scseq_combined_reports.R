#' Save PDF of combined scRNA-seq reports
#'
#' For integrated data with multiple samples.
#'
#' @param scseq \code{Seurat} object with \code{scseq$orig_clusters} providing the originally identified clusters for
#'  each sample. This permits visualizing how integration shift the TSNE coordinates and cluster assignments.
#' @inheritParams save_scseq_reports
#'
#' @return Saves report to \code{fname}.
#' @export
#'
#' @examples
save_combined_scseq_reports <- function(scseq, markers, orig_clusters, fname, point_size = 3) {


  if (class(scseq) == 'Seurat') scseq <- srt_to_sce(scseq, 'SCT')
  group_idents <- scseq$orig.ident
  groups <- unique(group_idents)


  # get scale that uses all possible levels
  all_levels <- c(names(markers), levels(unlist(orig_clusters)))
  all_levels <- unique(all_levels)
  fill_scale <- scale_fill_names(all_levels)


  # new TSNE coords with original clusters (first page)
  orig_plots <- list()
  for (i in seq_along(groups)) {
    group <- groups[i]
    group_scseq <- scseq[, group_idents == group]
    group_scseq$cluster <- factor(orig_clusters[[group]], levels = all_levels)

    suppressMessages(orig_plots[[i]] <- plot_tsne_cluster(group_scseq, legend_title = toupper(group)) +
                       theme_dimgray() + theme_no_xaxis() + theme_no_yaxis() + fill_scale)
  }


  # put cluster levels in order of markers
  if (!all(names(markers) %in% levels(scseq$cluster)))
    stop ('markers are missing for some clusters')
  scseq$cluster <- factor(scseq$cluster, levels = names(markers))

  # plot if percentage cells
  pcells_barplot <- plot_scseq_barplot(scseq)

  # new TSNE coords with new clusters (first page)
  suppressMessages(combined_plot <- plot_tsne_cluster(scseq, legend_title = 'COMBINED') +
                     theme_dimgray() + theme_no_xaxis() + theme_no_yaxis() + fill_scale)

  # get titles
  orig_title <- cowplot::ggdraw() +
    cowplot::draw_label('Original Clusters', x = 0, hjust = 0)

  label <- cowplot::ggdraw() +
    cowplot::draw_label('TSNE PLOTS', x=1, hjust = 1,  colour = 'dimgray', size = 12)

  orig_title <- cowplot::plot_grid(orig_title, label, ncol=2)

  combined_title <- cowplot::ggdraw() +
    cowplot::draw_label('Combined Clusters', x = 0, hjust = 0)

  combined_title <- cowplot::plot_grid(combined_title,  label, ncol=2)

  # combine plots with original clusters
  orig_grid <- cowplot::plot_grid(orig_title, orig_plots[[1]], orig_plots[[2]],
                                  ncol=1, rel_heights = c(.2,1,1))

  combined_grid <- cowplot::plot_grid(combined_title, combined_plot, pcells_barplot,
                                      ncol=1, rel_heights = c(.2,1,1.5), align = 'v', axis = 'l')

  # create pdf
  pdf(file = fname, paper = 'US', width = 8.50, height = 11.0, title = 'combined cluster markers')
  plot(orig_grid)
  plot(combined_grid)

  for (i in seq_along(markers)) {
    reports <- plot_combined_scseq_report(scseq, markers = markers[i], point_size = point_size)
    for (report in reports) plot(report)
  }
  dev.off()
}

#' Plot scRNA-seq barplot
#'
#' Plots percentage of cells within each \code{scseq$cluster} split by \code{scseq$orig.ident}.
#'
#' @param scseq \code{SingleCellExperiment}
#'
#' @return \code{ggplot}
#' @export
#' @keywords internal
#'
#' @examples
plot_scseq_barplot <- function(scseq) {

  group_idents <- scseq$orig.ident
  groups <- unique(group_idents)

  # percent of cells per group in new cluster
  combined_pcells <- list()
  for (i in seq_along(groups)) {
    group <- groups[i]
    group_scseq <- scseq[, group_idents == group]
    combined_pcells[[i]] <- data.frame(table(group_scseq$cluster) / ncol(group_scseq) * 100,
                                       Group = toupper(group))
  }
  combined_pcells <- do.call(rbind, combined_pcells)

  pcells_plot <- ggplot2::ggplot(combined_pcells, ggplot2::aes(x=Var1, y=Freq, fill=Group)) +
    ggplot2::geom_bar(stat='identity',position="dodge") +
    cowplot::theme_cowplot() +
    scale_fill_names(toupper(groups)) +
    ggplot2::xlab('') +
    ggplot2::guides(fill=ggplot2::guide_legend(title="Percent Cells")) +
    theme_dimgray() +
    ggplot2::scale_y_continuous(expand = c(0, 0), position = 'right', sec.axis = ggplot2::dup_axis()) +
    ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   axis.title.y = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_line(),
                   axis.text.y.left = ggplot2::element_blank(),
                   axis.ticks.y.left = ggplot2::element_blank(),
                   plot.margin=ggplot2::unit(c(20, 5.5, 5.5, 5.5), "points"))

  return(pcells_plot)
}


#' Plot grid of combined scRNA-seq clusters and gene markers together
#'
#' For integrated data with multiple samples.
#'
#' @inheritParams plot_scseq_report
#'
#' @return List of \code{ggplot} objects
#' @export
#' @keywords internal
#'
#' @examples
plot_combined_scseq_report <- function(scseq, markers, point_size = 3) {
  if (class(scseq) == 'Seurat') scseq <- srt_to_sce(scseq, 'SCT')
  selected_group <- names(markers)
  genes <- markers[[1]]

  # make orig.ident the clusters so that can highlight
  group_idents <- scseq$orig.ident
  groups <- unique(group_idents)
  scseq$orig.ident <- scseq$cluster

  # get grid of gene plots by group
  gene_plots <- list()
  for (i in 1:length(genes)) {
    for (j in 1:length(groups)) {
      idx <- length(gene_plots) + 1

      gene_plots[[idx]] <-
        plot_tsne_gene(scseq, genes[i], hide_mask = group_idents != groups[j], point_size = point_size) +
        theme_dimgray() +
        theme_no_xaxis() +
        theme_no_yaxis() +
        ggplot2::theme(legend.position = 'none', plot.title=ggplot2::element_text(size=12, hjust = 0)) +
        ggplot2::ggtitle(genes[i])
    }
  }

  # show clusters by group ----
  cluster_plots <- list()
  for (j in seq_along(groups)) {
    group <- groups[[j]]
    group_scseq <- scseq[, group_idents == group]

    ncells <- sum(group_scseq$orig.ident == selected_group)
    pcells <- round(ncells / ncol(group_scseq) * 100)
    subtitle <- paste0(ncells, ' cells - ', pcells, '%')

    cluster_plots[[j]] <- plot_tsne_cluster(group_scseq, selected_groups = selected_group, point_size = point_size) +
      ggplot2::theme(plot.margin = ggplot2::unit(c(20, 5.5, 20, 5.5), "points")) +
      theme_dimgray() +
      theme_no_xaxis() +
      theme_no_yaxis() +
      ggplot2::theme(legend.position = 'none',
                     plot.title = ggplot2::element_text(size=14, color = 'black', hjust = 0)) +
      ggplot2::ggtitle(toupper(group), subtitle)
  }

  cluster_grid <- cowplot::plot_grid(plotlist = cluster_plots, align='vh', ncol = 2)


  # max two genes per page
  gene_grids <- list()
  gene_nrows <- c()
  ngplots <- length(gene_plots)
  for (i in seq(1, ngplots, 4)) {
    idx <- length(gene_grids) + 1
    plotlist <- head(gene_plots[i:ngplots], 4)
    gene_grids[[idx]] <- cowplot::plot_grid(plotlist = plotlist, ncol = 2, align = 'vh')
    gene_nrows <- c(gene_nrows, ceiling(length(plotlist) / 2))
  }

  # title
  title <- cowplot::ggdraw() +
    cowplot::draw_label(selected_group, x = 0, hjust = 0)

  label <- cowplot::ggdraw() +
    cowplot::draw_label('TSNE PLOTS', x=1, hjust = 1,  colour = 'dimgray', size = 12)

  title <- cowplot::plot_grid(title, label)

  # plot each page with a title, cluster plots, and up to two rows of gene plots
  pages <- list()
  for (i in seq_along(gene_grids)) {
    pages[[i]] <- cowplot::plot_grid(title, cluster_grid, gene_grids[[i]],
                                     align = 'vh', ncol = 1, rel_heights = c(0.1, 1, gene_nrows[i]))
  }
  return(pages)
}

