#' Plot TSNE coloured by cluster
#'
#'
#' @return \code{ggplot}
#' @export
plot_tsne_cluster <- function(scseq, selected_clusters = levels(scseq$cluster), legend_title = 'Cluster', cols = NULL) {

  colnames(scseq) <- make.unique(colnames(scseq))
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
#' @param pt.size Numeric scalar, specifying the size of the points. Defaults to 3.
#'
#' @return \code{ggplot}
#' @export
plot_tsne_gene <- function(scseq, gene) {

  colnames(scseq) <- make.unique(colnames(scseq))
  scseq <- Seurat::as.Seurat(scseq, counts = NULL)
  Idents(scseq) <- scseq$cluster

  pt.size <- min(6000/ncol(scseq), 2)

  gene_plot <- Seurat::FeaturePlot(scseq, gene, reduction = 'TSNE', pt.size = pt.size, order = TRUE) +
    theme_no_axis_vals() +
    ggplot2::xlab('TSNE1') +
    ggplot2::ylab('TSNE2') +
    ggplot2::theme(plot.title = ggplot2::element_blank()) +
    theme_dimgray(with_nums = FALSE) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour = 'black'))

  return(gene_plot)
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
    ggplot2::geom_point(stat = "identity", size = 4, colour="#333333", pch=21) +
    ggplot2::scale_fill_distiller(palette = 'Reds', name = '', direction = 1) +
    ggplot2::xlab('') +
    ggplot2::ylab('') +
    ggplot2::ggtitle('BioGPS Human Gene Atlas Expression') +
    ggplot2::coord_flip() +
    ggpubr::theme_pubr()  +
    theme_dimgray() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 16, color = '#333333', margin = ggplot2::margin(b = 25)),
                   axis.text.y = ggplot2::element_text(color = '#333333'),
                   axis.text = ggplot2::element_text(size = 14),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_line(size = 0),
                   panel.grid.major.y = ggplot2::element_line(linetype = 'longdash', size = 0.1),
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


