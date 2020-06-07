#' Plot TSNE coloured by cluster
#'
#'
#' @return \code{ggplot}
#' @export
plot_tsne_cluster <- function(scseq, legend = FALSE, cols = NULL, title = NULL, ...) {

  levs <- levels(scseq$cluster)
  if (is.null(cols)) cols <- get_palette(levs)

  # dynamic label/point size
  pt.size <- min(6000/ncol(scseq), 2)
  nc <- length(levs)
  # shorten labels
  if (nc > 30)
    levels(scseq$cluster) <- shorten_cluster_labels(levs)

  num <- suppressWarnings(!anyNA(as.numeric(levs)))
  label.size <- if(num) 6 else if (nc>30) 5.2 else if(nc > 17) 5.5 else 6

  pl <- DimPlot(scseq, reduction = 'TSNE', cols = cols, pt.size = pt.size, label.size = label.size, ...) +
    theme_no_axis_vals() +
    ggplot2::xlab('TSNE1') +
    ggplot2::ylab('TSNE2') +
    theme_dimgray(with_nums = FALSE)


  if (!legend)
    pl <- pl + ggplot2::theme(legend.position = 'none', text = ggplot2::element_text(color = 'dimgray'))

  if (!is.null(title))
    pl <- pl + ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title.position = "plot",
                   plot.title = ggplot2::element_text(color = '#333333', hjust = 0, size = 16, face = 'plain', margin = ggplot2::margin(b = 25)))

  return(pl)
}

#' Shorted cluster labels when there are lots of clusters
#'
#' @param labels Character vector of labels
#'
#' @return \code{labels} where each word is shortened. trailing underscore numbers are preserved.
#' @export
#' @keywords internal
#' @examples
#' labels <- c('Non-Classic Mono', 'B-cell', 'B-cell_1', 'B-cell_12')
#' shorten_cluster_labels(labels)
#'
shorten_cluster_labels <- function(labels) {

  # keep unique numbers
  has.nums <- grepl('_\\d+$', labels)
  nums <- gsub('^.+?(_\\d+)$', '\\1', labels)
  new <- gsub('_\\d+$', '', labels)

  # keep at most 15 characters
  new <- stringr::str_trunc(new, 15)

  # restore underscore_nums
  new[has.nums] <- paste0(new[has.nums], nums[has.nums])

  return(new)
}

#' Plot TSNE coloured by gene or QC metric
#'
#' @param scseq \code{SingleCellExperiment} object.
#' @param feature Character vector specifying feature to colour cells by.
#' @param pt.size Numeric scalar, specifying the size of the points. Defaults to 3.
#'
#' @return \code{ggplot}
#' @export
plot_tsne_feature <- function(scseq, feature, cols = c('lightgray', 'blue'), reverse_scale = feature %in% c('ribo_percent', 'log10_sum', 'log10_detected'), title = paste0('Expression by Cell: ', feature)) {

  pt.size <- min(6000/ncol(scseq), 2)

  if (reverse_scale) cols <- rev(cols)

  FeaturePlot(scseq, feature, reduction = 'TSNE', pt.size = pt.size, order = TRUE, cols = cols) +
    theme_no_axis_vals() +
    ggplot2::xlab('TSNE1') +
    ggplot2::ylab('TSNE2') +
    ggplot2::theme(plot.title = ggplot2::element_blank()) +
    theme_dimgray(with_nums = FALSE) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour = 'black'),
                   plot.title.position = "plot",
                   plot.title = ggplot2::element_text(color = '#333333', hjust = 0, size = 16, face = 'plain', margin = ggplot2::margin(b = 25)))


}


#' Generate plotly for cluster origins of integrated dataset
#'
#'
#' @param scseq \code{SingleCellExperiment}
#' @param clust A cluster in \code{scseq} to generate plotly for
#' @param sc_dir Directory containing original datasets used when integrating \code{scseq}
#'
#' @return plotly
#' @export
#' @keywords internal
#'
plot_cluster_labels <- function(scseq, clust, sc_dir) {

  df <- tibble::as_tibble(scseq@colData) %>%
    dplyr::filter(cluster == clust) %>%
    dplyr::group_by(batch) %>%
    dplyr::mutate(nsample = n()) %>%
    dplyr::group_by(batch, orig.cluster) %>%
    dplyr::summarize(ncell = n(),
                     group = orig.ident[1],
                     nsample = nsample[1],
                     pcell = n() / nsample[1] * 100,
                     boc = paste(batch[1], orig.cluster[1], sep='_')) %>%
    dplyr::arrange(group, batch, pcell) %>%
    dplyr::select(-orig.cluster)

  # get mapping between original cluster index and name
  annots <- c()
  batches <- unique(df$batch)

  for (i in seq_along(batches)) {
    batch <- batches[i]
    annot_path <- file.path(sc_dir, batch, 'annot.rds')
    annot <- readRDS(annot_path)
    names(annot) <- paste(batch, seq_along(annot), sep='_')
    annots <- c(annots, annot)
  }

  df$oc <- annots[df$boc]
  df$customdata <- paste(df$ncell, 'of', df$nsample, 'cells')

  (pl <- plotly::plot_ly(data = df,
                         y = ~boc,
                         x = ~pcell,
                         color = ~batch,
                         customdata = ~customdata,
                         height = (nrow(df)*30) + 140,
                         text = ~nsample,
                         type = 'scatter',
                         mode = 'markers',
                         marker = list(size = 9, line = list(width = 1, color = '#333333')),
                         hoverlabel = list(bgcolor = '#000000', align = 'left'),
                         hovertemplate = paste0(
                           '<span style="color: crimson; font-weight: bold; text-align: left;">From Sample</span>: %{customdata}<br>',
                           '<extra></extra>')
  ) %>%
      plotly::layout(
        hovermode= 'closest',
        margin = list(t = 65, r = 20, l = 0, pad = 10),
        title = list(text = 'For Each Sample: Percent of Cells From Original Clusters', x = 0, y = .99, yanchor = 'top', font = list(size = 16, color = '#333333')),
        legend = list(traceorder = 'reversed', itemclick = FALSE, font = list(size = 14, color = '#333333')),
        xaxis = list(title = '', range = c(-4, 104), fixedrange=TRUE, side = 'top', zeroline = FALSE, tickfont = list(size = 14, color = '#333333')),
        yaxis = list(title = '', fixedrange=TRUE, tickvals = df$boc, ticktext = df$oc, tickmode = 2, gridwidth = 1, gridcolor = 'gray', ticklen = 10, tickcolor = 'white', tickfont = list(size = 14, color = '#333333'))
      ) %>%
      plotly::config(displayModeBar = 'hover',
                     displaylogo = FALSE,
                     doubleClick = 0,
                     showAxisDragHandles = FALSE,
                     showAxisRangeEntryBoxes = FALSE,
                     showTips = FALSE,
                     modeBarButtonsToRemove = c('zoom2d',
                                                'pan2d',
                                                'autoScale2d',
                                                'resetScale2d',
                                                'hoverClosestCartesian',
                                                'hoverCompareCartesian',
                                                'select2d',
                                                'lasso2d',
                                                'zoomIn2d',
                                                'zoomOut2d',
                                                'toggleSpikelines'),
                     toImageButtonOptions = list(format = "png")
      ))

  return(pl)
}

#' Format gene plots for sample comparison for drugseqr app
#'
#' @param group Level in \code{scseq$orig.ident} to show cells for. Either \code{'ctrl'} or \code{'test'}
#' @param scseq \code{SingleCellExperiment} object.
#'
#' @return \code{plot} formatted for drugseqr app
#' @export
#' @keywords internal
plot_tsne_feature_sample <- function(gene, scseq, group = 'test') {

  plot <- plot_tsne_feature(scseq, gene)
  gene <- make.names(gene)

  # the min and max gene expression value
  lims <- range(plot$data[[gene]])

  # show selected group only
  sel.cells <- colnames(scseq)[scseq$orig.ident == group]
  plot$data <- plot$data[row.names(plot$data) %in% sel.cells, ]

  # add selected group as title
  plot <- plot + ggplot2::ggtitle(toupper(group)) +
    ggplot2::theme(plot.title = ggplot2::element_text(color = 'black'))

  # use the same scale for each sample so that comparable
  suppressMessages(plot <- plot + ggplot2::scale_color_continuous(low ="lightgray", high = "blue", limits = lims))

  # remove control plot labels and legend
  if (group == 'ctrl')
    plot <- plot + ggplot2::xlab('') + ggplot2::ylab('') +
    ggplot2::theme(legend.position = 'none')

  return(plot)
}


#' Plot single-cell ridgeline plots
#'
#' @param feature Feature name to generate ridge plot for. Either a row or \code{colData} of \code{scseq}.
#' @param scseq \code{SingleCellExperiment}.
#' @param selected_cluster Name of the selected cluster.
#' @param by.sample if \code{TRUE} plot \code{feature} ridge for each \code{scseq$batch}. Default (\code{FALSE})
#'  will plot \code{feature} for each \code{scseq$cluster}.
#' @param with.height Whether to return height of plot. See value for details.
#' @param decreasing if \code{TRUE}, ridgelines with smaller mean values of \code{feature} will show up on top.
#'  Used to show features where smaller values indicate potential QC issues.
#'
#' @return ggplot object
#' @export
#' @keywords internal
plot_ridge <- function(feature, scseq, selected_cluster, by.sample = FALSE, with.height = FALSE, decreasing = feature %in% c('ribo_percent', 'log10_sum', 'log10_detected')) {
  if (is.null(selected_cluster)) selected_cluster <- ''

  # for one vs one comparisons
  selected_cluster <- strsplit(selected_cluster, '-vs-')[[1]]

  # get color for selected cluster
  clus_levs <- levels(scseq$cluster)
  seli <- as.numeric(selected_cluster)
  sel <- clus_levs[seli]
  nsel <- length(sel)
  colors <- get_palette(clus_levs)
  color  <- colors[seli]

  # either highlight test group or selected cluster
  if (by.sample) {
    if (!isTruthy(selected_cluster)) scseq <- scseq[row.names(scseq) != feature, ]
    scseq <- scseq[, scseq$cluster %in% sel]
    y <- factor(scseq$batch)
    hl <- scseq$orig.ident
    title <- paste('Expression by Sample:', feature, 'in', sel)

  } else {
    y <- scseq$cluster

    is.char <- suppressWarnings(is.na(as.numeric(clus_levs)))
    if (any(is.char)) {
      clus_levs <- paste0(seq_along(clus_levs), ': ', clus_levs)
      levels(y) <- clus_levs
      sel <- clus_levs[seli]
    }

    hl <- as.character(y)
    hl[!hl %in% sel] <- 'out'
    hl <- factor(hl, levels = c(sel, 'out'))
    title <- paste('Expression by Cluster:', feature)
  }


  if (feature %in% row.names(scseq))
    x <- as.numeric(SingleCellExperiment::logcounts(scseq[feature, ]))
  else
    x <- scseq[[feature]]

  # errors if boolean/factor/NULL
  if (is.null(x) || !is.numeric(x)) {
    pl <- NULL
    if (with.height) pl <- list(plot = plot, height = 453)
    return(pl)
  }

  m <- tapply(x, y, mean)
  y <- factor(y, levels = levels(y)[order(m, decreasing = decreasing)])
  df <- data.frame(x, hl, y) %>%
    dplyr::add_count(y) %>%
    dplyr::filter(n > 2)


  (pl <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = hl, alpha = hl, color = hl, height = ..density..)) +
      ggplot2::scale_fill_manual(values = c(color, 'gray')) +
      ggplot2::scale_alpha_manual(values = c(rep(0.95, nsel), 0.25)) +
      ggplot2::scale_color_manual(values = c(rep('black', nsel), 'gray')) +
      ggridges::geom_density_ridges(scale = 3, rel_min_height = 0.001, stat = "density", trim = TRUE) +
      ggridges::theme_ridges(center_axis_labels = TRUE) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::coord_cartesian(clip = "off") +
      theme_dimgray() +
      ggplot2::xlab('') +
      ggplot2::ggtitle(title) +
      ggplot2::theme(legend.position = 'none',
                     plot.title.position = "plot", panel.grid.major.x = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(color = '#333333', size = 16, face = 'plain', margin = ggplot2::margin(b = 25) ),
                     axis.text.y = ggplot2::element_text(color = '#333333', size = 14),
                     axis.text.x = ggplot2::element_text(color = '#333333', size = 14),
                     axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank()
      ))


  ncells <- tapply(df$x, as.character(df$y), length)
  ncells <- ncells[levels(df$y)]
  xlim <- ggplot2::ggplot_build(pl)$layout$panel_scales_x[[1]]$range$range[2]



  for (i in seq_along(ncells)) {
    pl <- pl + ggplot2::annotation_custom(
      grob = grid::textGrob(label = paste(ncells[i], 'cells'), vjust = -0.3, just = 'right',
                            gp = grid::gpar(fontsize = 14, col = '#333333')),
      ymax = i,
      ymin = i,
      xmax = xlim,
      xmin = xlim
    )
  }

  if (with.height) pl <- list(plot = pl, height = max(453, length(ncells) * 50))
  return(pl)
}


#' Plot dprimes for integrated single-cell dataset with replicates.
#'
#' @param gene a row name within data.frames of \code{tts} to plot dprimes for.
#' @param annot Charactor vector with cluster names for \code{tts}.
#' @param selected_cluster Name of \code{tts} indicating selected cluster.
#' @param tts List of topTables, one for each cluster.
#' @param exclude_ambient If true, ambient dprimes are greyed out along with non-significant genes.
#'
#' @return plotly object.
#' @export
#' @keywords internal
plot_scseq_dprimes <- function(gene, annot, selected_cluster, tts, exclude_ambient) {

  tt <- lapply(tts, function(x) x[gene, ])
  tt <- do.call(rbind, tt)
  tt <- tt[!is.na(tt$t), ]
  if (nrow(tt) == 0) return(NULL)
  path_df <- get_path_df(tt, path_id = '')

  seli <- which(row.names(tt) == selected_cluster)
  clusters <- as.numeric(row.names(tt))

  # highlight clusters where this gene is significant
  row.names(tt) <- row.names(path_df) <- path_df$Gene <- annot[clusters]

  path_df$ambient <- tt$ambient
  link <- as.character(path_df$Gene)
  path_df$Link <- paste0('<span style="color: dimgray">', link, '</span>')
  path_df$Link[seli] <- gsub('dimgray', 'black', path_df$Link[seli])
  path_df$color <- ifelse(tt$adj.P.Val < 0.05, 'black', 'gray')

  if (exclude_ambient) path_df$color[path_df$ambient] <- 'gray'

  # plot trips up if numbered clusters
  is.number <-  suppressWarnings(!is.na(as.numeric(link)))
  path_df$Gene[is.number] <- paste('Cluster', link[is.number])

  path_df <- path_df[order(abs(tt$dprime), decreasing = TRUE), ]
  plot_dprimes(path_df, drugs = FALSE)
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
                   plot.title.position = "plot",
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_line(size = 0),
                   panel.grid.major.y = ggplot2::element_line(linetype = 'longdash', size = 0.5),
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
  if (nlevs == 2) {
    values <- c("#729ECE", "#FF9E4A")

  } else if (nlevs <= 10) {
    values <- head(tableau10medium, nlevs)

  } else if (nlevs <= 20) {
    values <- head(tableau20, nlevs)

  } else {
    set.seed(0)
    values <- randomcoloR::distinctColorPalette(nlevs)

  }
  return(values)
}
