#' Plot TSNE coloured by cluster
#'
#' @param scseq SingleCellExperiment
#' @param legend Should legend be kept? Default is \code{FALSE}.
#' @param title Title of plot
#' @param ... additional argument to \code{DimPlot}
#' @inheritParams DimPlot
#'
#'
#' @return \code{ggplot}
#' @keywords internal
plot_cluster <- function(scseq = NULL, reduction = 'TSNE', plot_data = NULL, legend = FALSE, cols = NULL, title = NULL, title_color ='#333333', title_hjust = 0, ...) {


  # dynamic label/point size
  if (!is.null(scseq)) {
    levs <- levels(scseq$cluster)
    pt.size <- min(6000/ncol(scseq), 2)
    nc <- length(levs)
    if (nc > 30) levels(scseq$cluster) <- shorten_cluster_labels(levs)

    reds <- SingleCellExperiment::reducedDimNames(scseq)
    reds <- reds[reds %in% c('UMAP', 'TSNE')]
    if (!reduction %in% reds) reduction <- reds[1]

  } else {
    levs <- levels(plot_data$cluster)
    pt.size <- min(6000/nrow(plot_data), 2)
    nc <- length(levs)
    if (nc > 30) {
      levels(plot_data$cluster) <- shorten_cluster_labels(levs)
      add_warn <- TRUE
    }
    reduction <- gsub('1$', '', colnames(plot_data)[1])
  }

  if (is.null(cols))
    cols <- get_palette(levs, with_all = TRUE)

  num <- suppressWarnings(!anyNA(as.numeric(levs)))
  label.size <- if(num) 6 else if (nc>30) 5.2 else if(nc > 17) 5.5 else 6

  pl <- DimPlot(scseq, plot_data, reduction = reduction, cols = cols, pt.size = pt.size, label.size = label.size, ...) +
    theme_no_axis_vals() +
    ggplot2::xlab(paste0(reduction, '1')) +
    ggplot2::ylab(paste0(reduction, '2')) +
    theme_dimgray(with_nums = FALSE)



  if (!legend)
    pl <- pl + ggplot2::theme(legend.position = 'none', text = ggplot2::element_text(color = 'dimgray'))

  if (nc > 30) {
    title <- 'Labels Omitted'
    title_color <- 'lightgray'
    title_hjust <- 1
  }

  if (!is.null(title))
    pl <- pl + ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title.position = "plot",
                   plot.title = ggplot2::element_text(color = title_color, hjust = title_hjust, size = 16, face = 'plain', margin = ggplot2::margin(b = 15)))

  return(pl)
}

update_cluster_plot <- function(plot, annot, hl) {
  annot <- format_ridge_annot(annot)
  if (nrow(plot$layers[[2]]$data) != length(annot)) return(NULL)

  plot$layers[[2]]$data$cluster <- annot
  cols <- rep('black', length(annot))

  if (!is.null(hl)) {
    cols <- rep('#5f5f5f', length(annot))
    cols[hl] <- 'black'
  }
  plot$layers[[2]]$aes_params$colour <- cols

  return(plot)
}

#' Shorted cluster labels when there are lots of clusters
#'
#' @param labels Character vector of labels
#'
#' @return \code{labels} where each word is shortened. trailing underscore numbers are preserved.
#'
#' @keywords internal
#' @examples
#' labels <- c('Non-Classic Mono', 'B-cell', 'B-cell_1', 'B-cell_12')
#' dseqr:::shorten_cluster_labels(labels)
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
#' @param reverse_scale Should the scale be reversed (lower values of feature
#'   are darker)?
#' @param title Title to add to plot.
#' @inheritParams FeaturePlot
#'
#'
#' @return \code{ggplot}
#' @keywords internal
plot_feature <- function(scseq,
                         feature,
                         reduction = 'TSNE',
                         cols = c('lightgray', 'blue'),
                         reverse_scale = feature %in% c('ribo_percent', 'log10_sum', 'log10_detected'),
                         title = paste0('Expression by Cell: ', feature),
                         pt.size = min(6000/ncol(scseq), 2)) {


  if (reverse_scale) cols <- rev(cols)

  reds <- SingleCellExperiment::reducedDimNames(scseq)
  reds <- reds[reds %in% c('UMAP', 'TSNE')]
  if (!reduction %in% reds) reduction <- reds[1]

  FeaturePlot(scseq, feature, reduction = reduction, pt.size = pt.size, order = TRUE, cols = cols) +
    theme_no_axis_vals() +
    ggplot2::xlab(paste0(reduction, '1')) +
    ggplot2::ylab(paste0(reduction, '2')) +
    ggplot2::theme(plot.title = ggplot2::element_blank()) +
    theme_dimgray(with_nums = FALSE) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour = 'black'),
                   plot.title.position = "plot",
                   plot.title = ggplot2::element_text(color = '#333333', hjust = 0, size = 16, face = 'plain', margin = ggplot2::margin(b = 25)))


}


#' Format gene plots for sample comparison for dseqr app
#'
#' @param group Level in \code{scseq$orig.ident} to show cells for. Either \code{'ctrl'} or \code{'test'}
#' @param scseq \code{SingleCellExperiment} object.
#'
#' @return \code{plot} formatted for dseqr app
#'
#' @keywords internal
plot_feature_sample <- function(feature, scseq, group, reduction = 'TSNE', plot = NULL, cols = c('lightgray', 'blue')) {

  if (is.null(plot)) plot <- plot_feature(scseq, feature, reduction)

  is.group <- scseq$orig.ident == group
  ncells <- sum(is.group)
  col <- make.names(feature)
  plot$layers[[1]]$aes_params$size <-  min(6000/(ncells), 2)

  # the min and max feature expression value
  lims <- range(plot$data[[col]])

  # show selected group only
  sel.cells <- colnames(scseq)[is.group]
  plot$data <- plot$data[row.names(plot$data) %in% sel.cells, ]

  # add selected group as title
  fname <- ifelse(group == 'test', 'Test', 'Control')
  plot <- plot + ggplot2::ggtitle(paste("Expression in", fname, 'Samples:', feature))

  # use the same scale for each sample so that comparable
  suppressMessages(plot <- plot + ggplot2::scale_color_continuous(low = cols[1], high = cols[2], limits = lims))

  return(plot)
}

#' Update Feature Plot
#'
#' Used to change feature with fast re-plot
#'
#' @param plot ggplot2 object
#' @param feature_data Named numeric vector of feature values. Names are cell
#'   names and must be a subset of \code{row.names(plot$data)}
#' @param feature Name of feature
#' @param reverse_scale Should low values be emphasized? Used for features where
#'  low values indicate possible quality issues.
#'
#' @return ggplot2 object
#' @keywords internal
#'
update_feature_plot <- function(plot,
                                feature_data,
                                feature,
                                reverse_scale = feature %in% c('ribo_percent', 'log10_sum', 'log10_detected')) {

  prev <- tail(colnames(plot$data), 1)
  col <- make.names(feature)
  plot$data <- plot$data[names(feature_data), ]
  plot$data[[prev]] <- NULL
  plot$data[[col]] <- feature_data[row.names(plot$data)]

  # order by expression
  ord <- order(plot$data[[col]], decreasing = reverse_scale)
  plot$data <- plot$data[ord, ]

  # update mapping to new gene
  plot$layers[[1]]$mapping$colour <- rlang::quo_set_expr(
    plot$layers[[1]]$mapping$colour, rlang::ensym(col))

  if (reverse_scale) suppressMessages(
    plot <- plot +
      ggplot2::scale_color_gradientn(
        colors = c('blue', 'lightgray'),
        guide = "colorbar"
      ))

  plot <- plot + ggplot2::ggtitle(paste('Expression by Cell:', feature))
  return(plot)

}

#' Downsample cells within each cluster across a group
#'
#' Used for sample comparison plots to show equal number of test and control
#' cells for each cluster
#'
#' @param scseq SingleCellExperiment
#'
#' @return \code{scseq} with equal number
#' @keywords internal
#'
downsample_group <- function(scseq, group = 'orig.ident') {
  ncells <- tibble::tibble(cluster=scseq$cluster, group=scseq[[group]], cell = colnames(scseq))

  keep <- ncells %>%
    dplyr::group_by(group, cluster) %>%
    dplyr::add_count() %>%
    dplyr::ungroup() %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(nmin = min(n)) %>%
    dplyr::group_by(group, cluster) %>%
    dplyr::group_modify(~dplyr::sample_n(.x, nmin)) %>%
    dplyr::pull(cell)

  return(scseq[, keep])

}

#' Downsample SingleCellExperiment clusters
#'
#' Gets a maximum number of cells in each cluster. Used for generate mini
#' datasets for faster label transfer.
#'
#' @param scseq SingleCellExperiment
#' @param max.cells maximum number of cells in each \code{scseq$cluster}
#'
#' @return scseq
#' @keywords internal
#'
downsample_clusters <- function(scseq, max.cells = 200) {
  cells <- data.frame(id = colnames(scseq),
                      cluster = scseq$cluster)

  set.seed(0)
  keep <- cells %>%
    dplyr::group_by(cluster) %>%
    dplyr::sample_n(min(max.cells, dplyr::n())) %>%
    dplyr::pull(id)

  scseq[, keep]
}


#' Get data for single-cell ridgeline plots
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
#' @return list used by \link{VlnPlot}
#'
#' @keywords internal
get_ridge_data <- function(feature, scseq, selected_cluster, by.sample = FALSE, decreasing = feature %in% c('ribo_percent', 'log10_sum', 'log10_detected'), with_all = FALSE, dgrlogs = NULL) {
  n <- NULL

  if (isTruthy(selected_cluster)) {
    # for one vs one comparisons
    selected_cluster <- strsplit(selected_cluster, '-vs-')[[1]]

    # get color for selected cluster
    seli <- as.numeric(selected_cluster)

  } else {
    selected_cluster <- ''
    seli <- 0
  }

  clus_levs <- levels(scseq$cluster)

  # specific args if selecting 'All Clusters'
  maxi <- length(clus_levs)+1
  if (seli[1] > maxi) {
    # invalid 'flip' state
    return(NULL)

  } else if (seli[1] == maxi) {
    sel <- c(clus_levs, 'All Clusters')
    nsel <- 1

  } else {
    sel <- clus_levs[seli]
    nsel <- length(sel)
  }

  colors <- get_palette(clus_levs, with_all=with_all)
  color  <- colors[seli]
  colors_dark <- get_palette(clus_levs, dark=TRUE, with_all=with_all)
  color_dark  <- colors_dark[seli]


  # either highlight test group or selected cluster
  if (by.sample) {
    keep <- scseq$cluster %in% sel
    y <- factor(scseq$batch[keep])
    hl <- scseq$orig.ident[keep]

    # cluster name get's appended by VlnPlot
    title <- paste('Expression by Sample:', feature, 'in')

  } else {
    hl <- as.integer(scseq$cluster)
    y <- factor(hl)
    hl[!hl %in% seli] <- 'out'
    hl <- factor(hl, levels = c(seli, 'out'))
    title <- paste('Expression by Cluster:', feature)
  }


  if (feature %in% row.names(scseq)) {
    if (!is.null(dgrlogs)) x <- fast_dgr_row(dgrlogs, feature)
    else x <- as.numeric(SingleCellExperiment::logcounts(scseq[feature, ]))
    if (exists('keep')) x <- x[keep]

  } else {
    x <- scseq[[feature]]
  }

  m <- tapply(x, y, mean)
  clus_ord <- order(m, decreasing = decreasing)
  y <- factor(y, levels = levels(y)[clus_ord])
  df <- data.frame(x, hl, y) %>%
    dplyr::add_count(y)

  ncells <- tapply(df$x, as.character(df$y), length)
  ncells <- ncells[levels(df$y)]

  # down sample to reduce plot size
  df <- df %>%
    dplyr::group_by(y) %>%
    dplyr::mutate(nsamp = min(n, 1000)) %>%
    dplyr::sample_n(nsamp)

  res <- list(
    df = df,
    nsel = nsel,
    seli = seli,
    color = color,
    color_dark = color_dark,
    ncells = ncells,
    clus_ord = clus_ord,
    clus_levs = clus_levs,
    title = title,
    by.sample = by.sample
  )

  return(res)
}

#' Add cluster numbers to annotation
#'
#' Adds cluster number to non-numeric cluster names.
#' E.g. if \code{annot = c('Monocytes', '2')} then
#' \code{annot = c('1: Monocytes', '2')} is returned
#'
#' @param annot Character vector of cluster names
#'
#' @return \code{annot} with cluster numbers pre-pended to non-numeric values.
#' @keywords internal
#'
format_ridge_annot <- function(annot) {

  is.char <- suppressWarnings(is.na(as.numeric(annot)))
  if (any(is.char)) {
    annot[is.char] <- paste0(seq_along(annot)[is.char], ': ', annot[is.char])
  }
  return(annot)
}

VlnPlot <- function(feature = NULL,
                    scseq = NULL,
                    selected_cluster = NULL,
                    by.sample = FALSE,
                    with_all = FALSE,
                    with.height = FALSE,
                    ridge_data = NULL,
                    is_mobile = FALSE,
                    decreasing = feature %in% c('ribo_percent',
                                                'log10_sum',
                                                'log10_detected')) {

  if (is.null(ridge_data)) {
    ridge_data <- get_ridge_data(
      feature, scseq, selected_cluster, by.sample, decreasing, with_all)
  }

  list2env(ridge_data, envir = environment())

  if (by.sample) {
    title <- paste(title, c(clus_levs, 'All Clusters')[seli])
    if (is_mobile) df <- shorten_y(df)
  } else {
    annot <- format_ridge_annot(clus_levs)
    levels(df$y) <- annot[clus_ord]
    if (seli[1]) levels(df$hl) <- c(annot[seli], 'out')
  }

  # errors if boolean/factor/NULL
  if (is.null(df$x) || !is.numeric(df$x)) {
    pl <- NULL
    if (with.height) pl <- list(plot = pl, height = 453)
    return(pl)
  }

  data <- as.data.frame(df[,'x'])

  pl <- SingleExIPlot(data,
                      idents = df$y,
                      hl = df$hl,
                      pt.size = 0.001,
                      title = title,
                      color = color,
                      color_dark = color_dark,
                      nsel = nsel,
                      ncells = ncells)


  if (with.height) pl <- list(plot = pl, height = max(235, length(ncells) * 45))
  return(pl)

}

shorten_y <- function(df) {
  short <- df %>%
    group_by(y) %>%
    dplyr::slice(1) %>%
    group_by(hl) %>%
    arrange(as.character(y)) %>%
    mutate(new = paste0(toupper(hl), seq_along(hl)))

  levels(df$y) <- short$new
  return(df)
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
#'
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
  row.names(tt) <- row.names(path_df) <- path_df$Gene <- c(annot, 'All Clusters')[clusters]

  path_df$ambient <- tt$ambient
  link <- as.character(path_df$Gene)
  is.num <-  suppressWarnings(!is.na(as.numeric(link)))
  add.num <- !is.num & link != 'All Clusters'

  path_df$Link <- paste0('<span style="color: dimgray">', link, '</span>')
  path_df$Link[seli] <- gsub('dimgray', 'black', path_df$Link[seli])
  path_df$Link[add.num] <- paste0(path_df$Link[add.num], ' [', seq_along(link)[add.num], ']')
  path_df$color <- ifelse(tt$adj.P.Val < 0.05, 'black', 'gray')

  if (exclude_ambient) path_df$color[path_df$ambient] <- 'gray'

  # plot trips up if numbered clusters
  path_df$Gene[is.num] <- paste('Cluster', link[is.num])

  path_df <- path_df[order(abs(tt$dprime), decreasing = TRUE), ]

  # move 'All Clusters' to top for single-cell dprimes plot
  clusts <- row.names(path_df)
  if ('All Clusters' %in% clusts) {
    path_df <- rbind(path_df['All Clusters', ],
                     path_df[setdiff(clusts, 'All Clusters'), ])

    path_df$logfc[1] <- 'n/a'
    path_df$ambient[1] <- 'n/a'
  }

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
#'
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
    theme_pubr()  +
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
#'
#' @keywords internal
get_palette <- function(levs, dark = FALSE, with_all = FALSE) {
  if (with_all) levs <- c(levs, 'All Clusters')

  # modified palettes from scater
  tableau20 <- c("#1F77B4", "#FF7F0E", "#FFD8B1", "#FFBB78", "#2CA02C", "#66CDAA",
                 "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
                 "#8C564B", "#C49C94", "#E377C2", "#F7B6D2",
                 "#17BECF", "#9EDAE5", "#CDCC5D", "#FFFACD", "#000000")

  tableau20_dark <- c("#004771", "#B69A7E", "#984802", "#A06705", "#036003", "#499279",
                      "#33870E", "#821919", "#CF1701", "#5B3979", "#7D5E91",
                      "#55342D", "#7B574F", "#A22283", "#CE308A",
                      "#056F79", "#028491", "#777600", "#B6B392", "#333333")

  tableau10medium <- c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                       "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                       "#CDCC5D", "#6DCCDA")

  tableau10medium_dark <- c("#365D83", "#9C5800", "#33702A", "#A22616",
                            "#6D4B86", "#644741", "#B62B8A", "#5E5E5E",
                            "#777600", "#097984")


  set.seed(0)
  nlevs <- length(levs)
  if (nlevs == 2) {
    values <- c("#729ECE", "#FF9E4A")

  } else if (nlevs <= 10) {
    pal <- if(dark) tableau10medium_dark else tableau10medium
    values <- sample(pal, nlevs)

  } else if (nlevs <= 20) {
    pal <- if(dark) tableau20_dark else tableau20
    values <- sample(pal, nlevs)

  } else {
    pal <- grDevices::colors(distinct = TRUE)
    pal <- pal[grep('gr(a|e)y|white|ivory|beige|seashell|snow',
                    pal,
                    invert = TRUE)]
    pal <- sample(pal, nlevs)
    values <- col2hex(pal, dark)
  }
  return(values)
}

col2hex <- function(cname, dark) {
  colMat <- grDevices::col2rgb(cname)
  if (dark) colMat <- colMat/1.4
  grDevices::rgb(red = colMat[1, ]/255,
                 green = colMat[2, ]/255,
                 blue = colMat[3, ]/255)
}


plot_scseq_diff <- function(pt.dat, feature = 'abundance') {

  title <- 'Abundance in Test vs Control Samples'
  delta_name <- 'Δ cells'

  if (feature != 'abundance') {
    title <- paste('Expression in Test vs Control Samples:', feature)
    delta_name <- 'Δ gene'
  }

  # color per direction
  fill <- c('#A50F15', '#08519C', '#FFFFFF')
  names(fill) <- levels(pt.dat$direction)
  fill <- fill[names(fill) %in% pt.dat$direction]

  n.alpha <- length(unique(pt.dat$alpha))
  b.alpha <- c(.5, 0.1)
  l.alpha <- c('< .05', '≥ .05')

  if (n.alpha == 1) {
    b.alpha <- 0.1
    l.alpha <- '≥ .05'
  }

  ggplot2::ggplot(ggplot2::aes(x=x, y=y, fill=direction, alpha=alpha), data=pt.dat) +
    ggplot2::geom_tile(color='lightgray') +
    cowplot::theme_cowplot() +
    theme_no_axis_vals() +
    theme_dimgray(with_nums = FALSE) +
    ggplot2::scale_fill_manual(name = delta_name, values = fill) +
    ggplot2::scale_alpha_identity(name = 'p-val', breaks = b.alpha, labels = l.alpha, guide = "legend") +
    ggplot2::theme(legend.position = 'right',
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.title = ggplot2::element_text(size=12),
                   legend.text = ggplot2::element_text(size=10),
                   plot.title.position = "plot",
                   plot.title = ggplot2::element_text(color = "#333333", hjust = 0, size = 16, face = 'plain', margin = ggplot2::margin(b = 15))) +
    ggplot2::guides(fill = ggplot2::guide_legend(order = 1),
                    alpha = ggplot2::guide_legend(order = 2)) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab('') +
    ggplot2::ylab('')
}


get_gene_diff <- function(gene, tts, grid) {

  tts <- lapply(tts, function(x) x[gene,, drop=FALSE])
  tts <- tts[!is.na(tts)]
  tt <- data.table::rbindlist(tts, fill = TRUE)
  tt <- as.data.frame(tt)
  row.names(tt) <- names(tts)
  tt <- tt[!is.na(tt$logFC), c('P.Value', 'logFC')]
  tt$cluster <- row.names(tt)

  # get points in grid with cells
  pt.dat <- grid %>%
    add_count(grid_xi, grid_yi) %>%
    left_join(tt) %>%
    distinct() %>%
    filter(n > 0) %>% # at least n cells
    dplyr::rename(pval = P.Value, direction = logFC) %>%
    dplyr::mutate(pval = replace(pval, is.na(pval), 1),
                  direction = case_when(is.na(direction) ~ '?',
                                        direction > 0 ~ '↑',
                                        direction < 0 ~ '↓'),
                  alpha = case_when(pval >= 0.05 ~ 0.1,
                                    pval <  0.05 ~ range02(-log10(pval)))) %>%
    dplyr::mutate(direction = factor(direction, levels = c('↑', '↓', '?'), ordered = TRUE))

  return(pt.dat)
}

get_abundance_diff <- function(scseq, group = scseq$orig.ident, sample = scseq$batch) {
  reds <- SingleCellExperiment::reducedDimNames(scseq)
  red <- reds[reds %in% c('UMAP', 'TSNE')]

  if (red == 'UMAP') nx=120; ny=60
  if (red == 'TSNE') nx=100; ny=50

  red.mat <- SingleCellExperiment::reducedDim(scseq, red)

  dat <- data.frame(x=red.mat[,1], y=red.mat[,2], group=group, sample=sample)

  # downsample within group so that each sample has same number of cells
  dsamp <- dat %>%
    dplyr::add_count(sample) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(nmin=min(n)) %>%
    dplyr::group_by(sample) %>%
    dplyr::sample_n(nmin) %>%
    dplyr::ungroup() %>%
    dplyr::select(-n, -nmin)

  # downsample so that each group has same number of cells
  nmin <- min(table(dsamp$group))
  dsamp <- dsamp %>%
    dplyr::group_by(group) %>%
    dplyr::sample_n(nmin) %>%
    dplyr::ungroup()

  # create grid on non-downsampled
  x <- seq(min(dat$x), max(dat$x), length.out = nx)
  y <- seq(min(dat$y), max(dat$y), length.out = ny)

  # in nx*ny grid: get xy bin for downsampled points
  points <- cbind(dsamp$x, dsamp$y)
  binxy <- data.frame(x = findInterval(points[,1], x),
                      y = findInterval(points[,2], y))

  binxy$x <- factor(binxy$x, levels = 1:nx)
  binxy$y <- factor(binxy$y, levels = 1:ny)

  # bin by group for color (logFC sign)
  binxy$group <- factor(dsamp$group)
  gtab <- table(binxy)
  gtab <- as.data.frame.table(gtab)

  # in nx*ny grid: get xy bin for all points
  points <- cbind(dat$x, dat$y)
  binxy <- data.frame(x = findInterval(points[,1], x),
                      y = findInterval(points[,2], y))

  binxy$x <- factor(binxy$x, levels = 1:nx)
  binxy$y <- factor(binxy$y, levels = 1:ny)

  # bin by sample for alpha (pvals)
  binxy$group <- factor(dat$sample)
  stab <- table(binxy)
  stab <- as.data.frame.table(stab)

  stab <- tidyr::pivot_wider(stab, names_from=group, values_from=Freq)
  rns <- paste(stab$x, stab$y, sep='-')
  stab$x <- stab$y <- NULL
  stab <- as.matrix(stab)
  row.names(stab) <- rns
  stab <- stab[rowSums(stab) != 0, ]

  uniq <- !duplicated(dat$sample)
  levs <- dat$group[uniq]
  names(levs) <- dat$sample[uniq]
  orig.ident <- levs[colnames(stab)]
  abd <- diff_abundance(stab, orig.ident=orig.ident, filter = FALSE)

  # get difference in number of cells between groups in nx*ny grid
  is.in <- gtab$group == 'test'
  d <- gtab[is.in, ]
  row.names(d) <- paste(d$x, d$y, sep='-')
  d$pval <- 1
  d[row.names(abd), 'pval'] <- abd$P.Value
  d$tots <- d$Freq + gtab$Freq[!is.in]
  d$diff <- d$Freq - gtab$Freq[!is.in]

  # scale difference by total number of cells
  not.zero <- d$tots != 0
  d$diff[not.zero] <- d$diff[not.zero]/d$tots[not.zero]

  xx <- x[-length(x)] + 0.5*diff(x)
  d$x <- xx[d$x]
  yy <- y[-length(y)] + 0.5*diff(y)
  d$y <- yy[d$y]

  # alpha 0.1: not significant
  # alpha 0.5-1: significant
  pt.dat <- d[d$tots > 0, ] # greater than 0 cells
  pt.dat$alpha <- 0.1
  is.sig <- pt.dat$pval < 0.05
  pt.dat$alpha[is.sig] <- range02(-log10(pt.dat$pval)[is.sig])

  pt.dat$direction <- "="
  pt.dat$direction[pt.dat$diff > 0] <- "↑"
  pt.dat$direction[pt.dat$diff < 0] <- "↓"
  pt.dat$direction <- factor(pt.dat$direction, levels = c('↑', '↓', '='), ordered = TRUE)

  return(pt.dat)
}


range02 <- function(x, newMax=1, newMin=0.5){
  (x - min(x))/(max(x)-min(x)) * (newMax - newMin) + newMin
}
