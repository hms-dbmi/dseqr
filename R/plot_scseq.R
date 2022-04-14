
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
    dplyr::group_by(.data$cluster) %>%
    dplyr::sample_n(min(max.cells, dplyr::n())) %>%
    dplyr::pull(.data$id)

  scseq[, keep]
}


#' Get data for single-cell violin plots
#'
#' @param feature Feature name to generate violin plot for. Either a row or \code{colData} of \code{scseq}.
#' @param scseq \code{SingleCellExperiment}.
#' @param selected_cluster Name of the selected cluster.
#' @param by.sample if \code{TRUE} plot \code{feature} violin for each \code{scseq$batch}. Default (\code{FALSE})
#'  will plot \code{feature} for each \code{scseq$cluster}.
#' @param with.height Whether to return height of plot. See value for details.
#' @param decreasing if \code{TRUE}, violinlines with smaller mean values of \code{feature} will show up on top.
#'  Used to show features where smaller values indicate potential QC issues.
#'
#' @return list used by \link{VlnPlot}
#'
#' @keywords internal
get_violin_data <- function(feature, scseq, selected_cluster, by.sample = FALSE, decreasing = feature %in% c('ribo_percent', 'log10_sum', 'log10_detected'), with_all = FALSE, h5logs = NULL) {

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

  is.gene <- feature %in% row.names(scseq)

  title.type <- ifelse(is.gene, 'Expression', 'Value')

  # either highlight test group or selected cluster
  if (by.sample) {
    keep <- scseq$cluster %in% sel
    y <- factor(scseq$batch[keep])
    hl <- scseq$orig.ident[keep]

    # cluster name get's appended by VlnPlot
    title <- paste(title.type, 'by Sample:', feature, 'in')

  } else {
    hl <- as.integer(scseq$cluster)
    y <- factor(hl)
    hl[!hl %in% seli] <- 'out'
    hl <- factor(hl, levels = c(seli, 'out'))
    title <- paste(title.type, 'by Cluster:', feature)
  }


  if (is.gene) {
    if (!is.null(h5logs)) {
      x <- h5logs[feature, ]
      x <- x[colnames(scseq)]
    } else {
      x <- as.numeric(SingleCellExperiment::logcounts(scseq[feature, ]))
    }
    if (exists('keep')) x <- x[keep]

  } else {
    x <- scseq[[feature]]
  }

  mean.x <- tapply(x, y, mean)
  clus_ord <- order(mean.x, decreasing = decreasing)
  y <- factor(y, levels = levels(y)[clus_ord])

  df <- data.frame(x, hl, y) %>%
    dplyr::add_count(.data$y)

  ncells <- pct.cells <- NULL

  # show percent > 0 for genes
  # show number of cells in cluster otherwise
  if (is.gene) {
    get.pct <- function(x) round(sum(x > 0) / length(x) * 100)
    pct.cells <- tapply(df$x, as.character(df$y), get.pct)
    pct.cells <- pct.cells[levels(df$y)]
  } else {
    ncells <- tapply(df$x, as.character(df$y), length)
    ncells <- ncells[levels(df$y)]
  }


  # down sample to reduce plot size
  df <- df %>%
    dplyr::group_by(.data$y) %>%
    dplyr::mutate(nsamp = min(.data$n, 1000)) %>%
    dplyr::sample_n(.data$nsamp)

  res <- list(
    df = df,
    nsel = nsel,
    seli = seli,
    color = color,
    color_dark = color_dark,
    ncells = ncells,
    pct.cells = pct.cells,
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
#' @param pad_left Should cluster numbers be padded on the left? Used to align
#' numbers for violin plot.
#'
#' @return \code{annot} with cluster numbers pre-pended to non-numeric values.
#' @export
#'
#' @examples
#'
#' annot <- c('Monocytes', '2')
#' add_cluster_numbers(annot)
#'
add_cluster_numbers <- function(annot, pad_left = FALSE) {

  nums <- seq_along(annot)
  # use figure space
  if (pad_left) nums <- stringr::str_pad(nums, max(nchar(nums)), pad = '\U2007')

  is.char <- suppressWarnings(is.na(as.numeric(annot)))
  if (any(is.char)) {
    text <- annot[is.char]
    annot[is.char] <- paste0(nums[is.char], ': ', text)
  } else {
    annot <- nums
  }

  return(annot)
}

#' Violin Plot
#'
#' @param feature Character, gene name.
#' @param scseq \code{SingleCellExperiment} object.
#' @param selected_cluster Integer indicating which \code{level(scseq$cluster)} to highlight
#' (if \code{by.sample} is \code{FALSE}) or to subset data to (if \code{by.sample} is \code{TRUE}).
#' To plot by sample without subsetting to a cluster, set to \code{nlevels(scseq$cluster) + 1}.
#' @param by.sample Boolean indicating if violin plots should be by sample or by cluster (default).
#' @param with_all Set to \code{TRUE} if plotting by sample to get same colors as Shiny app.
#' @param with.height if \code{TRUE} returns a list object with plot and appropriate height.
#' @param violin_data Used for Shiny app. Result of \code{get_violin_data}.
#' Used when logcounts are seperate from \code{SingleCellExperiment} as in the Shiny app.
#' @param is_mobile If \code{TRUE} shortens up y-axis labels to fit mobile better.
#' @param hide_pct If \code{TRUE}, barplots indicating number of cells are hidden.
#' @param hide_ncells If \code{TRUE}, barplots indicating percent expression are hidden.
#' @param decreasing if \code{FALSE} (default) violins are ordered, bottom to top,
#' by increasing mean expression.
#'
#' @return \code{ggplot} or list
#'
plot_violin <- function(feature = NULL,
                        scseq = NULL,
                        selected_cluster = NULL,
                        by.sample = FALSE,
                        with_all = FALSE,
                        with.height = FALSE,
                        violin_data = NULL,
                        is_mobile = FALSE,
                        hide_pct = is_mobile,
                        hide_ncells = is_mobile,
                        decreasing = feature %in% c('ribo_percent',
                                                    'log10_sum',
                                                    'log10_detected')) {

  selected_cluster <- as.character(selected_cluster)

  # global variable bindings created by list2env
  df = nsel = seli = color = color_dark = ncells = pct.cells = clus_ord = clus_levs = title = NULL

  if (is.null(violin_data)) {
    violin_data <- get_violin_data(
      feature, scseq, selected_cluster, by.sample, decreasing, with_all)
  }

  list2env(violin_data, envir = environment())

  if (by.sample) {
    title <- paste(title, c(clus_levs, 'All Clusters')[seli])
    if (is_mobile) df <- shorten_y(df)

  } else {
    if (is_mobile) clus_levs <- seq_along(clus_levs)
    annot <- add_cluster_numbers(clus_levs, pad_left = TRUE)
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
  height <- length(levels(df$y)) * 38  + 130

  if (hide_pct) pct.cells <- NULL
  if (hide_ncells) ncells <- NULL

  # will color violins for all if not enough in highlight group
  n1.hl <- table(df$hl)
  n1.hl <- n1.hl[names(n1.hl) != 'out'] <= 1
  if (any(n1.hl)) {
    n1.names <- names(n1.hl)[n1.hl]
    df$hl[df$hl %in% n1.names] <- 'out'
    color <- color[!n1.hl]
    color_dark <- color_dark[!n1.hl]
    nsel <- max(0, nsel - sum(n1.hl))
  }

  pt.size <- df |>
    dplyr::group_by(.data$y) |>
    dplyr::mutate(ngt0 = sum(.data$x > 0)) |>
    dplyr::mutate(pt.size = ifelse(.data$ngt0 > 100, 0.01, 1)) |>
    dplyr::pull(.data$pt.size)

  pl <- SingleViolinIPlot(
    data,
    idents = df$y,
    hl = df$hl,
    pt.size = pt.size,
    title = title,
    color = color,
    color_dark = color_dark,
    nsel = nsel,
    ncells = ncells,
    pct.cells = pct.cells)


  if (with.height) pl <- list(plot = pl, height = height)
  return(pl)

}

# shorten labels for sample violin plot on mobile
shorten_y <- function(df) {
  short <- df %>%
    dplyr::group_by(.data$y) %>%
    dplyr::slice(1) %>%
    dplyr::group_by(.data$hl) %>%
    dplyr::arrange(as.character(.data$y)) %>%
    dplyr::mutate(new = paste0(toupper(.data$hl), seq_along(.data$hl))) %>%
    as.data.frame()

  row.names(short) <- short$y

  levels(df$y) <- short[levels(df$y), 'new']
  return(df)
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
  gene_dat <- sort(gene_dat, decreasing = TRUE)[seq(15)]
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
#' @param dark Use dark palettes? Default is \code{FALSE}.
#' @param with_all if \code{TRUE}, adds an additional level to get colors for.
#'  Default is \code{FALSE}
#'
#' @return Character vector with colour codes of \code{length(levs)}.
#'
#' @export
#' @examples
#' levs <- c('CD14 Mono', 'CD16 Mono', 'NK Cells')
#' get_palette(levs)
#'
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

    # remove shades of white
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


get_grid_expression <- function(gene, tts, grid) {

  tts <- lapply(tts, function(x) x[gene,, drop=FALSE])
  tts <- tts[!is.na(tts)]
  tt <- data.table::rbindlist(tts, fill = TRUE)
  tt <- as.data.frame(tt)
  row.names(tt) <- names(tts)
  tt <- tt[!is.na(tt$logFC), c('P.Value', 'logFC')]
  tt$cluster <- row.names(tt)

  # get points in grid with cells
  pt.dat <- grid %>%
    dplyr::add_count(.data$xi, .data$yi) %>%
    dplyr::left_join(tt) %>%
    dplyr::distinct() %>%
    dplyr::filter(.data$n > 3) %>% # at least n cells
    dplyr::mutate(logFC = ifelse(is.na(.data$logFC), 0, .data$logFC)) %>%
    dplyr::rename(pval = .data$P.Value, diff = .data$logFC) %>%
    dplyr::mutate(pval = replace(.data$pval, is.na(.data$pval), 1))

  return(pt.dat)
}

# checks multiple grid sizes and picks the one closest to have an average of
# target cells in non-empty cells
get_grid_size <- function(red.mat, target = 20) {

  # grid sizes to check
  dims <- list(
    c(120, 60),
    c(100, 50),
    c(80, 40),
    c(60, 30),
    c(40, 20)
  )

  xrange <- range(red.mat[,1])
  yrange <- range(red.mat[,2])

  df <- data.frame(x = red.mat[,1], y = red.mat[,2])
  avgs <- c()

  for (dim in dims) {
    nx <- dim[1]
    ny <- dim[2]

    # get number of cells in each non-empty box of grid
    counts <- df %>% dplyr::mutate(
      cut_x = cut(.data$x,
                  breaks = seq(xrange[1], xrange[2], length.out = nx),
                  include.lowest = TRUE),
      cut_y = cut(.data$y,
                  breaks = seq(yrange[1], yrange[2], length.out = ny),
                  include.lowest = TRUE)
    ) %>%
      dplyr::count(.data$cut_x, .data$cut_y) %>%
      dplyr::filter(.data$n > 0) %>%
      dplyr::pull(.data$n)


    avgs <- c(avgs, mean(counts))
  }

  # pick grid with closest to target as average
  chosen <- which.min(abs(avgs - target))
  return(dims[[chosen]])
}


get_grid_abundance <- function(scseq, group = scseq$orig.ident, sample = scseq$batch) {
  reds <- SingleCellExperiment::reducedDimNames(scseq)
  red <- reds[reds %in% c('UMAP', 'TSNE')]

  red.mat <- SingleCellExperiment::reducedDim(scseq, red)
  grid_size <- get_grid_size(red.mat)
  nx <- grid_size[1]
  ny <- grid_size[2]

  dat <- data.frame(x=red.mat[,1], y=red.mat[,2], group=group, sample=sample)

  # downsample within group so that each sample has same number of cells
  dsamp <- dat %>%
    dplyr::add_count(.data$sample) %>%
    dplyr::group_by(.data$group) %>%
    dplyr::mutate(nmin=min(.data$n)) %>%
    dplyr::group_by(.data$sample) %>%
    dplyr::sample_n(.data$nmin) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$n, -.data$nmin)

  # downsample so that each group has same number of cells
  nmin <- min(table(dsamp$group))
  dsamp <- dsamp %>%
    dplyr::group_by(.data$group) %>%
    dplyr::sample_n(nmin) %>%
    dplyr::ungroup()

  # create grid on non-downsampled
  x <- seq(min(dat$x), max(dat$x), length.out = nx)
  y <- seq(min(dat$y), max(dat$y), length.out = ny)

  # in nx*ny grid: get xy bin for downsampled points
  points <- cbind(dsamp$x, dsamp$y)
  binxy <- data.frame(x = findInterval(points[,1], x),
                      y = findInterval(points[,2], y))

  binxy$x <- factor(binxy$x, levels = seq_len(nx))
  binxy$y <- factor(binxy$y, levels = seq_len(ny))

  # bin by group for color (logFC sign)
  binxy$group <- factor(dsamp$group)
  gtab <- table(binxy)
  gtab <- as.data.frame.table(gtab)

  # in nx*ny grid: get xy bin for all points
  points <- cbind(dat$x, dat$y)
  binxy <- data.frame(x = findInterval(points[,1], x),
                      y = findInterval(points[,2], y))

  binxy$x <- factor(binxy$x, levels = seq_len(nx))
  binxy$y <- factor(binxy$y, levels = seq_len(ny))

  # bin by sample for alpha (pvals)
  binxy$group <- factor(dat$sample)
  stab <- table(binxy)
  stab <- as.data.frame.table(stab)

  stab <- tidyr::pivot_wider(stab, names_from = .data$group, values_from = .data$Freq)
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

  # get xy coords that define polygons in grid
  diff.x <- diff(x[c(1, 2)])
  diff.y <- diff(y[c(1, 2)])

  x2 <- x + diff.x
  d$x1 <- x[d$x]
  d$x2 <- x2[d$x]

  y2 <- y + diff.y
  d$y1 <- y[d$y]
  d$y2 <- y2[d$y]

  grid_abundance <- d[d$tots > 0, ] # greater than 0 cells
  return(grid_abundance)
}

add_grid_colors <- function(data) {

  # alpha 0.1: not significant
  # alpha 0.5-1: significant
  alpha <- rep(0.1, nrow(data))
  is.sig <- data$pval < 0.05
  alpha[is.sig] <- range02(-log10(data$pval)[is.sig])

  # for pval = 0
  alpha[is.infinite(alpha)] <- 1

  data$color <- '#FFFFFF'
  data$color[data$diff > 0] <- 'red'
  data$color[data$diff < 0] <- 'blue'
  data$color <- add_alpha(data$color, alpha)

  return(data)
}

add_alpha <- function(colors, alpha){
  if(missing(alpha)) stop("provide a value for alpha between 0 and 1")
  rgb <- grDevices::col2rgb(colors, alpha=TRUE)
  rgb[4,] <- round(rgb[4,]*alpha)
  new.colors <- rgb(rgb[1,], rgb[2,], rgb[3,], rgb[4,], maxColorValue = 255)
  return(new.colors)
}


range02 <- function(x, newMax=1, newMin=0.5) {
  scales::rescale(x, to = c(newMin, newMax))
}
