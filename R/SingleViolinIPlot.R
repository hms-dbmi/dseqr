#' Violin plot by identity
#'
#' @param data Data to plot
#' @param idents Idents to use
#' @param hl Factor used to highlight specific \code{idents} via
#' fill, alpha, and color geoms.
#' @param title Plot title
#' @param sort Sort identity classes (on the x-axis) by the average
#' expression of the attribute being potted
#' @param y.max Maximum Y value to plot
#' @param adjust Adjust parameter for geom_violin
#' @param pt.size size parameter for geom_jitter
#' @param pt.shape shape parameter for geom_jitter
#' @param pct.cells Vector specifying percent of cells that express the marker,
#' once for each \code{levels(idents)}. If specified, an annotated bar plot is drawn.
#' @param ncells Vector specifying number of cells in each \code{levels(idents)}.
#' If specified, an annotated bar plot is drawn to indicate number of cells.
#' @param color Colors to use for fill with length equal to one less than
#' \code{length(levels(hl))}. The last level of \code{hl} is filled with gray.
#' @param color_dark Darker versions of \code{color} used for points and lines.
#' @param nsel Number of selected \code{idents}. Default is 1.
#' @param log plot Y axis on log scale
#' @param seed.use Random seed to use. If NULL, don't set a seed
#'
#' @return A ggplot-based Violin-by-Identity plot
#' @importFrom rlang %||% .data
#'
SingleViolinIPlot <- function(
  data,
  idents,
  hl = NULL,
  title = NULL,
  sort = FALSE,
  y.max = NULL,
  adjust = 1,
  pt.size = 1,
  pt.shape = 20,
  seed.use = 42,
  log = FALSE,
  color = NULL,
  color_dark = NULL,
  nsel = 1,
  ncells = NULL,
  pct.cells = NULL
) {
  if (!is.null(seed.use)) {
    set.seed(seed.use)
  }
  if (!is.data.frame(data) || ncol(data) != 1) {
    stop("'SingleViolinIPlot requires a data frame with 1 column")
  }
  feature <- 'x'
  data$ident <- idents
  data$hl <- hl
  if ((is.character(sort) && nchar(sort) > 0) || sort) {
    data$ident <- factor(
      x = data$ident,
      levels = names(x = rev(x = sort(
        x = tapply(
          X = data[, feature],
          INDEX = data$ident,
          FUN = mean
        ),
        decreasing = grepl(pattern = paste0('^', tolower(x = sort)), x = 'decreasing')
      )))
    )
  }
  if (log) {
    noise <- stats::rnorm(n = length(data[, feature])) / 200
    data[, feature] <- data[, feature] + 1
  } else {
    noise <- stats::rnorm(n = length(data[, feature])) / 100000
  }
  if (all(data[, feature] == data[, feature][1])) {
    warning("All cells have the same value of ", feature, ".")
  } else{
    data[, feature] <- data[, feature] + noise
  }
  axis.label <- ''
  y.max <- y.max %||% max(data[, feature][is.finite(x = data[, feature])])
  vln.geom <- ggplot2::geom_violin
  fill <- 'hl'

  x <- 'ident'
  y <- paste0("`", feature, "`")
  xlab <- ''
  ylab <- axis.label
  geom <- list(
    vln.geom(scale = 'width', adjust = adjust, trim = TRUE),
    ggplot2::scale_fill_manual(values = c(color, 'gray')),
    ggplot2::scale_alpha_manual(values = c(rep(0.4, nsel), 0.25)),
    ggplot2::scale_color_manual(values = c(color_dark, 'gray')),
    ggplot2::theme(legend.position = 'none',
                   plot.title.position = "plot",
                   axis.text.y = ggplot2::element_text(color = '#333333', size = 14, hjust = 0),
                   axis.text.x = ggplot2::element_text(color = '#333333', size = 14),
                   panel.grid.major.y = ggplot2::element_line(colour='dimgray', size=0.1),
                   plot.title = ggplot2::element_text(color = '#333333', size = 16, face = 'plain', margin = ggplot2::margin(b = 25))),
    ggplot2::coord_flip(),
    theme_dimgray()
  )
  jitter_color <- c(color_dark, 'dimgray')[as.numeric(droplevels(hl))]
  jitter <- ggplot2::geom_jitter(height = 0,
                                 size = pt.size,
                                 shape = pt.shape,
                                 width = 0.1,
                                 color = jitter_color,
                                 alpha = ifelse(pt.size < 0.1, 0.3, 0.5),
                                 show.legend = FALSE)
  log.scale <- ggplot2::scale_y_log10()
  axis.scale <- ggplot2::ylim


  (plot <- ggplot2::ggplot(
    data = data,
    mapping = ggplot2::aes_string(x = x, y = y, fill = fill, alpha=fill, color=fill)
  ) +
      ggplot2::labs(x = xlab, y = ylab, title = title, fill = NULL) +
      cowplot::theme_cowplot())
  (plot <- do.call(what = '+', args = list(plot, geom)))
  plot <- plot + if (log) {
    log.scale
  } else {
    axis.scale(min(data[, feature]), y.max)
  }
  if (pt.size[1] > 0) {
    plot <- plot + jitter
  }

  # set y limit and add mean value
  lim <- max(data$x)
  lim <- lim+lim*0.2

  plot <- plot +
    ggplot2::ylim(c(NA, lim)) +
    ggplot2::stat_summary(
      fun = "mean",
      geom = "crossbar",
      width = 0.3,
      colour = "black",
      size = .3)

  # add pct.cells/ncells labels or return
  if (!is.null(pct.cells)) {
    fraction.rect <- pct.cells/max(pct.cells)
    labels <- pct.cells
    subtitle <- 'Percent Expressed'

  } else if (!is.null(ncells)) {
    fraction.rect <- ncells/max(ncells)
    labels <- formatC(ncells, format="d", big.mark=",")
    subtitle <- 'Number of Cells'

  } else {
    return(plot)
  }

  # numbers
  nlabs <- length(labels)
  plot <- plot +
    ggplot2::annotate('text',
                      label=labels,
                      y=rep(lim, nlabs),
                      x=seq_len(nlabs),
                      vjust = -0.5,
                      hjust = 'right',
                      color = '#ACACAC',
                      size = 4)

  # bars to indicate fraction
  xmin <- diff(range(data$x))*0.1
  hl.levels <- levels(hl)[seq_along(color_dark)]
  data <- data[order(hl), ]
  hl.idents <- unique(data$ident[data$hl %in% hl.levels])
  hl.idx <- match(hl.idents, levels(idents))
  col.bars <- rep('dimgray', nlabs)
  col.bars[hl.idx] <- color_dark

  plot <- plot +
    ggplot2::annotate('rect',
                      ymin=lim-(fraction.rect*xmin),
                      ymax=rep(lim, nlabs),
                      xmin=seq_len(nlabs),
                      xmax=seq_len(nlabs)+0.07,
                      color = col.bars,
                      fill = col.bars)

  # use subtitle to indicate what numbers are
  plot <- plot +
    ggplot2::labs(subtitle = subtitle) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 1, color = 'dimgray', size = 14))

  return(plot)
}
