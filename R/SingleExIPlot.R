#' Plot a single expression by identity on a plot
#'
#' @param type Make either a 'ridge' or 'violin' plot
#' @param data Data to plot
#' @param idents Idents to use
#' @param sort Sort identity classes (on the x-axis) by the average
#' expression of the attribute being potted
#' @param y.max Maximum Y value to plot
#' @param adjust Adjust parameter for geom_violin
#' @param cols Colors to use for plotting
#' @param log plot Y axis on log scale
#' @param seed.use Random seed to use. If NULL, don't set a seed
#'
#' @return A ggplot-based Expression-by-Identity plot
#'
SingleExIPlot <- function(
  data,
  idents,
  hl = NULL,
  title = NULL,
  split = NULL,
  type = 'violin',
  sort = FALSE,
  y.max = NULL,
  adjust = 1,
  pt.size = 0.05,
  cols = NULL,
  seed.use = 42,
  log = FALSE,
  color = NULL,
  color_dark = NULL,
  nsel = 1,
  ncells = NULL
) {
  if (!is.null(seed.use)) {
    set.seed(seed.use)
  }
  if (!is.data.frame(data) || ncol(data) != 1) {
    stop("'SingleExIPlot requires a data frame with 1 column")
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
    warning(paste0("All cells have the same value of ", feature, "."))
  } else{
    data[, feature] <- data[, feature] + noise
  }
  axis.label <- ''
  y.max <- y.max %||% max(data[, feature][is.finite(x = data[, feature])])
  if (type == 'violin' && !is.null(split)) {
    data$split <- split
    vln.geom <- geom_violin
    fill <- 'split'
  } else if (type == 'splitViolin' && !is.null(x = split )) {
    data$split <- split
    vln.geom <- geom_split_violin
    fill <- 'split'
    type <- 'violin'
  } else {
    vln.geom <- ggplot2::geom_violin
    fill <- 'hl'
  }

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
  if (is.null(split)) {
    jitter_color <- c(color_dark, 'dimgray')[as.numeric(droplevels(hl))]
    jitter <- ggplot2::geom_jitter(height = 0,
                                   size = pt.size,
                                   shape = '.',
                                   color = jitter_color,
                                   alpha = 0.5,
                                   show.legend = FALSE)
  } else {
    jitter <- ggplot2::geom_jitter(
      position = ggplot2::position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9),
      size = pt.size,
      show.legend = FALSE
    )
  }
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
  if (pt.size > 0) {
    plot <- plot + jitter
  }
  if (!is.null(cols)) {
    if (!is.null(x = split)) {
      idents <- unique(x = as.vector(x = data$ident))
      splits <- unique(x = as.vector(x = data$split))
      labels <- if (length(x = splits) == 2) {
        splits
      } else {
        unlist(x = lapply(
          X = idents,
          FUN = function(pattern, x) {
            x.mod <- gsub(
              pattern = paste0(pattern, '.'),
              replacement = paste0(pattern, ': '),
              x = x,
              fixed = TRUE
            )
            x.keep <- grep(pattern = ': ', x = x.mod, fixed = TRUE)
            x.return <- x.mod[x.keep]
            names(x = x.return) <- x[x.keep]
            return(x.return)
          },
          x = unique(x = as.vector(x = data$split))
        ))
      }
      if (is.null(x = names(x = labels))) {
        names(x = labels) <- labels
      }
    } else {
      labels <- levels(x = droplevels(data$ident))
    }
    plot <- plot + ggplot2::scale_fill_manual(values = cols, labels = labels)
  }

  # add n cells
  lim <- max(data$x)
  lim <- lim+lim*0.2
  plot <- plot + ggplot2::ylim(c(NA, lim)) +
    ggplot2::annotate('text',
                      label=paste(ncells, 'cells'),
                      y=rep(lim, length(ncells)),
                      x=seq_along(ncells),
                      vjust = -0.3,
                      hjust = 'right',
                      color = 'dimgray',
                      size = 5)
  return(plot)
}
