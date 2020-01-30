#' Dimensional reduction plot
#'
#' Graphs the output of a dimensional reduction technique on a 2D scatter plot where each point is a
#' cell and it's positioned based on the cell embeddings determined by the reduction technique. By
#' default, cells are colored by their identity class (can be changed with the group.by parameter).
#'
#' @param object SingleCellExperiment object
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param cells Vector of cells to plot (default is all cells)
#' @param cols Vector of colors, each color corresponds to an identity class. This may also be a single character
#' or numeric value corresponding to a palette as specified by \code{\link[RColorBrewer]{brewer.pal.info}}.
#' By default, ggplot2 assigns colors
#' @param pt.size Adjust point size for plotting
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca
#' @param group.by Name of one or more metadata columns to group (color) cells by
#' (for example, orig.ident); pass 'ident' to group by identity class
#' @param split.by Name of a metadata column to split plot by;
#' see \code{\link[Seurat]{FetchData}} for more details
#' @param shape.by If NULL, all points are circles (default). You can specify any
#' cell attribute (that can be pulled with FetchData) allowing for both
#' different colors and different shapes on cells
#' @param order Specify the order of plotting for the idents. This can be
#' useful for crowded plots if points of interest are being buried. Provide
#' either a full list of valid idents or a subset to be plotted last (on top)
#' @param label Whether to label the clusters
#' @param label.size Sets size of labels
#' @param repel Repel labels
#' @param cells.highlight A list of character or numeric vectors of cells to
#' highlight. If only one group of cells desired, can simply
#' pass a vector instead of a list. If set, colors selected cells to the color(s)
#' in \code{cols.highlight} and other cells black (white if dark.theme = TRUE);
#' will also resize to the size(s) passed to \code{sizes.highlight}
#' @param cols.highlight A vector of colors to highlight the cells as; will
#' repeat to the length groups in cells.highlight
#' @param sizes.highlight Size of highlighted cells; will repeat to the length
#' groups in cells.highlight
#' @param na.value Color value for NA points when using custom scale
#' @param combine Combine plots into a single gg object; note that if TRUE; themeing will not work when plotting multiple features
#' @param ncol Number of columns for display when combining plots
#' @param ... Extra parameters passed on to \code{\link{CombinePlots}}
#'
#' @return A ggplot object
#'
#' @importFrom rlang !! %||%
#'
#' @export
#'
#' @note For the old \code{do.hover} and \code{do.identify} functionality, please see
#' \code{HoverLocator} and \code{CellSelector}, respectively.
#'
#' @aliases TSNEPlot PCAPlot ICAPlot
#' @seealso \code{\link{FeaturePlot}} \code{\link{HoverLocator}}
#' \code{\link{CellSelector}} \code{link{FetchData}}
#'
#' @examples
#' DimPlot(object = pbmc_small)
#' DimPlot(object = pbmc_small, split.by = 'ident')
#'
DimPlot <- function(
  object,
  dims = c(1,2),
  cells = NULL,
  cols = NULL,
  pt.size = NULL,
  reduction = 'TSNE',
  group.by = NULL,
  split.by = NULL,
  shape.by = NULL,
  order = NULL,
  label = FALSE,
  label.size = 4,
  repel = FALSE,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50',
  combine = TRUE,
  ncol = NULL,
  label.highlight = NULL,
  ...
) {


  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }

  cells <- cells %||% colnames(x = object)
  data <- SingleCellExperiment::reducedDim(object, reduction)
  data <- as.data.frame(x = data)
  dims <- paste0(reduction, dims)
  group.by <- group.by %||% 'cluster'
  data[, group.by] <- SummarizedExperiment::colData(object)[cells, group.by]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by]]
  }
  if (!is.null(x = split.by)) {
    data[, split.by] <- object[[split.by]]
  }

  plots <- lapply(
    X = group.by,
    FUN = function(x) {
      plot <- SingleDimPlot(
        data = data[, c(dims, x, split.by, shape.by)],
        dims = dims,
        col.by = x,
        cols = cols,
        pt.size = pt.size,
        order = order,
        label = FALSE,
        cells.highlight = cells.highlight,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        na.value = na.value
      )
      if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = x,
          clusters = levels(data$cluster),
          repel = repel,
          size = label.size,
          split.by = split.by,
          label.highlight = label.highlight

        )
      }
      if (!is.null(x = split.by)) {
        plot <- plot + FacetTheme() +
          facet_wrap(
            facets = vars(!!sym(x = split.by)),
            ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
              length(x = unique(x = data[, split.by]))
            } else {
              ncol
            }
          )
      }
      return(plot)
    }
  )
  if (combine) {
    plots <- CombinePlots(
      plots = plots,
      ncol = if (!is.null(x = split.by) && length(x = group.by) > 1) {
        1
      } else {
        ncol
      },
      ...
    )
  }
  return(plots)
}

#' Plot a single dimension
#'
#' @param data Data to plot
#' @param dims A two-length numeric vector with dimensions to use
#' @param pt.size Adjust point size for plotting
#' @param col.by ...
#' @param cols Vector of colors, each color corresponds to an identity class. This may also be a single character
#' or numeric value corresponding to a palette as specified by \code{\link[RColorBrewer]{brewer.pal.info}}.
#' By default, ggplot2 assigns colors
#' @param shape.by If NULL, all points are circles (default). You can specify any cell attribute
#' (that can be pulled with FetchData) allowing for both different colors and different shapes on
#' cells.
#' @param order Specify the order of plotting for the idents. This can be useful for crowded plots if
#' points of interest are being buried. Provide either a full list of valid idents or a subset to be
#' plotted last (on top).
#' @param label Whether to label the clusters
#' @param repel Repel labels
#' @param label.size Sets size of labels
#' @param cells.highlight A list of character or numeric vectors of cells to
#' highlight. If only one group of cells desired, can simply
#' pass a vector instead of a list. If set, colors selected cells to the color(s)
#' in \code{cols.highlight} and other cells black (white if dark.theme = TRUE);
#'  will also resize to the size(s) passed to \code{sizes.highlight}
#' @param cols.highlight A vector of colors to highlight the cells as; will
#' repeat to the length groups in cells.highlight
#' @param sizes.highlight Size of highlighted cells; will repeat to the length
#' groups in cells.highlight
#' @param na.value Color value for NA points when using custom scale.
#'
SingleDimPlot <- function(
  data,
  dims,
  col.by = NULL,
  cols = NULL,
  pt.size = NULL,
  shape.by = NULL,
  order = NULL,
  label = FALSE,
  repel = FALSE,
  label.size = 4,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50'
) {

  pt.size <- pt.size %||% AutoPointSize(data = data)
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  if (!is.data.frame(x = data)) {
    data <- as.data.frame(x = data)
  }
  if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
    stop("Cannot find dimensions to plot in data")
  } else if (is.numeric(x = dims)) {
    dims <- colnames(x = data)[dims]
  }
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = sizes.highlight %||% pt.size,
      cols.highlight = cols.highlight,
      col.base = cols[1] %||% '#C3C3C3',
      pt.size = pt.size
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- 'highlight'
    pt.size <- highlight.info$size
    cols <- highlight.info$color
  }
  if (!is.null(x = order) && !is.null(x = col.by)) {
    if (typeof(x = order) == "logical") {
      if (order) {
        data <- data[order(data[, col.by]), ]
      }
    } else {
      order <- rev(x = c(
        order,
        setdiff(x = unique(x = data[, col.by]), y = order)
      ))
      data[, col.by] <- factor(x = data[, col.by], levels = order)
      new.order <- order(x = data[, col.by])
      data <- data[new.order, ]
      if (length(x = pt.size) == length(x = new.order)) {
        pt.size <- pt.size[new.order]
      }
    }
  }
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find ", col.by, " in plotting data, not coloring plot")
    col.by <- NULL
  } else {
    # col.index <- grep(pattern = col.by, x = colnames(x = data), fixed = TRUE)
    col.index <- match(x = col.by, table = colnames(x = data))
    if (grepl(pattern = '^\\d', x = col.by)) {
      # Do something for numbers
      col.by <- paste0('x', col.by)
    } else if (grepl(pattern = '-', x = col.by)) {
      # Do something for dashes
      col.by <- gsub(pattern = '-', replacement = '.', x = col.by)
    }
    colnames(x = data)[col.index] <- col.by
  }
  if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
    warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
  }
  plot <- ggplot2::ggplot(data = data) +
    ggplot2::geom_point(
      mapping = ggplot2::aes_string(
        x = dims[1],
        y = dims[2],
        color = paste0("`", col.by, "`"),
        shape = shape.by
      ),
      size = pt.size
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3))) +
    ggplot2::labs(color = NULL)
  if (label && !is.null(x = col.by)) {
    plot <- LabelClusters(
      plot = plot,
      id = col.by,
      repel = repel,
      size = label.size
    )
  }
  if (!is.null(x = cols)) {
    plot <- plot + if (length(x = cols) == 1 && (is.numeric(x = cols) || cols %in% rownames(x = brewer.pal.info))) {
      ggplot2::scale_color_brewer(palette = cols, na.value = na.value)
    } else {
      ggplot2::scale_color_manual(values = cols, na.value = na.value)
    }
  }
  plot <- plot + cowplot::theme_cowplot()
  return(plot)
}


#' Label clusters on a ggplot2-based scatter plot
#'
#' @param plot A ggplot2-based scatter plot
#' @param id Name of variable used for coloring scatter plot
#' @param clusters Vector of cluster ids to label
#' @param labels Custom labels for the clusters
#' @param split.by Split labels by some grouping label, useful when using
#' \code{\link[ggplot2]{facet_wrap}} or \code{\link[ggplot2]{facet_grid}}
#' @param repel Use \code{geom_text_repel} to create nicely-repelled labels
#' @param ... Extra parameters to \code{\link[ggrepel]{geom_text_repel}}, such as \code{size}
#'
#' @return A ggplot2-based scatter plot with cluster labels
#'
#' @export
#'
#' @seealso \code{\link[ggrepel]{geom_text_repel}} \code{\link[ggplot2]{geom_text}}
#'
#' @examples
#' plot <- DimPlot(object = pbmc_small)
#' LabelClusters(plot = plot, id = 'ident')
#'
LabelClusters <- function(
  plot,
  id,
  clusters = NULL,
  labels = NULL,
  split.by = NULL,
  repel = TRUE,
  label.highlight = NULL,
  ...
) {

  xynames <- unlist(x = GetXYAesthetics(plot = plot), use.names = TRUE)
  if (!id %in% colnames(x = plot$data)) {
    stop("Cannot find variable ", id, " in plotting data")
  }
  if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
    warning("Cannot find splitting variable ", id, " in plotting data")
    split.by <- NULL
  }
  data <- plot$data[, c(xynames, id, split.by)]
  possible.clusters <- as.character(x = na.omit(object = unique(x = data[, id])))
  groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[, id])))
  if (any(!groups %in% possible.clusters)) {
    stop("The following clusters were not found: ", paste(groups[!groups %in% possible.clusters], collapse = ","))
  }
  labels.loc <- lapply(
    X = groups,
    FUN = function(group) {
      data.use <- data[data[, id] == group, , drop = FALSE]
      data.medians <- if (!is.null(x = split.by)) {
        do.call(
          what = 'rbind',
          args = lapply(
            X = unique(x = data.use[, split.by]),
            FUN = function(split) {
              medians <- apply(
                X = data.use[data.use[, split.by] == split, xynames, drop = FALSE],
                MARGIN = 2,
                FUN = median,
                na.rm = TRUE
              )
              medians <- as.data.frame(x = t(x = medians))
              medians[, split.by] <- split
              return(medians)
            }
          )
        )
      } else {
        as.data.frame(x = t(x = apply(
          X = data.use[, xynames, drop = FALSE],
          MARGIN = 2,
          FUN = median,
          na.rm = TRUE
        )))
      }
      data.medians[, id] <- group
      return(data.medians)
    }
  )
  labels.loc <- do.call(what = 'rbind', args = labels.loc)
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop("Length of labels (", length(x = labels),  ") must be equal to the number of clusters being labeled (", length(x = labels.loc), ").")
  }
  names(x = labels) <- groups
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }

  label_colors <- 'black'
  if (!is.null(label.highlight)) {
    label_colors <- rep('dimgray', nrow(labels.loc))
    label_colors[label.highlight] <- 'black'
  }

  geom.use <- ifelse(test = repel, yes = ggrepel::geom_text_repel, no = ggplot2::geom_text)
  plot <- plot + geom.use(
    data = labels.loc,
    mapping = ggplot2::aes_string(x = xynames['x'], y = xynames['y'], label = id),
    color = label_colors,
    ...
  )
  return(plot)
}


#' Get X and Y aesthetics from a plot for a certain geom
#'
#' @param plot A ggplot2 object
#' @param geom Geom class to filter to
#' @param plot.first Use plot-wide X/Y aesthetics before geom-specific aesthetics
#'
#' @return A named list with values 'x' for the name of the x aesthetic and 'y' for the y aesthetic
#'
GetXYAesthetics <- function(plot, geom = 'GeomPoint', plot.first = TRUE) {
  geoms <- sapply(
    X = plot$layers,
    FUN = function(layer) {
      return(class(x = layer$geom)[1])
    }
  )
  geoms <- which(x = geoms == geom)
  if (length(x = geoms) == 0) {
    stop("Cannot find a geom of class ", geom)
  }
  geoms <- min(geoms)
  if (plot.first) {
    x <- as.character(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)[2]
    y <- as.character(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)[2]
  } else {
    x <- as.character(x = plot$layers[[geoms]]$mapping$x %||% plot$mapping$x)[2]
    y <- as.character(x = plot$layers[[geoms]]$mapping$y %||% plot$mapping$y)[2]
  }
  return(list('x' = x, 'y' = y))
}


#' Automagically calculate a point size for ggplot2-based scatter plots
#'
#' It happens to look good
#'
#' @param data A data frame being passed to ggplot2
#'
#' @return The "optimal" point size for visualizing these data
#'
#' @examples
#' df <- data.frame(x = rnorm(n = 10000), y = runif(n = 10000))
#' AutoPointSize(data = df)
#'
AutoPointSize <- function(data) {
  return(min(1583 / nrow(x = data), 1))
}


#' Combine ggplot2-based plots into a single plot
#'
#' @param plots A list of gg objects
#' @param ncol Number of columns
#' @param legend Combine legends into a single legend
#' choose from 'right' or 'bottom'; pass 'none' to remove legends, or \code{NULL}
#' to leave legends as they are
#' @param ... Extra parameters passed to plot_grid
#'
#' @return A combined plot
#'
#' @export
#'
#' @examples
#' pbmc_small[['group']] <- sample(
#'   x = c('g1', 'g2'),
#'   size = ncol(x = pbmc_small),
#'   replace = TRUE
#' )
#' plots <- FeaturePlot(
#'   object = pbmc_small,
#'   features = c('MS4A1', 'FCN1'),
#'   split.by = 'group',
#'   combine = FALSE
#' )
#' CombinePlots(
#'   plots = plots,
#'   legend = 'none',
#'   nrow = length(x = unique(x = pbmc_small[['group', drop = TRUE]]))
#' )
#'
CombinePlots <- function(plots, ncol = NULL, legend = NULL, ...) {
  plots.combined <- if (length(x = plots) > 1) {
    if (!is.null(x = legend)) {
      if (legend != 'none') {
        plot.legend <- cowplot::get_legend(plot = plots[[1]] + ggplot2::theme(legend.position = legend))
      }
      plots <- lapply(
        X = plots,
        FUN = function(x) {
          return(x + NoLegend())
        }
      )
    }
    plots.combined <- cowplot::plot_grid(
      plotlist = plots,
      ncol = ncol,
      align = 'hv',
      ...
    )
    if (!is.null(x = legend)) {
      plots.combined <- switch(
        EXPR = legend,
        'bottom' = cowplot::plot_grid(
          plots.combined,
          plot.legend,
          ncol = 1,
          rel_heights = c(1, 0.2)
        ),
        'right' = cowplot::plot_grid(
          plots.combined,
          plot.legend,
          rel_widths = c(3, 0.3)
        ),
        plots.combined
      )
    }
    plots.combined
  } else {
    plots[[1]]
  }
  return(plots.combined)
}



#' Removes the legend
#'
#' @param ... Extra parameters to be passed to \code{theme}
#'
#' @export
#' @examples
#' # Generate a plot with no legend
#' library(ggplot2)
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' p + NoLegend()
#'
NoLegend <- function(...) {
  no.legend.theme <- ggplot2::theme(
    # Remove the legend
    legend.position = 'none',
    # Validate the theme
    validate = TRUE,
    ...
  )
  return(no.legend.theme)
}

