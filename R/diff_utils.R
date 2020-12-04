#' Save lmfit result to disk
#'
#' @param lm_fit Result of run_limma or run_limma_scseq
#' @param dataset_dir directory to save results in
#' @param numsv Number of surrogate variables modeled. Default is 0.
#' @param anal_suffix suffix to append to saved name.
#'
#' @return NULL
#' @export
#' @keywords internal
save_lmfit <- function(lm_fit, dataset_dir, numsv = 0, anal_suffix = '') {

  if (nchar(anal_suffix)) anal_suffix <- paste0(anal_suffix, '_')
  fit_name <- paste0('lm_fit_', anal_suffix, paste0(numsv, 'svs.rds'))
  fit_path <- file.path(dataset_dir, fit_name)
  saveRDS(lm_fit, fit_path)
}



format_scaling <- function(scaling, adj, group, exprs) {

  scaling %>%
    dplyr::rename('MDS1' = V1, 'MDS2' = V2) %>%
    dplyr::mutate(Sample = row.names(exprs)) %>%
    dplyr::mutate(Group = group) %>%
    dplyr::mutate(Group =  dplyr::recode(Group, ctrl = 'Control', test = 'Test')) %>%
    dplyr::mutate(Title = ifelse(adj, 'adjusted', 'not adjusted'))
}



#' Get scalings for MDS plots
#'
#' For interactive MDS plot of expression values with and without surrogate variable analysis.
#'
#' @param exprs \code{matrix} of expression values.
#' @param adj \code{matrix} of expression values with surrogate variables/pairs regressed out.
#' @param group Character vector with values \code{'control'} and \code{'test'} indicating group membership.
#' @importFrom magrittr "%>%"
#'
#' @return List of tibbles with MDS scalings with and without SVA
#' @export
get_mds <- function(exprs, adj, group) {

  # get_dist acts on rows
  exprs <- t(exprs[complete.cases(exprs), ])
  adj <- t(adj[complete.cases(adj), ])

  dist <- get_dist(exprs, method = 'spearman')
  dist_adj <- get_dist(adj, method = 'spearman')

  # sammon scaling for dimensionality reduction
  capture.output({
    scaling <- tibble::as_tibble(MASS::sammon(dist, k = 2)$points) %>%
      format_scaling(adj = FALSE, group = group, exprs = exprs)

    scaling_adj <- tibble::as_tibble(MASS::sammon(dist_adj, k = 2)$points) %>%
      format_scaling(adj = TRUE, group = group, exprs = exprs)
  })

  return(list(scaling = scaling, scaling_adj = scaling_adj))
}

# file MASS/R/sammon.R
# copyright (C) 1994-2005 W. N. Venables and B. D. Ripley
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
sammon <- function(d, y = cmdscale(d, k), k = 2, niter = 100, trace = TRUE, magic = 0.2, tol = 1e-4) {
  call <- match.call()
  if(any(is.infinite(d))) stop("Infs not allowed in 'd'")
  if(any(is.na(d)) && missing(y))
    stop("an initial configuration must be supplied if there are NAs in 'd'")
  if(!is.matrix(y)) stop("'y' must be a matrix")

  if(is.null(n <- attr(d, "Size"))) {
    x <- as.matrix(d)
    if((n <- nrow(x)) != ncol(x))
      stop("distances must be result of 'dist' or a square matrix")
    rn <- rownames(x)
  } else {
    x <- matrix(0, n, n)
    x[row(x) > col(x)] <- d
    x <- x + t(x)
    rn <- attr(d, "Labels")
  }
  n <- as.integer(n)
  if(is.na(n)) stop("invalid size")
  ab <- x[row(x) < col(x)] <= 0
  if (any(ab, na.rm = TRUE)) {
    ab <- !is.na(ab) & ab
    aa <- cbind(as.vector(row(x)), as.vector(col(x)))[row(x) < col(x),]
    aa <- aa[ab, , drop=FALSE]
    stop(gettextf("zero or negative distance between objects %d and %d",
                  aa[1,1], aa[1,2]), domain = NA)
  }
  nas <- is.na(x)
  diag(nas) <- FALSE  # diag never used
  if(any(rowSums(!nas) < 2)) stop("not enough non-missing data")

  if(any(dim(y) != c(n, k)) ) stop("invalid initial configuration")
  if(any(!is.finite(y))) stop("initial configuration must be complete")
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"
  z <- .C(VR_sammon,
          x = x,
          n,
          as.integer(k),
          y = y,
          as.integer(niter),
          e = double(1),
          as.integer(trace),
          as.double(magic),
          as.double(tol),
          NAOK = TRUE)
  points <- z$y
  dimnames(points) <- list(rn, NULL)
  list(points=points, stress=z$e, call=call)
}

#' Enhanced Distance Matrix Computation and Visualization
#' @description Clustering methods classify data samples into groups of similar
#'   objects. This process requires some methods for measuring the distance or
#'   the (dis)similarity between the observations. Read more:
#'   \href{http://www.sthda.com/english/wiki/clarifying-distance-measures-unsupervised-machine-learning}{STHDA
#'    website - clarifying distance measures.}. \itemize{ \item get_dist():
#'   Computes a distance matrix between the rows of a data matrix. Compared to
#'   the standard \code{\link[stats]{dist}}() function, it supports
#'   correlation-based distance measures including "pearson", "kendall" and
#'   "spearman" methods. \item fviz_dist(): Visualizes a distance matrix }
#' @param x a numeric matrix or a data frame.
#' @param method the distance measure to be used. This must be one of
#'   "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski",
#'   "pearson", "spearman" or "kendall".
#' @param stand logical value; default is FALSE. If TRUE, then the data will be
#'   standardized using the function scale(). Measurements are standardized for
#'   each variable (column), by subtracting the variable's mean value and
#'   dividing by the variable's standard deviation.
#' @param ... other arguments to be passed to the function dist() when using get_dist().
#' @return \itemize{ \item get_dist(): returns an object of class "dist". \item
#'   fviz_dist(): returns a ggplot2 }
#' @seealso \code{\link[stats]{dist}}
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @examples
#' data(USArrests)
#' res.dist <- get_dist(USArrests, stand = TRUE, method = "pearson")
#'
#' fviz_dist(res.dist,
#'    gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
#' @name dist
#' @rdname dist
#' @export
get_dist <- function(x, method = "euclidean",  stand = FALSE, ...){

  if(stand) x <- scale(x)
  if(method %in% c("pearson", "spearman", "kendall")){
    res.cor <- stats::cor(t(x),  method = method)
    res.dist <- stats::as.dist(1 - res.cor, ...)
  }
  else res.dist <- stats::dist(x, method = method, ...)

  res.dist
}

#' Plot MDS plotlys
#'
#' @param scaling tibble with coordinate columns MDS1 and MDS2 calculated from expression data without correction for surrogate variables
#' @param scaling_adj Optional. Same as \code{scaling} but using expression data adjusted for surrogate variables.
#'   If omitted, an MDS plot is created only for \code{scaling}.
#' @param group_colors colors to use, one for each unique groups in \code{scaling$Group}.
#' @param title Plot title.
#'
#' @return plotly object
#' @export
plotlyMDS <- function(scaling, scaling_adj = NULL, group_colors = c('#337ab7', '#e6194b'), adjusted = FALSE, title = 'Sammon MDS plots') {

  if(is.null(scaling)) return(NULL)
  # make x and y same range
  xrange <- range(c(scaling$MDS1, scaling_adj$MDS1))
  yrange <- range(c(scaling$MDS2, scaling_adj$MDS2))

  addx <- diff(xrange) * .10
  addy <- diff(yrange) * .10

  xrange <- xrange + c(-addx, addx)
  yrange <- yrange + c(-addy, addy)


  shapes <- list(
    type = "rect",
    x0 = 0,
    x1 = 1,
    xref = "paper",
    y0 = 0,
    y1 = 16,
    yanchor = 1,
    yref = "paper",
    ysizemode = "pixel",
    fillcolor = '#cccccc',
    line = list(color = "#cccccc"))

  margin <- list(t = 60, l = 10, r = 10,  b = 10)

  xaxis <-list(title = 'MDS 1', zeroline = FALSE, showticklabels = FALSE, range = xrange,
                        linecolor = '#cccccc', mirror = TRUE, linewidth = 1)

  yaxis <- list(title = 'MDS 2', zeroline = FALSE, showticklabels = FALSE, range = yrange,
               linecolor = '#cccccc', mirror = TRUE, linewidth = 1)


  if (!adjusted) {

    pl <- plotly::plot_ly(scaling,
                          x = ~MDS1,
                          y = ~MDS2,
                          customdata = ~Group,
                          color = ~Group,
                          colors = group_colors,
                          showlegend = FALSE,
                          hovertemplate = paste0(
                            '<b>Group</b>: %{customdata}<br>',
                            '<b>Sample</b>: %{text}',
                            '<extra></extra>')) %>%
      plotly::add_markers(text = ~Sample, hoverinfo = 'text') %>%
      plotly::layout(
        title = list(text = title, x = 0.08, y = 0.98),
        margin = margin,
        shapes = shapes,
        xaxis = xaxis,
        yaxis = yaxis,
        annotations = list(
          list(x = 0.5 , y = 1.055, text = "Not Adjusted", showarrow = F, xref='paper', yref='paper'))
        ) %>%
      plotly::config(displaylogo = FALSE, displayModeBar = FALSE)

  } else {

    pl <- plotly::plot_ly(scaling_adj,
                          x = ~MDS1,
                          y = ~MDS2,
                          customdata = ~Group,
                          color = ~Group,
                          colors = group_colors,
                          showlegend = FALSE,
                          hovertemplate = paste0(
                            '<b>Group</b>: %{customdata}<br>',
                            '<b>Sample</b>: %{text}',
                            '<extra></extra>')) %>%
      plotly::add_markers(text = ~Sample, hoverinfo = 'text') %>%
      plotly::layout(
        margin = margin,
        shapes = shapes,
        xaxis = xaxis,
        yaxis = yaxis,
        annotations = list(
          list(x = 0.5 , y = 1.055, text = "Adjusted", showarrow = F, xref='paper', yref='paper'))
      ) %>%
      plotly::config(displaylogo = FALSE,
                     displayModeBar = 'hover',
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
                                                'toggleSpikelines'))
  }

  return(pl)
}
