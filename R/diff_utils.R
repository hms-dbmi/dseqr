#' Linear model fitting of eset with limma.
#'
#' After selecting control and test samples for a contrast, surrogate variable
#' analysis (\code{\link[sva]{sva}}) and linear model fitting with \link[limma]{lmFit} is performed.
#'
#'
#' If analyses need to be repeated, previous results can be reloaded with \code{\link[base]{readRDS}}
#' and supplied to the \code{prev_anal} parameter. In this case, previous selections will be reused.
#'
#' @param eset Annotated eset created by \code{\link{load_seq}}.
#' @param data_dir String specifying directory of GSE folders.
#' @param annot String, column name in fData. For duplicated
#'   values in this column, the row with the highest interquartile range
#'   across selected samples will be kept. Appropriate values are \code{"SYMBOL"} (default - for gene level analysis)
#'   or \code{"ENTREZID_HS"} (for probe level analysis).
#' @param prev_anal Previous result of \code{run_limma}. If present, previous group
#'   selections will be reused.
#'
#' @export
#'
#' @return List with:
#'   \item{fit}{result of \code{\link[limma]{lmFit}}.}
#'   \item{mod}{\code{model.matrix} used for \code{fit}}
#'
#' @examples
#'
#' eset <- load_seq(data_dir)
#' lm_fit <- run_limma(eset)
#'
run_limma <- function (eset, annot = "SYMBOL", svobj = list('sv' = NULL), numsv = 0, prev_anal = NULL, filter = FALSE) {

  # check for annot column
  if (!annot %in% colnames(Biobase::fData(eset)))
    stop(annot, " column in fData(eset) is missing.")

  # determine if this is rna seq data
  rna_seq <- 'norm.factors' %in% colnames(Biobase::pData(eset))

  if (!is.null(prev_anal)) {
    # re-use selections from previous analysis
    prev_names <- row.names(prev_anal$pdata)
    stopifnot(all(prev_names %in% colnames(eset)))

    eset <- eset[, prev_names]
    Biobase::pData(eset)$group <- prev_anal$pdata$group

  } else {
    # select/add contrast
    eset <- select_contrast(eset)
  }

  # filtering low counts (as in tximport vignette)
  if (filter & rna_seq) eset <- filter_genes(eset)

  # add vsd element for cleaning
  eset <- add_vsd(eset, rna_seq = rna_seq)

  # add surrogate variable/pair adjusted ("clean") expression matrix for iqr_replicates
  eset <- add_adjusted(eset, svobj, numsv = numsv)

  # remove rows with duplicated/NA annot (SYMBOL or ENTREZID)
  eset <- iqr_replicates(eset, annot)

  # lmFit
  lm_fit <- fit_lm(eset = eset,
                   svobj = svobj,
                   numsv = numsv,
                   rna_seq = rna_seq)
  return (lm_fit)
}

#' Filter genes in RNA-seq ExpressionSet
#'
#' @param eset ExpressionSet with 'counts' assayDataElement and group column in pData
#'
#' @return filtered \code{eset}
#' @export
#' @keywords internal
#'
filter_genes <- function(eset) {
  counts <- Biobase::assayDataElement(eset, 'counts')
  keep <- edgeR::filterByExpr(counts, group = eset$group)
  eset <- eset[keep, ]
  if (!nrow(eset)) stop("No genes with reads after filtering")
  return(eset)
}

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

#' Run limma analysis.
#'
#' Runs limma differential expression analysis on all contrasts selected by
#' \code{add_contrast}. Analysis performed with and without surrogate
#' variables discovered by \code{diff_setup}. Also prints MDS plot and saves
#' results.
#'
#' @param eset Annotated eset created by \code{load_raw}. Replicate features and
#'   non-selected samples removed by \code{iqr_replicates}.
#'
#' @return list with slots:
#'   * \code{fit} Result of \link[limma]{lmFit}.
#'   * \code{mod} model matrix used for fit.
#'
fit_lm <- function(eset, svobj = list(sv = NULL), numsv = 0, rna_seq = TRUE){

  # setup model matrix with surrogate variables
  group <- Biobase::pData(eset)$group
  mod <- model.matrix(~0 + group)
  colnames(mod) <- gsub('^group', '', colnames(mod))
  svind <- seq_len(numsv)
  svmod <- svobj$sv[, svind, drop = FALSE]
  if (length(svind)) colnames(svmod) <- paste0('SV', svind)
  mod <- cbind(mod, svmod)

  lm_fit <- run_lmfit(eset, mod, rna_seq)

  # add enids for go/kegg pathway analyses
  lm_fit$fit$genes <- Biobase::fData(eset)[, 'ENTREZID', drop = FALSE]

  return(lm_fit)
}

#' Get top table
#'
#' @param lm_fit Result of \link{run_limma}
#' @param groups Test and Control group as strings.
#'
#' @return result of \link[limma]{toptable}
#' @export
get_top_table <- function(lm_fit, groups = c('test', 'ctrl'), with.es = TRUE, robust = TRUE) {
  contrast <- paste(make.names(groups[1]), make.names(groups[2]), sep = '-')

  ebfit <- fit_ebayes(lm_fit, contrast, robust = robust)
  tt <- limma::topTable(ebfit, coef = contrast, n = Inf, sort.by = 'p')
  if (with.es) tt <- add_es(tt, ebfit, groups = groups)

  return(tt)
}

#' Add expression data adjusted for pairs/surrogate variables
#'
#' @param eset ExpressionSet
#' @param svobj surrogate variable object
#' @param numsv Number of surrogate variables to adjust for
#' @param adj_path Path to file to restore/save adjuted data from/to for future speed
#'
#' @return eset with \code{adjusted} element added
#' @export
add_adjusted <- function(eset, svobj = list(sv = NULL), numsv = 0, adj_path = NULL) {


  if (!is.null(adj_path) && file.exists(adj_path)) {
    Biobase::assayDataElement(eset, 'adjusted') <- readRDS(adj_path)
    return(eset)
  }

  # get mods with group and pair effects
  mods <- get_mods(eset@phenoData)
  mod <- mods$mod
  mod0 <- mods$mod0

  # remove pairs from alternative model so that get cleaned
  pair_cols <- colnames(mod0)[-1]

  svs <- svobj$sv[, seq_len(numsv), drop = FALSE]

  mod <- mod[, !colnames(mod) %in% pair_cols]
  mod.clean <- cbind(mod0[, pair_cols], svs)

  # used DESeq::vsd transformed counts for RNA-Seq
  # for microarray this will be exprs slot
  y <- Biobase::assayDataElement(eset, 'vsd')

  adj <- clean_y(y, mod, mod.clean)
  Biobase::assayDataElement(eset, 'adjusted') <- adj
  if (!is.null(adj_path)) saveRDS(adj, adj_path)
  return(eset)
}



#' Get model matrices for surrogate variable analysis
#'
#' Used by \code{add_adjusted} to create model matrix with surrogate variables.
#'
#' @param eset Annotated eset with samples selected during \code{add_contrasts}.
#'
#' @return List with model matrix(mod) and null model matrix (mod0) used for \code{sva}.
#'
#' @export
get_mods <- function(pdata) {

  # make full and null model matrix
  group_levels <- setdiff(unique(pdata$group), NA)
  group <- factor(pdata$group, levels = group_levels)
  pair <- factor(pdata$pair)

  contrasts.fun <- function(l)lapply(l, contrasts, contrasts = FALSE)

  # if pairs include in alternative and null model
  if (length(pair)) {
    mod <- stats::model.matrix(~0 + group + pair, contrasts.arg = contrasts.fun(list(group=group, pair=pair)))
    mod0 <- stats::model.matrix(~1 + pair, data = group, contrasts.arg = contrasts.fun(list(pair=pair)))

  } else {
    mod <- stats::model.matrix(~0 + group)
    mod0 <- stats::model.matrix(~1, data = group)
  }

  # rename group columns
  colnames(mod)[1:length(group_levels)] <- group_levels

  if (length(pair)) {
    # remove non matched pairs
    pair_cols <- colnames(mod0)[-1]
    has.pair <- colSums(mod0[, pair_cols, drop = FALSE]) >= 2
    has.pair <- names(which(has.pair))

    # if multiple pair columns then remove first
    if (length(has.pair) > 1) has.pair <- has.pair[-1]

    mod <- mod[, c(group_levels, has.pair), drop = FALSE]
    mod0 <- mod0[, c("(Intercept)", has.pair), drop = FALSE]
  }

  return(list("mod" = mod, "mod0" = mod0))
}

model.matrix.dummy <- function(group, pair) {

}




#' Run surrogate variable analysis
#'
#' @param mods result of \code{get_mods}
#' @param eset Expression set.
#' @param rna_seq Boolean is it RNA-seq?
#'
#' @export
run_sva <- function(mods, eset, rna_seq = TRUE) {

  # remove duplicated rows (from 1:many PROBE:SYMBOL) as affect sva
  if (rna_seq) {
    PROBE <- Biobase::fData(eset)[,1]
  } else {
    PROBE <- Biobase::fData(eset)$PROBE
  }

  expr <- unique(data.table::data.table(Biobase::exprs(eset), PROBE))[, PROBE := NULL]
  expr <- as.matrix(expr)

  # sva or svaseq
  sva_fun <-ifelse(rna_seq, sva::svaseq, sva::sva)

  set.seed(100)
  svobj <- sva_fun(expr, mods$mod, mods$mod0)
  return(svobj)
}




#' Removes features with replicated annotation.
#'
#' For rows with duplicated annot, highested IQR retained.
#'
#' @param mod Model matrix without surrogate variables. generated by \code{diff_setup}.
#' @param svobj Result from \code{sva} function called during \code{diff_setup}.
#' @param annot feature to use to remove replicates.
#' @param rm.dup remove duplicates (same measure, multiple ids)?
#'
#' @return Expression set with unique features at probe or gene level.
#' @export
iqr_replicates <- function(eset, annot = "SYMBOL", rm.dup = FALSE, keep_path = NULL) {

  # for R CMD check
  iqrange = SYMBOL = NULL

  # do less work if possible as can take seconds
  fdata <- Biobase::fData(eset)
  annot.all <- Biobase::fData(eset)[, annot]
  annot.na  <- is.na(annot.all)
  annot.dup <- duplicated(annot.all[!annot.na])

  if (!is.null(keep_path) && file.exists(keep_path)) {
    keep <- readRDS(keep_path)
    eset <- eset[keep, ]
    Biobase::featureNames(eset) <- fdata[keep, annot]

  } else if (!any(annot.dup)) {
    eset <- eset[!annot.na, ]
    Biobase::featureNames(eset) <- fdata[!annot.na, annot]

  } else {
    adj <- Biobase::assayDataElement(eset, 'adjusted')

    # add inter-quartile ranges, row, and feature data to exprs data
    data <- as.data.frame(adj)
    data$iqrange <- matrixStats::rowIQRs(adj)
    data$row <- 1:nrow(data)
    data[, colnames(fdata)] <- fdata

    # remove rows with NA annot (occurs if annot is SYMBOL)
    data <- data[!is.na(data[, annot]), ]

    # for rows with same annot, keep highest IQR
    data <- data.table::data.table(data)
    data <- data[data[, .I[which.max(iqrange)], by = eval(annot)]$V1]

    # use row number to keep selected features
    if (!is.null(keep_path)) saveRDS(data$row, keep_path)
    eset <- eset[data$row, ]

    # use annot for feature names
    Biobase::featureNames(eset) <- Biobase::fData(eset)[, annot]
  }

  if (rm.dup) {
    not.dup <- !duplicated(Biobase::assayDataElement(eset, 'adjusted'))
    eset <- eset[not.dup, ]
  }

  return(eset)
}


format_scaling <- function(scaling, adj, group, exprs) {

  scaling %>%
    dplyr::rename('MDS1' = V1, 'MDS2' = V2) %>%
    dplyr::mutate(Sample = row.names(exprs)) %>%
    dplyr::mutate(Group = group) %>%
    dplyr::mutate(Group =  dplyr::recode(Group, ctrl = 'Control', test = 'Test')) %>%
    dplyr::mutate(Title = ifelse(adj, 'adjusted', 'not adjusted'))
}




#' Fit ebayes model
#'
#' @param lm_fit Result of call to \link{run_limma}
#' @param contrasts Character vector of contrasts to fit.
#'
#' @return result of \link[limma]{eBayes}
#' @export
fit_ebayes <- function(lm_fit, contrasts, robust = TRUE) {
  colnames(lm_fit$fit$coefficients) <- colnames(lm_fit$mod) <- make.names(colnames(lm_fit$mod))
  contrast_matrix <- limma::makeContrasts(contrasts = contrasts, levels = lm_fit$mod)

  eb_fit <- limma::contrasts.fit(lm_fit$fit, contrast_matrix)
  eb_fit <- limma::eBayes(eb_fit, robust = robust)
  return (eb_fit)
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
#' @param scaling tibble with columns MDS1 and MDS2 corresponding to differential expression without SVA
#' @param scaling_adj tibble with columns MDS1 and MDS2 corresponding to differential expression with SVA/pairs
#'
#' @return plotly object
#' @export
plotlyMDS <- function(scaling, scaling_adj, group_colors = c('#337ab7', '#e6194b'), adjusted = FALSE) {

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
        title = list(text = 'Sammon MDS plots', x = 0.08, y = 0.98),
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


#' Perform lmFit analysis from limma.
#'
#' If paired samples, runs \code{\link{duplicateCorrelation}} to estimate intra-patient variance.
#'
#' @param eset Annotated eset created by \code{load_raw}. Non-selected samples
#'    and duplicate features removed by \code{add_contrasts} and
#'    \code{iqr_replicates}.
#' @param mod Model matrix generated by \code{diff_setup}. With
#'   or without surrogate variables.
#' @param rna_seq Is this an RNA-seq \code{eset}? Default is \code{TRUE}.
#'
#' @return result from call to limma \code{lmFit}.
run_lmfit <- function(eset, mod, rna_seq = TRUE) {

  pdata <- Biobase::pData(eset)
  pair <- pdata$pair
  y <- Biobase::exprs(eset)
  if (rna_seq) lib.size <- pdata$lib.size * pdata$norm.factors

  if (length(pair) & rna_seq) {
    # rna-seq paired
    # see https://support.bioconductor.org/p/110780/ for similar

    # first round
    v <- limma::voomWithQualityWeights(y, mod, lib.size = lib.size)
    corfit <- limma::duplicateCorrelation(v, mod, block = pair)

    # if couldn't estimate within-block correlation, model pair as fixed effect
    if (is.nan(corfit$consensus.correlation)) {
      fit_fun <- function() {
        v <- limma::voomWithQualityWeights(y, mod, lib.size = lib.size)
        limma::lmFit(v, design = mod)
      }

      mod <- get_mods(eset@phenoData)$mod
      fit <- fit_fun()

      # if no dof, drop pairs and retry
      if (fit$df.residual[1] == 0) {
        eset$pair <- NULL
        mod <- get_mods(eset@phenoData)$mod
        fit <- fit_fun()
      }

    } else {
      # second round
      v <- limma::voomWithQualityWeights(y, mod, lib.size = lib.size, block = pair, correlation = corfit$consensus.correlation, plot = TRUE)
      corfit <- limma::duplicateCorrelation(v, mod, block = pair)

      fit <- limma::lmFit(v, mod, correlation = corfit$consensus.correlation, block = pair)
    }

  } else if (rna_seq) {
    # rna-seq not paired
    # get normalized lib size and voom
    v <- limma::voomWithQualityWeights(y, mod, lib.size = lib.size)
    fit  <- limma::lmFit(v, design = mod)

  } else if (length(pair) & !rna_seq) {
    # microarray paired
    corfit <- limma::duplicateCorrelation(y, mod, block = pair)

    # if couldn't estimate within-block correlation, model pair as fixed effect
    if (is.nan(corfit$consensus.correlation)) {
      fit <- limma::lmFit(y, mod)

      # if no dof, drop pairs and retry
      if (fit$df.residual == 0) {
        eset$pair <- NULL
        mod <- get_mods(eset@phenoData)$mod
        fit <- limma::lmFit(y, mod)
      }

    } else {
      fit <- limma::lmFit(y, mod, correlation = corfit$consensus.correlation, block = pair)
    }

  } else {
    # microarray not paired
    fit <- limma::lmFit(y, mod)
  }

  return(list(fit = fit, mod = mod))
}


#' Adjusts expression data for surrogate variables.
#'
#' Factors out effect of surrogate variables discovered during surrogate variable
#' analysis.
#'
#' @param y Expression data of eset.
#' @param mod Full model matrix supplied to \code{sva}.
#' @param svs Surrogate variables returned by \code{sva} (svobj$sv).
#'
#' @return Expression data with effects of svs removed.
clean_y <- function(y, mod, mod.clean) {

  # if no factors to clean return original
  if (!ncol(mod.clean)) return(y)

  X = cbind(mod, mod.clean)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}
