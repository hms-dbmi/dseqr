#' Differential expression analysis of eset.
#'
#' After selecting control and test samples for a contrast, surrogate variable
#' analysis (\code{\link[sva]{sva}}) and differential expression analysis is performed.
#'
#' To specify a contrast, first select rows for control samples, and click the \emph{Add Control Rows}
#' button. Repeat for test samples. After control and test samples have been added
#' click the \emph{Done} button. If you make a mistake, click the \emph{Reset} button.
#'
#' Analysis results are saved in the corresponding in \code{data_dir} as "diff_expr.rds". If analyses
#' needs to be repeated, previous results can be reloaded with \code{\link[base]{readRDS}}
#' and supplied to the \code{prev_anal} parameter. In this case, previous selections will be reused.
#'
#' @param eset Annotated eset created by \code{\link{load_seq}}.
#' @param data_dir String specifying directory of GSE folders.
#' @param annot String, column name in fData. For duplicated
#'   values in this column, the row with the highest interquartile range
#'   across selected samples will be kept. Appropriate values are \code{"SYMBOL"} (default - for gene level analysis)
#'   or \code{"ENTREZID_HS"} (for probe level analysis).
#' @param prev_anal Previous result of \code{\link{diff_expr}}. If present, previous group
#'   selections will be reused.
#'
#' @export
#'
#' @return List of named lists, one for each GSE. Each named list contains:
#'   \item{pdata}{data.frame with phenotype data for selected samples.
#'      Columns \code{treatment} ('ctrl' or 'test'), \code{group}, and \code{pairs} are
#'      added based on user selections.}
#'   \item{top_tables}{List with results of \code{\link[limma]{topTable}} call (one per
#'      contrast). These results account for the effects of nuissance variables
#'      discovered by surrogate variable analysis.}
#'   \item{ebayes_sv}{Results of call to \code{\link[limma]{eBayes}} with surrogate
#'      variables included in the model matrix.}
#'   \item{annot}{Value of \code{annot} variable.}
#'
#' @examples
#'
#' data_dir <- system.file('extdata', 'IBD', package='drugseqr', mustWork = TRUE)
#' eset <- load_seq(data_dir)
#' anal <- diff_expr(eset, data_dir, anal_name = 'IBD')
run_limma <- function (eset, data_dir, annot = "SYMBOL", svobj = list('sv' = NULL), numsv = 0, prev_anal = NULL) {

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

  # add vsd element for cleaning
  eset <- add_vsd(eset)

  # add surrogate variable/pair adjusted ("clean") expression matrix for iqr_replicates
  eset <- add_adjusted(eset, svobj, numsv = numsv)

  # remove rows with duplicated/NA annot (SYMBOL or ENTREZID)
  eset <- iqr_replicates(eset, annot)

  # lmFit
  lm_fit <- fit_lm(eset = eset,
                   data_dir = data_dir,
                   svobj = svobj,
                   numsv = numsv,
                   rna_seq = rna_seq)
  return (lm_fit)
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
#' @param data_dir String, path to directory with raw data.
#' @param rna_seq Is the analysis for RNA-seq? Auto inferred in \code{diff_expr}.
#'
#' @return list with slots \itemize{
#'   \fit Result of \link[limma]{lmFit}.
#'   \mod model matrix used for fit.
#' }
fit_lm <- function(eset, data_dir, svobj = list(sv = NULL), numsv = 0, rna_seq = TRUE){

  # setup model matrix with surrogate variables
  group <- Biobase::pData(eset)$group
  mod <- model.matrix(~0 + group)
  colnames(mod) <- gsub('^group', '', colnames(mod))
  svind <- seq_len(numsv)
  svmod <- svobj$sv[, svind, drop = FALSE]
  if (length(svind)) colnames(svmod) <- paste0('SV', svind)
  mod <- cbind(mod, svmod)

  lm_fit <- run_lmfit(eset, mod, rna_seq)

  # save lm_fit result (slow and can re-use for other contrasts)
  # fit_ebayes is also input into goana/kegga
  fit_name <- paste('lm_fit', paste0(numsv, 'svs.rds'), sep = '_')
  fit_path <- file.path(data_dir, fit_name)
  if (!file.exists(fit_path)) saveRDS(lm_fit, fit_path)

  return(lm_fit)
}

#' Get top table
#'
#' @param lm_fit Result of \link{run_limma}
#' @param contrast String specifying contrast.
#'
#' @return result of \link[limma]{topTable}
#' @export
get_top_table <- function(lm_fit, groups = c('test', 'ctrl')) {
  contrast <- paste(groups[1], groups[2], sep = '-')
  ebfit <- fit_ebayes(lm_fit, contrast)
  tt <- limma::topTable(ebfit, coef = contrast, n = Inf)
  tt <- add_es(tt, ebfit, groups = groups)
  return(tt)
}

#' Add expression data adjusted for pairs/surrogate variables
#'
#' @param eset
#' @param mods
#' @param svobj
#' @param numsv
#'
#' @return eset with \code{adjusted} element added
#' @export
add_adjusted <- function(eset, svobj = list(sv = NULL), numsv = 0, adj_path = NULL) {


  if (!is.null(adj_path) && file.exists(adj_path)) {
    Biobase::assayDataElement(eset, 'adjusted') <- readRDS(adj_path)
    return(eset)
  }

  # get mods with group and pair effects
  mods <- get_mods(eset)
  mod <- mods$mod
  mod0 <- mods$mod0

  # remove pairs from alternative model so that get cleaned
  pair_cols <- colnames(mods$mod0)[-1]

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
#' Used by \code{diff_expr} to create model matrix with surrogate variables
#' in order to run \code{diff_anal}.
#'
#' @param eset Annotated eset with samples selected during \code{add_contrasts}.
#' @param svanal Perform surrogate variable analysis? Default is \code{TRUE}.
#' @param rna_seq Is this an RNA-seq experiment? Inferred from \code{eset} in \code{\link{diff_expr}}.
#'
#' @seealso \code{\link{add_contrast}}, \code{\link{diff_expr}}.
#' @return List with model matrix(mod), model matrix with surrogate
#'         variables(modsv), and result of \code{sva} function.
#'
#' @export
get_mods <- function(eset) {

  # make full and null model matrix
  pdata <- Biobase::pData(eset)
  group_levels <- setdiff(unique(pdata$group), NA)
  group <- factor(pdata$group, levels = group_levels)
  pair <- factor(pdata$pair)

  # if pairs include in alternative and null model
  if (length(pair)) {
    mod <- stats::model.matrix(~0 + group + pair)
    mod0 <- stats::model.matrix(~1 + pair, data = group)

  } else {
    mod <- stats::model.matrix(~0 + group)
    mod0 <- stats::model.matrix(~1, data = group)
  }

  # rename group columns
  colnames(mod)[1:length(group_levels)] <- group_levels

  if (length(pair)) {
    # remove non matched pairs
    pair_cols <- colnames(mod0)[-1]
    has.pair <- colSums(mod0[, pair_cols]) >= 2
    has.pair <- names(which(has.pair))

    mod <- mod[, c(group_levels, has.pair), drop = FALSE]
    mod0 <- mod0[, c("(Intercept)", has.pair), drop = FALSE]
  }

  return(list("mod" = mod, "mod0" = mod0))
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

  svobj <- sva_fun(expr, mods$mod, mods$mod0)
  return(svobj)
}




#' Removes features with replicated annotation.
#'
#' For rows with duplicated annot, highested IQR retained.
#'
#' @inheritParams diff_expr
#' @inheritParams diff_setup
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
fit_ebayes <- function(lm_fit, contrasts) {
  contrast_matrix <- limma::makeContrasts(contrasts = contrasts, levels = lm_fit$mod)
  eb_fit <- limma::contrasts.fit(lm_fit$fit, contrast_matrix)
  eb_fit <- limma::eBayes(eb_fit, robust = TRUE)
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

  dist <- factoextra::get_dist(exprs, method = 'spearman')
  dist_adj <- factoextra::get_dist(adj, method = 'spearman')

  # sammon scaling for dimensionality reduction
  capture.output({
    scaling <- tibble::as_tibble(MASS::sammon(dist, k = 2)$points) %>%
      format_scaling(adj = FALSE, group = group, exprs = exprs)

    scaling_adj <- tibble::as_tibble(MASS::sammon(dist_adj, k = 2)$points) %>%
      format_scaling(adj = TRUE, group = group, exprs = exprs)
  })

  return(list(scaling = scaling, scaling_adj = scaling_adj))
}

#' Plot MDS plotlys
#'
#' @param scaling tibble with columns MDS1 and MDS2 corresponding to differential expression without SVA
#' @param scaling_adj tibble with columns MDS1 and MDS2 corresponding to differential expression with SVA/pairs
#'
#' @return plotly object
#' @export
plotlyMDS <- function(scaling, scaling_adj, group_colors = c('#337ab7', '#e6194b')) {

  if(is.null(scaling)) return(NULL)
  # make x and y same range
  xrange <- range(c(scaling$MDS1, scaling_adj$MDS1))
  yrange <- range(c(scaling$MDS2, scaling_adj$MDS2))

  addx <- diff(xrange) * .10
  addy <- diff(yrange) * .10

  xrange <- xrange + c(-addx, addx)
  yrange <- yrange + c(-addy, addy)

  # legend styling
  l <- list(
    font = list(
      family = "sans-serif",
      size = 12,
      color = "#000"),
    bgcolor = "#f8f8f8",
    bordercolor = "#e7e7e7",
    borderwidth = 1)


  p1 <- plotly::plot_ly(scaling, x = ~MDS1, y = ~MDS2, color = ~Group, colors = group_colors, showlegend = FALSE) %>%
    plotly::add_markers(text = ~Sample, hoverinfo = 'text') %>%
    plotly::layout(
      xaxis = list(title = 'MDS 1', zeroline = FALSE, showticklabels = FALSE, range = xrange,
                   linecolor = '#cccccc', mirror = TRUE, linewidth = 1),
      yaxis = list(title = 'MDS 2', zeroline = FALSE, showticklabels = FALSE, range = yrange,
                   linecolor = '#cccccc', mirror = TRUE, linewidth = 1)
    )

  p2 <- plotly::plot_ly(scaling_adj, x = ~MDS1, y = ~MDS2, color = ~Group, colors = group_colors, showlegend = TRUE) %>%
    plotly::add_markers(text = ~Sample, hoverinfo = 'text') %>%
    plotly::layout(
      legend = l,
      xaxis = list(title = 'MDS 1', zeroline = FALSE, showticklabels = FALSE, range = xrange,
                   linecolor = '#cccccc', mirror = TRUE, linewidth = 1),
      yaxis = list(title = '', zeroline = FALSE, showticklabels = FALSE, range = yrange,
                   linecolor = '#cccccc', mirror = TRUE, linewidth = 1)
    )


  pl <- plotly::subplot(p1, p2, titleX = TRUE, titleY = TRUE) %>%
    plotly::layout(title = list(text = 'Sammon MDS plots', x = 0.08, y = 0.98),
                   xaxis = list(fixedrange=TRUE),
                   yaxis = list(fixedrange=TRUE),
                   margin = list(t = 60),
                   annotations = list(
                     list(x = 0.2 , y = 1.065, text = "Not Adjusted", showarrow = F, xref='paper', yref='paper'),
                     list(x = 0.8 , y = 1.065, text = "Adjusted", showarrow = F, xref='paper', yref='paper')),
                   shapes = list(
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
                     line = list(color = "#cccccc")))

  pl %>%
    plotly::config(displaylogo = FALSE,
                   displayModeBar = 'hover',
                   modeBarButtonsToRemove = c('toImage',
                                              'zoom2d',
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
    v <- limma::voomWithQualityWeights(y, mod, lib.size = lib.size, plot = TRUE)
    corfit <- limma::duplicateCorrelation(v, mod, block = pair)

    # second round
    v <- limma::voomWithQualityWeights(y, mod, lib.size = lib.size, block = pair, correlation = corfit$consensus.correlation, plot = TRUE)
    corfit <- limma::duplicateCorrelation(v, mod, block = pair)

    fit <- limma::lmFit(v, mod, correlation = corfit$consensus.correlation, block = pair)

  } else if (rna_seq) {
    # rna-seq not paired
    # get normalized lib size and voom
    v <- limma::voomWithQualityWeights(y, mod, lib.size = lib.size, plot = TRUE)
    fit  <- limma::lmFit(v, design = mod)

  } else if (length(pair) & !rna_seq) {
    # microarray paired
    corfit <- limma::duplicateCorrelation(y, mod, block = pair)
    fit <- limma::lmFit(y, mod, correlation = corfit$consensus.correlation, block = pair)

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
#' @seealso \code{\link{get_contrast_esets}}.
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
