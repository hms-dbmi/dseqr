
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
#' anal <- diff_expr(eset, data_dir)
diff_expr <- function (eset, data_dir, anal_name, annot = "SYMBOL", svanal = TRUE, prev_anal = NULL) {

  # check for annot column
  if (!annot %in% colnames(Biobase::fData(eset)))
    stop(annot, " column in fData(eset) is missing.")

  # determine if this is rna seq data
  rna_seq <- 'norm.factors' %in% colnames(Biobase::pData(eset))

  # are we re-using selections from previous analysis?
  if (!is.null(prev_anal)) {
    # retain selected only and add group to pdata
    eset <- eset[, row.names(prev_anal$pdata)]
    Biobase::pData(eset)$group <- prev_anal$pdata$group

  } else {
    # select/add contrast
    eset <- select_contrast(eset)
  }

  # setup for differential expression
  setup <- diff_setup(eset, svanal, rna_seq)

  # remove rows with duplicated/NA annot (SYMBOL or ENTREZID)
  dups <- tryCatch ({
    iqr_replicates(eset, setup$mod, setup$svobj, annot)

  }, error = function(err) {err$message <- "Couldn't fit model."; stop(err)})

  # differential expression

  anal <- diff_anal(eset = dups$eset,
                    anal_name = anal_name,
                    exprs_sva = dups$exprs_sva,
                    modsv = setup$modsv,
                    data_dir = data_dir,
                    annot = annot,
                    rna_seq = rna_seq)
  return (anal)
}




#' Generate model matrix with surrogate variables.
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
diff_setup <- function(eset, svanal = TRUE, rna_seq = TRUE){

  # incase svanal FALSE
  svobj <- list("sv" = NULL)

  # make full and null model matrix
  group_levels = c('ctrl', 'test')
  group <- factor(Biobase::pData(eset)$group, levels = group_levels)

  mod <- stats::model.matrix(~0 + group)
  mod0 <- stats::model.matrix(~1, data = group)

  colnames(mod)[1:length(group_levels)] <- group_levels


  if (svanal) {
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

    svobj <- tryCatch (
      {utils::capture.output(svobj <- sva_fun(expr, mod, mod0)); svobj},

      error = function(cond) {
        message("sva failed - continuing without.")
        return(list("sv" = NULL))
      })
  }

  if (is.null(svobj$sv) || svobj$n.sv ==  0) {
    svobj$sv <- NULL
    modsv <- mod
  } else {
    modsv <- cbind(mod, svobj$sv)
    colnames(modsv) <- c(colnames(mod), paste("SV", 1:svobj$n.sv, sep = ""))
  }
  return (list("mod" = mod, "modsv" = modsv, "svobj" = svobj))
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
#' @return List with:
#'    \item{eset}{Expression set with unique features at probe or gene level.}
#'    \item{exprs_sva}{Expression data from eset with effect of surrogate
#'       variable removed.}
iqr_replicates <- function(eset, mod = NULL, svobj = NULL, annot = "SYMBOL", rm.dup = FALSE) {

  # for R CMD check
  iqrange = SYMBOL = NULL

  if (length(svobj) > 0) {
    # get eset with surrogate variables modeled out
    exprs_sva <- clean_y(Biobase::exprs(eset), mod, svobj$sv)
  } else {
    exprs_sva <- Biobase::exprs(eset)
  }

  # add inter-quartile ranges, row, and feature data to exprs data
  data <- as.data.frame(exprs_sva)
  data$iqrange <- matrixStats::rowIQRs(exprs_sva)
  data$row <- 1:nrow(data)
  data[, colnames(Biobase::fData(eset))] <- Biobase::fData(eset)

  # remove rows with NA annot (occurs if annot is SYMBOL)
  data <- data[!is.na(data[, annot]), ]

  # for rows with same annot, keep highest IQR
  data <- data.table::data.table(data)
  data <- data[data[, .I[which.max(iqrange)], by = eval(annot)]$V1]

  # use row number to keep selected features
  eset <- eset[data$row, ]
  exprs_sva <- exprs_sva[data$row, ]

  # use annot for feature names
  Biobase::featureNames(eset) <- Biobase::fData(eset)[, annot]
  row.names(exprs_sva)   <- Biobase::fData(eset)[, annot]

  if (rm.dup) {
    not.dup <- !duplicated(exprs_sva)
    eset <- eset[not.dup, ]
    exprs_sva <- exprs_sva[not.dup, ]
  }

  return (list(eset = eset, exprs_sva = exprs_sva))
}


format_scaling <- function(scaling, with_sva, group, exprs) {

  scaling %>%
    dplyr::rename('MDS1' = V1, 'MDS2' = V2) %>%
    dplyr::mutate(Sample = row.names(exprs)) %>%
    dplyr::mutate(Group = factor(group)) %>%
    dplyr::mutate(Group =  dplyr::recode(Group, ctrl = 'Control', test = 'Test')) %>%
    dplyr::mutate(Title = ifelse(with_sva, 'with sva', 'without sva'))
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
#' @param exprs_sva Expression data with surrogate variables removed. Created by
#'    \code{iqr_replicates}
#' @param modsv Model matrix generated by \code{diff_setup}. With
#'   and without surrogate variables.
#' @param data_dir String, path to directory with raw data.
#' @param annot String, either "ENTREZID" or "SYMBOL" for probe or gene level
#'   analysis respectively. If "ENTREZID", appends "_entrezid.rds" to save name.
#' @param rna_seq Is the analysis for RNA-seq? Auto inferred in \code{diff_expr}.
#'
#' @seealso \code{\link{diff_expr}}.
#' @return List, final result of \code{diff_expr}. Used for subsequent
#'   meta-analysis.
diff_anal <- function(eset, anal_name, exprs_sva, modsv, data_dir, annot = "SYMBOL", rna_seq = TRUE){

  group_levels <- c('ctrl', 'test')
  contrasts <- 'test-ctrl'

  # differential expression (surrogate variables modeled and not)
  ebayes_sv <- fit_ebayes(eset, contrasts, modsv, rna_seq)

  # get results
  top_table <- limma::topTable(ebayes_sv, coef = 1, n = Inf)

  # voom (normalize/log) on exprs_sva for MDS plot
  if (rna_seq) {
    lib.size <- Biobase::pData(eset)$lib.size * Biobase::pData(eset)$norm.factors
    exprs_sva <- exprs_sva[edgeR::filterByExpr(exprs_sva), ]
    exprs_sva <- limma::voom(exprs_sva, modsv, lib.size)$E
  }

  # only store phenoData (exprs and fData large)
  pdata <- Biobase::pData(eset)

  # for MDS plot
  mds <- get_mds(Biobase::exprs(eset), exprs_sva, pdata$group)

  # save to disk
  diff_expr <- list(pdata = pdata, top_table = top_table, ebayes_sv = ebayes_sv, annot = annot, mds = mds)
  save_name <- paste("diff_expr", tolower(annot), anal_name, sep = "_")
  save_name <- paste0(save_name, ".rds")

  saveRDS(diff_expr, file.path(data_dir, save_name))
  return(diff_expr)
}

#' Get scalings for MDS plots
#'
#' For interactive MDS plot of expression values with and without surrogate variable analysis.
#'
#' @param exprs \code{matrix} of expression values.
#' @param exprs_sva \code{matrix} of expression values with surrogate variables regressed out.
#' @param group Character vector with values \code{'control'} and \code{'test'} indicating group membership.
#' @importFrom magrittr "%>%"
#'
#' @return List of tibbles with MDS scalings with and without SVA
#' @export
get_mds <- function(exprs, exprs_sva, group) {

  suggests <- c('factoextra', 'MASS', 'plotly')
  if (!is.installed(suggests, level = 'message')) return(NULL)

  # get_dist acts on rows
  exprs <- t(exprs[complete.cases(exprs), ])
  exprs_sva <- t(exprs_sva[complete.cases(exprs_sva), ])

  dist <- factoextra::get_dist(exprs, method = 'spearman')
  dist_sva <- factoextra::get_dist(exprs_sva, method = 'spearman')

  # sammon scaling for dimensionality reduction
  capture.output({
    scaling <- tibble::as_tibble(MASS::sammon(dist, k = 2)$points) %>%
      format_scaling(with_sva = FALSE, group = group, exprs = exprs)

    scaling_sva <- tibble::as_tibble(MASS::sammon(dist_sva, k = 2)$points) %>%
      format_scaling(with_sva = TRUE, group = group, exprs = exprs)
  })

  return(list(scaling = scaling, scaling_sva = scaling_sva))
}

#' Plot MDS plots
#'
#' @param scaling tibble with columns MDS1 and MDS2 corresponding to differential expression without SVA
#' @param scaling_sva tibble with columns MDS1 and MDS2 corresponding to differential expression with SVA
#'
#' @return plotly object
#' @export
plotMDS <- function(scaling, scaling_sva) {

  # make x and y same range
  xrange <- range(c(scaling$MDS1, scaling_sva$MDS1))
  yrange <- range(c(scaling$MDS2, scaling_sva$MDS2))

  addx <- diff(xrange) * .10
  addy <- diff(yrange) * .10

  xrange <- xrange + c(-addx, addx)
  yrange <- yrange + c(-addy, addy)

  p1 <- plotly::plot_ly(scaling, x = ~MDS1, y = ~MDS2, color = ~Group, colors = c('#337ab7', '#e6194b'), showlegend = FALSE) %>%
    plotly::config(displayModeBar = FALSE) %>%
    plotly::add_markers(text = ~Sample, hoverinfo = 'text') %>%
    plotly::layout(
      xaxis = list(title = 'MDS 1', zeroline = FALSE, showticklabels = FALSE, range = xrange),
      yaxis = list(title = 'MDS 2', zeroline = FALSE, showticklabels = FALSE, range = yrange)
    )

  p2 <- plotly::plot_ly(scaling_sva, x = ~MDS1, y = ~MDS2, color = ~Group, colors = c('#337ab7', '#e6194b'), showlegend = FALSE) %>%
    plotly::config(displayModeBar = FALSE) %>%
    plotly::add_markers(text = ~Sample, hoverinfo = 'text') %>%
    plotly::layout(
      xaxis = list(title = 'MDS 1', zeroline = FALSE, showticklabels = FALSE, range = xrange),
      yaxis = list(title = 'MDS 2', zeroline = FALSE, showticklabels = FALSE, range = yrange)
    )

  plotly::subplot(p1, p2, shareY = TRUE, shareX = TRUE, titleX = TRUE, titleY = TRUE) %>%
    plotly::layout(title = list(text = 'Sammon MDS plots', x = 0.08, y = 0.98),
                   margin = list(t = 60),
                   annotations = list(
                     list(x = 0.2 , y = 1.045, text = "Without SVA", showarrow = F, xref='paper', yref='paper'),
                     list(x = 0.8 , y = 1.045, text = "With SVA", showarrow = F, xref='paper', yref='paper')),
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
                     fillcolor = toRGB("gray80"),
                     line = list(color = "transparent")))

  return(pl)
}




#' Perform eBayes analysis from limma.
#'
#' Generates contrast matrix then runs eBayes analysis from limma.
#'
#' @param eset Annotated eset created by \code{load_raw}. Non-selected samples
#'    and duplicate features removed by \code{add_contrasts} and
#'    \code{iqr_replicates}.
#' @param contrasts Character vector of contrast names generated by
#'    \code{add_contrasts}.
#' @param mod Model matrix generated by \code{diff_setup}. With
#'   or without surrogate variables.
#' @param rna_seq Is this an RNA-seq \code{eset}? Default is \code{TRUE}.
#'
#' @return result from call to limma \code{eBayes}.

fit_ebayes <- function(eset, contrasts, mod, rna_seq = TRUE) {

  if (rna_seq) {
    # get normalized lib size and voom
    lib.size <- Biobase::pData(eset)$lib.size * Biobase::pData(eset)$norm.factors
    v <- limma::voomWithQualityWeights(Biobase::exprs(eset), design = mod, lib.size = lib.size, plot = TRUE)
    fit  <- limma::lmFit(v, design = mod)

  } else {
    fit <- limma::lmFit(Biobase::exprs(eset), mod)
  }

  contrast_matrix <- limma::makeContrasts(contrasts = contrasts, levels = mod)
  fit <- limma::contrasts.fit(fit, contrast_matrix)
  return (limma::eBayes(fit, robust = TRUE))
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
clean_y <- function(y, mod, svs) {

  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}
