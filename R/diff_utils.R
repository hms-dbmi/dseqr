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
#' data_dir <- system.file('extdata', 'IBD', package='drugseqr')
#' eset <- load_seq(data_dir)
#' anal <- diff_expr(eset, data_dir)

diff_expr <- function (eset, data_dir, annot = "SYMBOL", svanal = TRUE, prev_anal = NULL) {

  # check for annot column
  if (!annot %in% colnames(fData(eset)))
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
  anal <- diff_anal(dups$eset, dups$exprs_sva, setup$modsv, data_dir, annot, rna_seq)
  return (anal)
}



# ------------------------


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
  group_levels = c('control', 'test')
  group <- factor(pData(eset)$group, levels = group_levels)

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
    expr <- unique(data.table::data.table(exprs(eset), PROBE))[, PROBE := NULL]

    # sva or svaseq
    sva_fun <-ifelse(rna_seq, sva::svaseq, sva::sva)

    svobj <- tryCatch (
      {utils::capture.output(svobj <- sva_fun(as.matrix(expr), mod, mod0)); svobj},

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


# ------------------------



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
    exprs_sva <- clean_y(exprs(eset), mod, svobj$sv)
  } else {
    exprs_sva <- exprs(eset)
  }

  # add inter-quartile ranges, row, and feature data to exprs data
  data <- as.data.frame(exprs_sva)
  data$iqrange <- matrixStats::rowIQRs(exprs_sva)
  data$row <- 1:nrow(data)
  data[, colnames(fData(eset))] <- fData(eset)

  # remove rows with NA annot (occurs if annot is SYMBOL)
  data <- data[!is.na(data[, annot]), ]

  # for rows with same annot, keep highest IQR
  data <- data.table::data.table(data)
  data <- data[data[, .I[which.max(iqrange)], by = eval(annot)]$V1]

  # use row number to keep selected features
  eset <- eset[data$row, ]
  exprs_sva <- exprs_sva[data$row, ]

  # use annot for feature names
  featureNames(eset) <- fData(eset)[, annot]
  row.names(exprs_sva)   <- fData(eset)[, annot]

  if (rm.dup) {
    not.dup <- !duplicated(exprs_sva)
    eset <- eset[not.dup, ]
    exprs_sva <- exprs_sva[not.dup, ]
  }

  return (list(eset = eset, exprs_sva = exprs_sva))
}


# ------------------------


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


diff_anal <- function(eset, exprs_sva, modsv, data_dir, annot = "SYMBOL", rna_seq = TRUE){

  group_levels <- c('control', 'test')
  contrasts <- 'test-control'

  # differential expression (surrogate variables modeled and not)
  ebayes_sv <- fit_ebayes(eset, contrasts, modsv, rna_seq)

  # get results
  top_table <- limma::topTable(ebayes_sv, coef = 1, n = Inf)
  num_sig <- sum(top_table$adj.P.Val < 0.05)
  cat(contrasts, "(# p < 0.05):", num_sig, "\n")

  # voom (normalize/log) on exprs_sva
  if (rna_seq) {
    lib.size <- pData(eset)$lib.size * pData(eset)$norm.factors
    exprs_sva <- exprs_sva[edgeR::filterByExpr(exprs_sva), ]
    exprs_sva <- limma::voom(exprs_sva, modsv, lib.size)$E
  }

  # setup plot items
  group <- factor(pData(eset)$group, levels = group_levels)
  palette <- RColorBrewer::brewer.pal(12, "Paired")
  colours <- palette[group]

  # Add extra space to right of plot area
  graphics::par(mai = c(1, 1, 1, 1.4))

  # plot MDS
  limma::plotMDS(exprs_sva, pch = 19, main='MDS plot (after sva)', col = colours)
  graphics::legend("topright", inset = c(-0.18, 0), legend = group_levels,
                   fill = unique(colours), xpd = TRUE, bty = "n", cex = 0.65)

  # only store phenoData (exprs and fData large)
  pdata <- Biobase::pData(eset)

  # save to disk
  diff_expr <- list(pdata = pdata, top_table = top_table, ebayes_sv = ebayes_sv, annot = annot)
  save_name <- paste("diff_expr", tolower(annot), sep = "_")
  save_name <- paste0(save_name, ".rds")

  saveRDS(diff_expr, file.path(data_dir, save_name))
  return(diff_expr)
}


# ------------------------


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
#'
#' @return result from call to limma \code{eBayes}.

fit_ebayes <- function(eset, contrasts, mod, rna_seq = TRUE) {

  if (rna_seq) {
    # get normalized lib size and voom
    lib.size <- pData(eset)$lib.size * pData(eset)$norm.factors
    v <- limma::voom(exprs(eset), mod, lib.size)
    fit  <- limma::lmFit(v)

  } else {
    fit <- limma::lmFit(exprs(eset), mod)
  }

  contrast_matrix <- limma::makeContrasts(contrasts = contrasts, levels = mod)
  fit <- limma::contrasts.fit(fit, contrast_matrix)
  return (limma::eBayes(fit))
}



# ------------------------


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

