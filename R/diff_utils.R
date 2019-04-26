#' Differential expression analysis of esets.
#'
#' After selecting control and test samples for a contrast, surrogate variable
#' analysis (\code{\link[sva]{sva}}) and differential expression analysis is performed.
#'
#' To specify a contrast, first select rows for control samples, type a group name
#' in the \emph{Control group name} text input box and click the \emph{Add Group}
#' button. Repeat for test samples. After control and test samples have been added
#' click the \emph{Done} button. If you make a mistake, click the \emph{Reset} button.
#'
#' For each GSE, analysis results are saved in the corresponding GSE
#' folder in \code{data_dir} that was created by \code{\link{get_raw}}. If analyses
#' needs to be repeated, previous results can be reloaded with \code{\link{load_diff}}
#' and supplied to the \code{prev_anals} parameter. In this case, previous
#' selections, names, and pairs will be reused.
#'
#'
#' @param esets List of annotated esets. Created by \code{\link{load_raw}}.
#' @param data_dir String specifying directory of GSE folders.
#' @param annot String, column name in fData common to all esets. For duplicated
#'   values in this column, the row with the highest interquartile range
#'   across selected samples will be kept. If meta-analysis will follow, appropriate
#'   values are "SYMBOL" (default - for gene level analysis) or, if all esets are
#'   from the same platform, "PROBE" (for probe level analysis).
#' @param prev_anals Previous result of \code{\link{diff_expr}}, which can
#'    be reloaded using \code{\link{load_diff}}. If present, previous
#'   selections, names, and pairs will be reused.
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
#' eset_path <- system.file('extdata', 'IBD', 'eset.rds', package='drugseqr')
#' eset <- readRDS(eset_path)
#'
#' anal <- diff_expr(eset)


diff_expr <- function (eset, data_dir = getwd(), annot = "SYMBOL", svanal = TRUE) {

  # check for annot column
  if (!annot %in% colnames(fData(eset)))
    stop(annot, " column in fData(eset) is missing.")

  # determine if this is rna seq data
  rna_seq <- 'norm.factors' %in% colnames(Biobase::pData(eset))

  # select/add contrast
  eset <- select_contrast(eset)

  # setup for differential expression
  setup <- diff_setup(eset, svanal, rna_seq)

  # remove rows with duplicated/NA annot (SYMBOL or ENTREZID)
  dups <- tryCatch ({
    iqr_replicates(eset, setup$mod, setup$svobj, annot)

    }, error = function(err) {err$message <- "Couldn't fit model."; stop(err)})

  # differential expression
  anal <- diff_anal(dups$eset, dups$exprs_sva, cons$contrasts, cons$levels,
                    setup$modsv, gse_dir, gse_name, annot, rna_seq)

  return (anal)
}



# ------------------------


#' Generate model matrix with surrogate variables.
#'
#' Used by \code{diff_expr} to create model matrix with surrogate variables
#' in order to run \code{diff_anal}.
#'
#' @param eset Annotated eset with samples selected during \code{add_contrasts}.
#' @param group_levels Character vector of unique group names created by
#'    \code{add_contrasts}.
#'
#' @seealso \code{\link{add_contrasts}}, \code{\link{diff_expr}}.
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
        message(gse_name, ": sva failed - continuing without.")
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



# Removes features with replicated annotation.
#
# For rows with duplicated annot, highested IQR retained.
#
# @inheritParams diff_expr
# @inheritParams diff_setup
# @param mod Model matrix without surrogate variables. generated by \code{diff_setup}.
# @param svobj Result from \code{sva} function called during \code{diff_setup}.
# @param annot feature to use to remove replicates.
# @param rm.dup remove duplicates (same measure, multiple ids)?
#
# @return List with:
#    \item{eset}{Expression set with unique features at probe or gene level.}
#    \item{exprs_sva}{Expression data from eset with effect of surrogate
#       variable removed.}

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


# Run limma analysis.
#
# Runs limma differential expression analysis on all contrasts selected by
# \code{add_contrasts}. Analysis performed with and without surrogate
# variables discovered by \code{diff_setup}. Also prints MDS plot and saves
# results.
#
# @param eset Annotated eset created by \code{load_raw}. Replicate features and
#   non-selected samples removed by \code{iqr_replicates}.
# @param exprs_sva Expression data with surrogate variables removed. Created by
#    \code{iqr_replicates}
# @param contrasts Character vector generated by \code{add_contrasts}.
# @param group_levels Character vector of group names generated by
#    \code{add_contrasts}.
# @param mod, modsv Model matrix generated by \code{diff_setup}. With
#   and without surrogate variables.
# @param svobj Result from \code{sva} function called during \code{diff_setup}.
# @param gse_dir String, path to directory with GSE folders.
# @param gse_name String, name of GSE.
# @param annot String, either "ENTREZID" or "SYMBOL" for probe or gene level
#   analysis respectively. If "ENTREZID", appends "_entrezid.rds" to save name.
#
# @seealso \code{\link{diff_expr}}.
# @return List, final result of \code{diff_expr}. Used for subsequent
#   meta-analysis.


diff_anal <- function(eset, exprs_sva, contrasts, group_levels,
                      modsv, gse_dir, gse_name, annot = "SYMBOL", rna_seq = FALSE){


  # differential expression (surrogate variables modeled and not)
  ebayes_sv <- fit_ebayes(eset, contrasts, modsv, rna_seq)

  # annotate/store results
  top_tables <- list()
  contrast_names <- paste(gse_name, contrasts, sep = "_")

  for (i in seq_along(contrast_names)) {
    top_genes <- limma::topTable(ebayes_sv, coef = i, n = Inf)
    num_sig <- sum(top_genes$adj.P.Val < 0.05)
    top_tables[[contrast_names[i]]] <- top_genes
  }

  # only store phenoData (exprs and fData large)
  pdata <- pData(eset)

  # save to disk
  diff_expr <- list(pdata = pdata, top_tables = top_tables, ebayes_sv = ebayes_sv, annot = annot)
  save_name <- paste(gse_name, "diff_expr", tolower(annot), sep = "_")
  save_name <- paste0(save_name, ".rds")

  saveRDS(diff_expr, file = paste(gse_dir, save_name, sep = "/"))
  return (diff_expr)
}


# ------------------------


# Perform eBayes analysis from limma.
#
# Generates contrast matrix then runs eBayes analysis from limma.
#
# @param eset Annotated eset created by \code{load_raw}. Non-selected samples
#    and duplicate features removed by \code{add_contrasts} and
#    \code{iqr_replicates}.
# @param contrasts Character vector of contrast names generated by
#    \code{add_contrasts}.
# @param mod Model matrix generated by \code{diff_setup}. With
#   or without surrogate variables.
#
# @return result from call to limma \code{eBayes}.

fit_ebayes <- function(eset, contrasts, mod, rna_seq = FALSE) {

  if (rna_seq) {
    # get normalize lib size and voom
    lib.size <- pData(eset)$lib.size * pData(eset)$norm.factors
    voom <- limma::voom(exprs(eset), mod, lib.size)
    fit  <- limma::lmFit(voom)

  } else {
    fit <- limma::lmFit(exprs(eset), mod)
  }

  contrast_matrix <- limma::makeContrasts(contrasts = contrasts, levels = mod)
  fit <- limma::contrasts.fit(fit, contrast_matrix)
  return (limma::eBayes(fit))
}



# ------------------------


# Adjusts expression data for surrogate variables.
#
# Factors out effect of surrogate variables discovered during surrogate variable
# analysis.
#
# @param y Expression data of eset.
# @param mod Full model matrix supplied to \code{sva}.
# @param svs Surrogate variables returned by \code{sva} (svobj$sv).
#
# @seealso \code{\link{get_contrast_esets}}.
# @return Expression data with effects of svs removed.

clean_y <- function(y, mod, svs) {

  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}


# ------------------------



#' Loads previous differential expression analyses.
#'
#' @param gse_names Character vector specifying GSE names to be loaded.
#' @param data_dir String specifying directory of GSE folders.
#' @param annot Level of previous analysis (e.g. "SYMBOL" or "PROBE").
#'
#' @export
#' @return Result of previous call to \code{\link{diff_expr}}.
#' @examples
#' library(lydata)
#'
#' data_dir <- system.file("extdata", package = "lydata")
#' gse_names <- c("GSE9601", "GSE34817")
#' prev <- load_diff(gse_names, data_dir)

load_diff <- function(gse_name, data_dir = getwd(), annot = "SYMBOL") {

  anal    <- list()
  gse_dir <- file.path(data_dir, gse_name)

  pattern <- paste0(gse_name, "_diff_expr_", tolower(annot), '.rds')

  # get paths
  path <- list.files(gse_dir, pattern, full.names = TRUE)
  if (length(path))
    anal[[gse_name]] <- readRDS(path)

  return (anal)
}

