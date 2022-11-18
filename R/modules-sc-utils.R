#' Get predicted annotation for label transfer
#'
#' Clusters with an average prediction scores below \code{min.score} retain their original labels.
#'
#' @param ref_preds data.frame generated in \code{labelTransferForm} on event \code{submit_transfer}
#' @param ref_name Name of reference analysis that labels are transfered from.
#' @param dataset_name Name of analysis that labels are transfered to.
#' @param sc_dir Directory containing folders with analyses for \code{ref_name} and \code{dataset_name}.

#' @return Character vector of predicted labels from \code{ref_name}.
#'
#' @keywords internal
get_pred_annot <- function(ref_preds, ref_name, dataset_name, sc_dir) {

  senv <- loadNamespace('SingleR')

  if (ref_name %in% ls(senv)) {
    pred_annot <- pretty.unique(ref_preds)

  } else if (ref_name == 'reset') {
    # reset annotation
    # load query annotation
    query_annot_path <- scseq_part_path(sc_dir, dataset_name, 'annot')
    query_annot <- qs::qread(query_annot_path)
    pred_annot <- as.character(seq_along(query_annot))

  } else {
    ref_annot_path <- scseq_part_path(sc_dir, ref_name, 'annot')
    ref_annot <- qs::qread(ref_annot_path)
    ref_annot <- remove.unique(ref_annot)

    ref_preds <- ref_annot[as.numeric(ref_preds)]
    pred_annot <- pretty.unique(ref_preds)
  }
  return(pred_annot)
}

remove.unique <- function(annot) {
  annot <- gsub('_\\d+$', '', annot)
  annot <- gsub(' \\(\\d+\\)$', '', annot)
  return(annot)
}

pretty.unique <- function(annot) {
  annot <- make.unique(annot, '_')

  is.dup <- grepl('_\\d+$', annot)
  dup.num <- gsub('^.+?_(\\d+)$', '\\1', annot[is.dup])
  dup.num <- as.integer(dup.num)
  annot <- gsub('_\\d+$', '', annot)
  annot[is.dup] <- paste0(annot[is.dup], ' (', dup.num+1, ')')
  return(annot)
}

#' Check that label predictions aren't from overwritten datasets
#'
#' @param preds List of predictions
#' @param sc_dir Directory with single-cell datasets
#'
#' @return \code{preds} that were generated from current datasets
#'
#' @keywords internal
validate_preds <- function(preds, sc_dir) {
  senv <- loadNamespace('celldex')
  external <- ls(senv)

  ref_names <- names(preds)
  dated <- c()

  dated <- sapply(ref_names, function(ref_name) {
    if (ref_name %in% external) return(FALSE)
    resoln_name <- get_resoln_name(sc_dir, ref_name)

    ref_path <- scseq_part_path(sc_dir, resoln_name, 'scseq_sample')
    ref_date <- file.info(ref_path)$ctime
    ref_date <- as.numeric(ref_date)
    res_date <- attr(preds[[ref_name]], 'ref_date')

    dated <- !identical(res_date, ref_date)

    return(dated)
  })

  return(preds[!dated])
}

#' Run differential abundance analysis
#'
#' @param obj \code{SingleCellExperiment} or count matrix
#' @param orig.ident factor of group identities with length \code{ncol(obj)}
#' @return \code{data.frame} with differential abundance results calculated by
#' \link[limma]{topTable}
#'
#' @keywords internal
diff_abundance <- function(obj, annot = NULL, pairs = NULL, orig.ident = NULL, filter = TRUE) {

  if (methods::is(obj, 'SingleCellExperiment')) {
    abundances <- table(obj$cluster, obj$batch)
    abundances <- unclass(abundances)
    row.names(abundances) <- annot
    genes <- data.frame(cluster = levels(obj$cluster))

    extra.info <- obj@colData[match(colnames(abundances), obj$batch),]

  } else {
    genes <- data.frame(row.names = row.names(obj))
    abundances <- obj
    extra.info <- data.frame(orig.ident)
  }

  y.ab <- edgeR::DGEList(abundances, samples=extra.info, genes = genes)

  # adjusted counts for calculating logFC if two samples
  adj <- edgeR::equalizeLibSizes(y.ab)$pseudo.counts
  adj <- round(adj)

  group <- y.ab$samples$orig.ident
  group <- stats::relevel(group, 'ctrl')
  y.ab$samples$group <- group

  # if only two samples, sort by logFC and return
  if (ncol(adj) == 2) {
    ord <- match(c('test', 'ctrl'), y.ab$samples$group)
    adj <- as.data.frame(adj[, ord])
    adj$logFC <- log2(adj[[1]]+1) - log2(adj[[2]]+1)
    adj <- adj[order(abs(adj$logFC), decreasing = TRUE), 'logFC', drop=FALSE]
    return(adj)
  }

  if (!is.null(pairs))
    y.ab$samples$pair <- factor(pairs[colnames(y.ab),])

  if (filter) {
    keep <- edgeR::filterByExpr(y.ab, group=y.ab$samples$orig.ident)
    y.ab <- y.ab[keep,]
  }

  # setup eset so that it will work with run_limma
  eset <- Biobase::ExpressionSet(y.ab$counts,
                                 Biobase::AnnotatedDataFrame(y.ab$samples),
                                 Biobase::AnnotatedDataFrame(y.ab$genes))

  Biobase::assayDataElement(eset, 'vsd') <- Biobase::exprs(eset)
  Biobase::fData(eset)[, c('SYMBOL', 'ENTREZID')] <- row.names(eset)

  prev <- list(pdata = Biobase::pData(eset))
  eset <- crossmeta::run_limma_setup(eset, prev)
  lm_fit <- crossmeta::run_limma(eset, filter = filter)

  tt <- crossmeta::get_top_table(lm_fit, with.es = FALSE)
  tt$ENTREZID <- NULL

  # add cluster number
  tt$cluster <- Biobase::fData(eset[row.names(tt), ])$cluster

  return(tt)
}



#' Get label transfer choices data.frame
#'
#' Used for Single Cell tab label transfer selectizeInput
#'
#' @param anal_options Names list of analysis options with names \code{'Individual'} and \code{'Integrated'}.
#' @param preds Named list of predicted cluster labels. Names are values in \code{anal_options} lists.
#'
#' @return data.frame with columns \code{value}, \code{label}, \code{type}, and \code{preds}.
#'
#' @keywords internal
get_label_transfer_choices <- function(anal_options, selected_anal, preds) {

  # omit previous and current
  is.recent <- anal_options$type == 'Previous Session'
  is.sel <- anal_options$name == selected_anal
  type   <- utils::tail(anal_options$type[is.sel], 1)
  type   <-  ifelse(type %in% c('Integrated', 'Individual'), selected_anal, type)

  anal_options <- anal_options[!is.sel & !is.recent, ]


  external <- c('BlueprintEncodeData',
                'DatabaseImmuneCellExpressionData',
                'HumanPrimaryCellAtlasData',
                'MonacoImmuneData',
                'MouseRNAseqData',
                'NovershternHematopoieticData')

  external_type <- rep('External Reference', length(external))
  external_is.int <- rep(TRUE, length(external))

  choices <- data.frame(
    label = c('Reset Labels', external, anal_options$name),
    value = c('reset', gsub(' ', '', external), anal_options$name),
    optionLabel = c('Reset Labels', external, anal_options$optionLabel),
    type = factor(c(type, external_type, anal_options$type), ordered = TRUE),
    is.int = c(TRUE, external_is.int, anal_options$is.int),
    stringsAsFactors = FALSE
  )

  choices$preds <- choices$value %in% names(preds)
  choices <- rbind(NA, choices)

  return(choices)
}


# used to set selected annotation reference safely when switch datasets
get_selected_from_transfer_name <- function(transfer_name, dataset_name) {
  if (is.null(transfer_name)) return(NULL)

  transfer_datasets <- strsplit(transfer_name, ' \U2192 ')[[1]]
  ref_name <- transfer_datasets[1]
  query_name <- transfer_datasets[2]

  if (query_name != dataset_name) return(NULL)

  return(ref_name)
}




#' Utility to generate filename for single cell download csv
#'
#' @param cluster Character vector of cluster names
#' @param anal result of diff_expr_scseq
#' @param comparison_type either 'samples' or 'clusters'
#'
#' @return File name string
#'
#' @keywords internal
sc_dl_filename <- function(cluster, anal, comparison_type) {

  # remove underscores and spaces for cluster(s) and analysis name
  cluster <- gsub('[_ ]', '-', cluster)
  anal <- gsub('[_ ]', '-', anal)

  if (comparison_type == 'samples')
    cluster <- paste0('test-vs-ctrl_', paste(cluster, collapse = '-'))


  paste('single-cell', anal, cluster, paste0(Sys.Date(), '.csv'), sep='_')
}


#' Get cluster choices data.frame for selectize dropdown
#'
#' @param clusters Character vector of cluster names
#' @param dataset_dir Directory with single cell dataset.
#' @param ... Named arguments to \code{\link{get_cluster_stats}}.
#' @param sample_comparison is this for test vs control comparion? Default is \code{FALSE}.
#'
#' @return data.frame with columns for rendering selectizeInput cluster choices
#'
#' @keywords internal
get_cluster_choices <- function(clusters, sample_comparison = FALSE, with_all = FALSE, ...) {

  testColor <- get_palette(clusters, with_all=with_all)
  testColor <- testColor[seq_along(clusters)]

  inds <- seq_along(clusters)
  not.inds <- clusters != inds & clusters != 'All Clusters'
  clusters[not.inds] <- paste0(inds[not.inds], ': ', clusters[not.inds])


  # cluster choices are the clusters themselves
  # value is original cluster number
  choices <- data.frame(name = stringr::str_trunc(clusters, 27),
                        value = seq_along(clusters),
                        label = clusters,
                        testColor,
                        row.names = NULL, stringsAsFactors = FALSE)

  cluster_stats <- get_cluster_stats(sample_comparison = sample_comparison, ...)

  if (sample_comparison) {
    # non-formatted stats for item/hover
    choices$ntest <- cluster_stats$ntest
    choices$nctrl <- cluster_stats$nctrl
    choices$nsig  <- cluster_stats$nsig
    choices$nbig  <- cluster_stats$nbig
    choices$ntest_each <- cluster_stats$ntest_each
    choices$nctrl_each <- cluster_stats$nctrl_each

    # formatted for options
    choices$nctrlf <- html_space(cluster_stats$nctrl)
    choices$nsigf  <- html_space(cluster_stats$nsig)
    choices$nbigf  <- html_space(cluster_stats$nbig)

    choices <- rbind(utils::tail(choices, 1), utils::head(choices, nrow(choices)-1))

  } else {
    # cluster stats has 1 too many for cluster comparison (for 'All Clusters')
    idx <- seq_len(nrow(choices))

    # show the cell numbers/percentages
    choices$ncells <- cluster_stats$ncells[idx]
    choices$pcells <- html_space(round(cluster_stats$pcells[idx]))
    choices$ncellsf <- html_space(cluster_stats$ncells[idx])
  }

  return(choices)
}

#' Check if folder has single-cell files
#'
#' @param dataset_names Character vector of folder names to check
#' @param sc_dir Path to directory with \code{dataset_names} folders
#'
#' @return Boolean vector with \code{length(dataset_names)} indicating folders that
#' have single-cell files
#' @export
#'
check_has_scseq <- function(dataset_names, sc_dir) {
  if (!length(dataset_names)) return(FALSE)

  sapply(dataset_names, function(dataset_name) {
    fnames <- list.files(file.path(sc_dir, dataset_name))
    any(fnames %in% c('scseq.qs', 'shell.qs'))
  })
}


#' Get data.frame of single-cell dataset choices
#'
#' @param sc_dir Directory containing single-cell datasets
#' @param prev name of previously selected dataset
#'
#' @return data.frame of single-cell dataset choices for selectizeInput
#'
#' @keywords internal
get_sc_dataset_choices <- function(sc_dir, prev = NULL) {


  # exclude missing from integrated (e.g. manual delete)
  integrated <- get_integrated_datasets(sc_dir)
  has.scseq <- check_has_scseq(integrated, sc_dir)
  integrated <- integrated[has.scseq]

  int_type <- get_scdata_type(integrated, sc_dir, none = 'Integrated')
  sub <- duplicated(int_type) | duplicated(int_type, fromLast = TRUE)
  int_type[!sub] <- 'Integrated'

  int_opt <- integrated
  int_opt[sub] <- stringr::str_replace(int_opt[sub], paste0(int_type[sub], '_'), '')

  dataset_names <- list.dirs(sc_dir, full.names = FALSE, recursive = FALSE)
  individual <- setdiff(dataset_names, integrated)

  # exclude individual without scseq (e.g. folder with fastq.gz files only)
  has.scseq <- check_has_scseq(individual, sc_dir)
  individual <- individual[unlist(has.scseq)]

  # founder for subsets as type
  ind_type <- get_scdata_type(individual, sc_dir, none = 'Individual')

  # exclude founder name from option label
  sub <- ind_type != 'Individual'
  ind_opt <- individual
  ind_opt[sub] <- stringr::str_replace(ind_opt[sub], paste0(ind_type[sub], '_?'), '')
  ind_opt[ind_opt == ""] <- individual[ind_opt == ""]


  label <- c(integrated, individual)
  type <- c(int_type, ind_type)
  opt <- c(int_opt, ind_opt)
  is.int <- c(rep(TRUE, length(integrated)),
              rep(FALSE, length(individual)))

  # set previously selected
  if (length(label)) {
    if (is.null(prev) || !prev %in% label) prev <- label[1]

    prev_type <- type[label == prev]
    founder <- get_founder(sc_dir, prev)
    if (prev_type == founder)
      prev <- label[type %in% founder]

    type <- c(rep('Previous Session', length(prev)), type)
    is.int_prev <- is.int[label %in% prev]
    is.int <- c(is.int_prev, is.int)
    label <- c(prev, label)
    opt <- c(prev, opt)
  }

  choices <- data.frame(value = seq_along(label),
                        name = label,
                        label = label,
                        type = type,
                        is.int = is.int,
                        optionLabel = opt,
                        stringsAsFactors = FALSE)

  idx <- index_last(type, last = c('Integrated', 'Individual'))
  choices <- choices[idx, ]
  return(choices)
}

index_last <- function(values, last) {
  is.last <- numeric()
  for (val in last)
    is.last <- c(is.last, which(values == val))


  idx <- seq_along(values)
  not.last <- setdiff(idx, is.last)
  c(not.last, is.last)
}

#' Convert single-cell datasets data.frame to list
#'
#' Used to allow client side updates of single-cell dataset choices. This
#' is needed to prevent de-selection of current dataset after subset, intgration,
#' and loading.
#'
#' @param datasets data.frame of datasets from \link{get_sc_dataset_choices}
#'
#' @return list
#' @keywords internal
#'
datasets_to_list <- function(datasets) {
  res <- datasets$value
  names(res) <- datasets$optionLabel
  types <- datasets$type
  res <- lapply(unique(types), function(type) res[types==type])
  names(res) <- unique(types)
  return(res)

}


get_integrated_datasets <- function(sc_dir) {
  dataset_names <- list.dirs(sc_dir, full.names = FALSE, recursive = FALSE)

  # will have args
  args_paths <- file.path(sc_dir, dataset_names, 'args.json')
  args_exists <- which(file.exists(args_paths))

  # check args for integration type
  integrated <- c()
  for (i in args_exists) {
    dataset_name <- dataset_names[i]
    args <- load_args(sc_dir, dataset_name)
    if ('integration_type' %in% names(args))
      integrated <- c(integrated, dataset_name)
  }

  return(integrated)
}


#' Get choices data.frame for custom metrics (for subsetting)
#'
#' @param scseq \code{SingleCellExperiment}
#'
#' @return data.frame
#' @keywords internal
get_metric_choices <- function(scseq) {

  metrics <- scseq@colData
  names <- colnames(metrics)
  names <- names[sapply(metrics, is.logical)]

  # exclude samples and group
  exclude <- c('is_ctrl', 'is_test', unique(metrics$batch))
  names <- setdiff(names, exclude)

  if (!length(names)) return(NULL)

  levels <- c(TRUE, FALSE)
  choices <- lapply(names, function(name) {
    scseq$cluster <- factor(metrics[[name]], levels = levels)
    get_cluster_choices(levels, scseq = scseq)[1, , drop=FALSE]
  })

  choices <- do.call(rbind, choices)
  choices$name <- stringr::str_trunc(names, 27)
  choices$value <- choices$label <- names
  choices$testColor <- ''
  return(choices)

}

#' Validate a custom metric
#'
#' @param metric String with custom metric to evaluate.
#' @param scseq \code{SingleCellExperiment} to evaluate \code{metric} for.
#'
#' @return Either a data.frame is successful or a object with class 'try-error'.
#'
#' @keywords internal
validate_metric <- function(metric, scseq) {

  allowed <- tryCatch(interpret(metric), error = function(e) return(FALSE))
  if (!allowed) return("Metric not permitted")

  ft <- get_metric_features(metric)
  if (length(ft) == 0) return('No features specified')

  have.ft  <- any(ft %in% c(row.names(scseq), colnames(scseq@colData)))
  if (!have.ft) return('No features specified')

  dat <- try(evaluate_custom_metric(metric, scseq), silent = TRUE)
  return(dat)
}

#' Get features from custom metric
#'
#' @param metric Custom metric to extract features from.
#'
#' @return Character vector of feature names extracted from \code{metric}
#'
#' @keywords internal
get_metric_features <- function(metric) {

  ft <- strsplit(metric, '[|><=\\)\\(&!*%/]')[[1]]
  ft <- gsub(' ', '', ft)
  ft <- ft[ft != '']
  not.num <- is.na(suppressWarnings(as.numeric(ft)))
  ft <- ft[not.num]
  is.quote <- grepl("\"|'", ft)
  ft[!is.quote]
}


#' Evaluate custom single-cell metric
#'
#' @param metric String with metric to evaluate
#' @param scseq \code{SingleCellExperiment} to evaluate \code{metric} for.
#'
#' @return data.frame with column name \code{metric} and result.
#'
#' @keywords internal
evaluate_custom_metric <- function(metric, scseq) {

  ft <- get_metric_features(metric)

  # make sure will run
  expr <- SingleCellExperiment::logcounts(scseq)
  genes <- row.names(expr)[row.names(expr) %in% ft]

  if (length(genes)) {
    expr <- expr[genes,, drop = FALSE]
    expr <- t(as.matrix(expr))
    row.names(expr) <- NULL

  } else {
    expr <- data.frame(expr = rep(NA, ncol(scseq)))
  }

  qcs <- scseq@colData
  row.names(qcs) <- NULL
  qcs <- qcs[, colnames(qcs) %in% ft, drop = FALSE]

  dat <- cbind(expr, qcs)
  colnames(dat) <- c(colnames(expr), colnames(qcs))
  ft <- ft[ft %in% colnames(dat)]

  dat <- dat[, match(ft, colnames(dat)), drop = FALSE]

  # convert features to generic (in case e.g. minus)
  seq <- seq_len(ncol(dat))
  ft.num <- paste0('ft', seq)
  colnames(dat) <- ft.num

  metric.num <- metric
  for (i in seq) metric.num <- gsub(paste0('\\b', ft[i], '\\b'), ft.num[i], metric.num)

  dat <- as.data.frame(dat)
  dat <- within(dat, metric <- eval(rlang::parse_expr(metric.num)))
  dat <- S4Vectors::DataFrame(dat, row.names = colnames(scseq))
  dat <- dat[, 'metric', drop = FALSE]
  colnames(dat) <- metric
  return(dat)
}

cbind.safe <- function(prev, res) {
  res <- cbind(prev, res)
  row.names(res) <- row.names(prev)
  return(res)
}


#' Check safety of evaluating string as R code.
#'
#' @param expr_str String to evaluate,
#' @param max_length maximum length of string to evaluate,
#' @param whitelist R functions to allow use of.
#'
#' @return \code{TRUE} if expression is safe, otherwise \code{FALSE}.
#'
#' @keywords internal
interpret <- function(expr_str,
                      max_length = 200,
                      whitelist = c("&", ">=", "<=", ">", "<", '|', '==', '!=', '*')) {
  safer_eval <- function(expr) {
    if (rlang::is_call(expr)) {
      fn_name <- rlang::call_name(expr)
      if (!fn_name %in% whitelist) return(FALSE)
      do.call(get(fn_name, baseenv()), Map(safer_eval, rlang::call_args(expr)))
    } else if (rlang::is_syntactic_literal(expr)) {
      expr
    }
  }
  if(length(expr_str) >= max_length) return(FALSE)
  parsed <- try(rlang::parse_expr(expr_str), silent = TRUE)
  if (methods::is(parsed, 'try-error')) return(FALSE)
  safer_eval(parsed)
  return(TRUE)
}


#' Format character vectors that will be recognized by HTML
#'
#' @inheritParams base::format
#'
#' @return Formatted \code{x} with &nbsp; substituted for each space.
#'
#' @keywords internal
#'
html_space <- function(x, justify = 'right') {
  if (is.null(x) | !length(x)) return(NULL)
  width <- max(nchar(x))*7.64
  sprintf('<span style="width:%spx;display:inline-block;">%s</span>', width, x)
}

#' Get and save cluster stats for single-cell related selectizeInputs
#'
#' @param resoln_dir Sub directory with single cell dataset info specific to resolution.
#' @param contrast_dir Sub directory with single cell dataset info specific to contrast.
#' @param scseq \code{SingleCellExperiment} object to get/save stats for.
#'   if \code{NULL} (Default), will be loaded.
#' @param top_tables List of \code{limma::topTable} results used by
#'   scSampleComparison for integrated datasets.
#' @param use_disk Should results by saved to disk? used by scSampleComparison
#'   to persist results to disk.
#' @inheritParams get_cluster_choices
#'
#' @return List with cluster stats
get_cluster_stats <- function(resoln_dir = NULL, scseq = NULL, top_tables = NULL, use_disk = FALSE, sample_comparison = FALSE, contrast_dir = NULL) {

  stats_dir <- resoln_dir
  if (!is.null(contrast_dir)) stats_dir <- contrast_dir

  stats_path <- file.path(stats_dir, 'cluster_stats.qs')
  if (file.exists(stats_path) && use_disk) return(qs::qread(stats_path))

  if (is.null(scseq)) {
    dataset_dir <- dirname(resoln_dir)
    scseq <- load_scseq_qs(dataset_dir)
    scseq <- attach_clusters(scseq, resoln_dir)
    if (sample_comparison) scseq <- scseq[, scseq$orig.ident %in% c('test', 'ctrl')]
  }

  ncells <- c(tabulate(scseq$cluster), ncol(scseq))
  pcells <- ncells / ncol(scseq) * 100
  stats <- list(ncells = ncells, pcells = pcells)

  if (sample_comparison) {

    # number of total test and ctrl cells (shown)
    nbins <- length(levels(scseq$cluster))
    is.test <- scseq$orig.ident == 'test'
    stats$ntest <- c(tabulate(scseq$cluster[is.test], nbins = nbins), sum(is.test))
    stats$nctrl <- c(tabulate(scseq$cluster[!is.test], nbins = nbins), sum(!is.test))

    # number of test and ctrl cells in each sample (title)
    test <- unique(scseq$batch[is.test])
    ctrl <- unique(scseq$batch[!is.test])

    neach <- tapply(scseq$cluster, list(scseq$batch, scseq$cluster), length)
    neach[is.na(neach)] <- 0
    neach <- cbind(neach, apply(neach, 1, sum))

    stats$ntest_each <- apply(neach[test,, drop = FALSE], 2, paste, collapse = '-')
    stats$nctrl_each <- apply(neach[ctrl,, drop = FALSE], 2, paste, collapse = '-')
  }

  # number of significant differentially expressed genes in each cluster (pseudobulk)
  samples <- scseq$batch[scseq$orig.ident %in% c('test', 'ctrl')]
  reps <- length(unique(samples)) > 2

  if (sample_comparison & reps) {
    nsig <- rep(0, nbins)
    names(nsig) <- seq_len(nbins)

    test_clusters <- names(top_tables)
    nsig[test_clusters] <- sapply(top_tables, function(tt) {sum(tt$adj.P.Val < 0.05)})
    stats$nsig <- nsig
  }

  # show number of significant with logFC > 1
  if (sample_comparison) {
    nbig <- rep(0, nbins)
    names(nbig) <- seq_len(nbins)

    test_clusters <- names(top_tables)
    nbig[test_clusters] <- sapply(top_tables, function(tt) {sum(abs(tt$logFC) > 1)})
    stats$nbig <- nbig
  }

  if (use_disk) qs::qsave(stats, stats_path)
  return(stats)
}



#' Get contrast choices data.frame for selectize dropdown
#'
#' @param clusters Character vector of cluster names
#' @param test Name of test contrast
#'
#' @return data.frame with columns for rendering selectizeInput contrast choices
#'
#' @keywords internal
get_contrast_choices <- function(clusters, test) {

  # group choices are as compared to other clusters
  test_name <- clusters[as.numeric(test)]
  ctrl_names <- clusters[clusters != test_name]
  ctrls <- setdiff(seq_along(clusters), test)

  colours <- get_palette(clusters, with_all = TRUE)
  colours <- colours[seq_along(clusters)]
  names(colours) <- clusters

  contrast_choices <- data.frame(test = test,
                                 ctrl = c('all', ctrls),
                                 name = test_name,
                                 value = c(test, paste0(test, '-vs-', ctrls)),
                                 testColor = colours[test_name],
                                 ctrlColor = c('white', colours[ctrl_names]), row.names = NULL, stringsAsFactors = FALSE)



  contrast_choices$title <- paste(test_name, 'vs', c('all', ctrl_names))

  # handle single cluster
  if (!length(ctrls))
    contrast_choices <- contrast_choices[1, ]

  return(contrast_choices)

}


#' Get markers for one cluster against one other cluster
#'
#' @param con String identifying clusters to compare e.g. \code{'1vs2'}
#' @param dataset_dir Path to folder for dataset
#'
#' @return Named list with data.frames for each comparison direction
#'
#' @keywords internal
#'
get_contrast_markers <- function(con, markers) {
  con <- strsplit(con, '-vs-')[[1]]

  test <- markers[[con[1]]]
  ctrl <- markers[[con[2]]]

  common <- intersect(test$feature, ctrl$feature)
  ctrl <- ctrl[match(common, ctrl$feature), ]
  test <- test[match(common, test$feature), ]

  test$auc_diff <- test$auc - ctrl$auc
  test$pct_diff <- test$pct_in - ctrl$pct_in
  test$pct_out <- ctrl$pct_in
  test$auc <- NULL
  test <- test[order(test$pct_diff, decreasing = TRUE), ]

  return(test)
}


get_qc_table <- function(qc_metrics = NULL) {

  is.score <- names(qc_metrics) == 'numeric'
  html_metrics <- stringr::str_trunc(qc_metrics, 40, 'left')
  html_metrics[is.score] <- sprintf(
    "<span class='input-swatch' style='background: linear-gradient(to right, %s, %s);'></span><span title='%s'>%s</span>",
    const$colors$qc[1], const$colors$qc[2], qc_metrics[is.score], html_metrics[is.score])

  html_metrics[!is.score] <- sprintf(
    "<span class='input-swatch' style='background: linear-gradient(to right, %s 50%%, %s 50%%);'></span><span title='%s'>%s</span>",
    const$colors$bool[1], const$colors$bool[2], qc_metrics[!is.score], html_metrics[!is.score])



  data.table::data.table(Feature = html_metrics,
                         feature = qc_metrics)
}

get_gene_table <- function(markers,
                           species = 'Homo sapiens',
                           tx2gene = NULL) {

  # check if markers is limma::topTable
  is.tt <- !any(c('auc_diff', 'auc') %in% colnames(markers))

  # top markers uses features as genes and group as type
  features <- if (is.tt) row.names(markers) else markers$feature

  # add gene description
  if (is.null(tx2gene)) tx2gene <- dseqr.data::load_tx2gene(species)

  idx <- match(features, tx2gene$gene_name)
  desc <- tx2gene$description[idx]
  desc[is.na(desc)] <- features[is.na(desc)]

  # gene choices
  html_features <- sprintf('<span title="%s"><a target="_blank" href="%s%s">%s</a></span>',
                           desc,
                           'https://www.genecards.org/cgi-bin/carddisp.pl?gene=',
                           features, stringr::str_trunc(features, width=14))

  # swatch indicates group that is marker for
  group <- markers$group
  levs <- setdiff(unique(group), '')
  if (length(levs)>1) {
    pal <- get_palette(levs)
    names(pal) <- levs
    has.grp <- group != ''
    colors <- pal[group][has.grp]
    html_features[has.grp] <- paste0(
      sprintf('<span>%s: </span><span class="input-swatch" style="background: %s;"></span>',
              names(colors), colors),
      html_features[has.grp]
    )
  }

  if (is.tt) {

    table <- data.table::data.table(
      Feature = html_features,
      'logFC' = markers$logFC,
      'FDR' = markers$adj.P.Val,
      'N<0.5' = markers$`N<0.5`,
      feature = features
    )

  } else {

    table <- data.table::data.table(
      Feature = html_features,
      '\U{0394}AUC' = markers$auc_diff,
      'AUC' = markers[['auc']],
      '\U{0394}%' = markers$pct_diff,
      '%IN' = markers$pct_in,
      '%OUT' = markers$pct_out,
      feature = features
    )

    is.con <- 'auc_diff' %in% colnames(markers)
    if (is.con) {
      data.table::setnames(table, '%IN', '%A')
      data.table::setnames(table, '%OUT', '%B')
    }
  }
  return(table)
}

get_leftover_table <- function(features, species = 'Homo sapiens', tx2gene = NULL) {

  # add gene description
  if (is.null(tx2gene)) tx2gene <- dseqr.data::load_tx2gene(species)

  idx <- match(features, tx2gene$gene_name)
  desc <- tx2gene$description[idx]
  desc[is.na(desc)] <- features[is.na(desc)]

  # gene choices
  html_features <- sprintf('<span title="%s"><a target="_blank" href="%s%s">%s</a></span>',
                           desc,
                           'https://www.genecards.org/cgi-bin/carddisp.pl?gene=',
                           features, features)
  data.table::data.table(
    Feature = html_features,
    feature = features
  )
}

construct_top_markers <- function(markers, scseq) {

  # top 10 markers per cluster (allow repeats)
  top <- lapply(markers, utils::head, 10)
  top <- data.table::rbindlist(top)

  # rest of genes
  rest <- data.table::data.table(group = '',
                                 feature = setdiff(row.names(scseq), top$feature))

  rbind(top, rest, fill = TRUE)
}



validate_up_meta <- function(res, ref) {
  msg <- NULL
  groups <- stats::na.exclude(res$group)

  if (length(unique(groups)) < 2) {
    msg <- 'At least two group names needed.'
  }

  return(msg)
}


#' Integrate previously saved SingleCellExperiments
#'
#' Performs integration and saves as a new analysis.
#' Used by \code{explore_scseq_clusters} shiny app.
#'
#' @param sc_dir Directory with saved single-cell datasets.
#' @param test Character vector of test analysis names.
#' @param ctrl Character vector of control analysis names.
#' @param exclude_clusters Charactor vector of clusters for excluding cells. Only included to save to args.
#' @param exclude_cells Character vector of cell names to exclude.
#' @param integration_name Name for new integrated analysis.
#' @param integration_type Charactor vector of one or more integration types.
#' @param subset_metrics Metrics to subset based on.
#' @param is_include Boolean - are cells that match \code{subset_metrics} included or excluded?
#' @param progress optional Shiny \code{Progress} object.
#'
#' @seealso \code{\link{run_fastmnn}} \code{\link{run_harmony}}
#'
#' @return TRUE is successful, otherwise FALSE
#'
#' @keywords internal
integrate_saved_scseqs <- function(
    sc_dir,
    integration_name,
    dataset_names = NULL,
    scseqs = NULL,
    integration_type = c('harmony', 'fastMNN', 'Azimuth', 'symphony'),
    exclude_clusters = NULL,
    exclude_cells = NULL,
    subset_metrics = NULL,
    is_include = NULL,
    founder = integration_name,
    pairs = NULL,
    hvgs = NULL,
    ref_name = NULL,
    npcs = 30,
    cluster_alg = 'leiden',
    resoln = 1,
    progress = NULL,
    value = 0,
    tx2gene_dir = NULL) {

  # for save_scseq_args
  args <- c(as.list(environment()))
  args$progress <- args$sc_dir <- NULL
  args$date <- Sys.time()
  args$value <- NULL
  args$scseqs <- NULL

  if (is.null(progress)) {
    progress <- list(set = function(value, message = '', detail = '') {
      cat(value, message, detail, '...\n')
    })
  }

  dataset_name <- paste(integration_name, integration_type, sep = '_')

  if (is.null(scseqs)) {
    progress$set(value+1, detail = 'loading')
    scseqs <- load_scseq_subsets(dataset_names,
                                 sc_dir,
                                 subset_metrics,
                                 is_include,
                                 with_logs = TRUE,
                                 with_counts = TRUE)

  }

  if (!length(scseqs)) {
    progress$set(value+1, detail = 'error: no cells left')
    Sys.sleep(3)
    return(FALSE)
  }

  # make sure all the same species
  species <- unique(sapply(scseqs, function(x) x@metadata$species))
  if(length(species) > 1) stop('Multi-species integration not supported.')

  progress$set(value+2, detail = 'integrating')
  combined <- integrate_scseqs(scseqs,
                               type = integration_type,
                               hvgs = hvgs,
                               ref_name = ref_name,
                               species = species,
                               tx2gene_dir = tx2gene_dir)

  # keep successfully combined
  scseqs <- scseqs[unique(combined$batch)]

  combined$project <- dataset_name

  # retain original QC metrics
  combined <- add_combined_metrics(combined, scseqs)

  # add ambience info
  combined <- add_combined_ambience(combined, scseqs)
  rm(scseqs); gc()

  combined@metadata$species <- species
  combined@metadata$npcs <- npcs

  # save what is stable with resolution change
  scseq_data <- list(species = species,
                     pairs = pairs,
                     founder = founder,
                     resoln = resoln)

  is_ref <- integration_type %in% c('Azimuth', 'symphony')
  if (is_ref) {
    # default resoln for each azimuth/symphony ref
    resoln <- get_ref_resoln(ref_name)

    scseq_data$resoln <- resoln
    scseq_data$ref_name <- ref_name
    scseq_data$scseq <- combined

    save_scseq_data(scseq_data, dataset_name, sc_dir)
    save_ref_clusters(combined@colData, dataset_name, sc_dir)

  } else {
    # UMAP/TSNE on corrected reducedDim
    progress$set(value+3, detail = 'reducing')
    combined <- run_reduction(combined, dimred = 'corrected')

    # get graph
    scseq_data$snn_graph <- get_snn_graph(combined, npcs = npcs)
    combined$cluster <- get_clusters(scseq_data$snn_graph, cluster_alg, resoln)

    scseq_data$scseq <- combined
    save_scseq_data(scseq_data, dataset_name, sc_dir)
  }

  # run things that can change with resolution change
  run_post_cluster(combined,
                   dataset_name,
                   sc_dir,
                   resoln,
                   progress,
                   value+4,
                   reset_annot = !is_ref)

  save_scseq_args(args, dataset_name, sc_dir)

  return(TRUE)
}

#' Post-Cluster Processing of SingleCellExperiment
#'
#' Performs the following:
#' - downsampling (used for label transfer)
#' - pseudobulk
#' - saving data
#'
#' @param scseq SingleCellExperiment
#' @param dataset_name Name of dataset
#' @param sc_dir Directory with single-cell datasets
#' @param resoln resolution parameter
#' @param progress progress
#' @param value value
#' @param reset_annot if \code{TRUE} then overwrites saved annotation with cluster indices.
#'   Set to \code{FALSE} when Azimuth is used as cluster annotation is defined.
#'
#' @return NULL
#'
run_post_cluster <- function(scseq, dataset_name, sc_dir, resoln = 1, progress = NULL, value = 0, reset_annot = TRUE) {

  if (is.null(progress)) {
    progress <- list(set = function(value, message = '', detail = '') {
      cat(value, message, detail, '...\n')
    })
  }

  # used for label transfer
  scseq_sample <- downsample_clusters(scseq)

  progress$set(value+1, detail = 'pseudobulk')

  summed <-  aggregate_across_cells(scseq)

  anal <- list(scseq_sample = scseq_sample,
               clusters = scseq$cluster,
               summed = summed,
               applied = TRUE)

  if (reset_annot) anal$annot <- levels(scseq$cluster)

  # save in subdirectory e.g. snn1
  progress$set(value+3, detail = 'saving')
  dataset_subname <- file.path(dataset_name, get_resoln_dir(resoln))
  save_scseq_data(anal, dataset_subname, sc_dir, overwrite = FALSE)

}


# fast pseudobulk using presto
aggregate_across_cells <- function(scseq) {
  counts_mat <- SingleCellExperiment::counts(scseq)
  batch <- scseq$batch
  cluster <- scseq$cluster
  y <- paste(cluster, batch, sep='_')
  levs <- unique(y)

  agg <- presto::sumGroups(counts_mat, factor(y, levs), 1); gc()
  agg <- Matrix::Matrix(agg, sparse = TRUE); gc()
  agg <- Matrix::t(agg); gc()
  row.names(agg) <- row.names(counts_mat)
  colnames(agg) <- levs

  # setup colData
  dup.y <- duplicated(y)
  cdata <- S4Vectors::DataFrame(
    cluster = cluster[!dup.y],
    batch = batch[!dup.y])

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = agg),
    rowData = SingleCellExperiment::rowData(scseq),
    colData = cdata
  )

  return(sce)
}


#' Integration Utility Function for Background Process
#'
#' Used as background process via \link[callr]{r_bg}
#'
#' @inheritParams integrate_saved_scseqs
#' @param integration_types Character vector of integration types to run.
#'
#' @return TRUE is successful, otherwise FALSE
#' @keywords internal
#'
run_integrate_saved_scseqs <- function(sc_dir,
                                       tx2gene_dir,
                                       dataset_names,
                                       integration_name,
                                       integration_types,
                                       ref_name) {

  for (i in seq_along(integration_types)) {

    type <- integration_types[i]
    ref <- NULL
    if (type == 'reference') {
      ref <- ref_name
      type <- get_ref_type(ref)
    }

    # run integration
    res <- integrate_saved_scseqs(sc_dir,
                                  dataset_names,
                                  integration_name = integration_name,
                                  integration_type = type,
                                  ref_name = ref,
                                  value = i*8-8,
                                  tx2gene_dir = tx2gene_dir)

    # stop subsequent integration types if error
    if (!res) return(FALSE)
  }
  return(TRUE)
}

get_ref_type <- function(ref_name) {
  refs[refs$name == ref_name, 'type']
}


#' Subset SingleCellExperiment
#'
#' @param sc_dir Directory with saved single-cell datasets.
#' @param founder Name of original founding dataset (may be ancestor of \code{from_dataset}).
#' @param from_dataset Name of parent dataset.
#' @param dataset_name Name for new subset dataset.
#' @param exclude_clusters Charactor vector of clusters for excluding cells.
#' @param exclude_cells Character vector of cells to exclude. Used for integrated datasets to exclude clusters.
#' @param subset_metrics Character vector of metrics for excluding cells.
#' @param progress Optional shiny progress object.
#'
#' @return TRUE is successful, otherwise FALSE
#'
#' @keywords internal
#'
subset_saved_scseq <- function(sc_dir,
                               founder,
                               from_dataset,
                               dataset_name,
                               exclude_clusters = NULL,
                               subset_metrics = NULL,
                               is_integrated = FALSE,
                               is_include = FALSE,
                               progress = NULL,
                               hvgs = NULL,
                               ref_name = NULL,
                               tx2gene_dir = NULL) {

  if (is.null(progress)) {
    progress <- list(set = function(value, message = '', detail = '') {
      cat(value, message, detail, '...\n')
    })
  }

  # load scseq and subset using metrics
  progress$set(1, detail = 'loading')
  scseq <- load_scseq_subsets(from_dataset, sc_dir, subset_metrics, is_include,
                              with_counts = TRUE, with_logs = TRUE)[[1]]

  if (is.null(scseq)) {
    progress$set(1, detail = "error: no cells")
    Sys.sleep(3)
    return(FALSE)
  }

  # exclude clusters
  scseq <- scseq[, !scseq$cluster %in% exclude_clusters]

  # make repeated subsets order independent:
  # ------
  #
  # for an integrated dataset: re-integrate
  # for a unisample dataset: re-process

  if (is_integrated) {
    args <- load_args(sc_dir, from_dataset)

    # use same integration type as previously unless reference given
    itype <- args$integration_type
    if (!is.null(ref_name)) itype <- get_ref_type(ref_name)

    # use harmony if previous was azimuth/symphony and no ref specified
    if (itype %in% c('Azimuth', 'symphony') & is.null(ref_name)) itype <- 'harmony'

    scseqs <- split_scseq(scseq)
    rm(scseq); gc()

    integrate_saved_scseqs(
      sc_dir,
      scseqs = scseqs,
      integration_name = dataset_name,
      integration_type = itype,
      exclude_clusters = exclude_clusters,
      subset_metrics = subset_metrics,
      founder = founder,
      hvgs = hvgs,
      ref_name = ref_name,
      progress = progress,
      value = 1,
      tx2gene_dir = tx2gene_dir)

  } else {

    # for save_scseq_args
    args <- c(as.list(environment()))
    args$progress <- args$sc_dir <- args$scseq <- NULL
    args$date <- Sys.time()

    # remove previous reduction
    SingleCellExperiment::reducedDims(scseq) <- NULL

    process_raw_scseq(
      scseq,
      dataset_name,
      sc_dir,
      hvgs = hvgs,
      progress = progress,
      value = 1,
      founder = founder,
      ref_name = ref_name,
      tx2gene_dir = tx2gene_dir)

    save_scseq_args(args, dataset_name, sc_dir)
  }

  return(TRUE)
}

load_args <- function(sc_dir, dataset_name) {
  jsonlite::read_json(file.path(sc_dir, dataset_name, 'args.json'), simplifyVector = TRUE)
}

split_scseq <- function(scseq) {

  dataset_names <- unique(scseq$batch)
  scseqs <- list()
  for (dataset_name in dataset_names) {
    scseqi <- scseq[, scseq$batch == dataset_name]
    if (!ncol(scseqi)) next()

    # store ambience for sample
    amb_col <- paste0(dataset_name, '_ambience')
    ambience <- SummarizedExperiment::rowData(scseqi)[[amb_col]]
    SummarizedExperiment::rowData(scseqi)$ambience <- ambience

    scseqs[[dataset_name]] <- scseqi
  }

  return(scseqs)
}


#' Determine the founder of a single-cell dataset
#'
#' Founder is not necesarily the parent, but rather the original ancestor.
#'
#' @param sc_dir Directory with single-cell datasets.
#' @param dataset_name Name of dataset to determine founder for.
#'
#' @return Name of founder.
#' @keywords internal
#'
get_founder <- function(sc_dir, dataset_name) {

  # check for founder of parent
  fpath <- file.path(sc_dir, dataset_name, 'founder.qs')
  founder <- qread.safe(fpath)
  if (is.null(founder)) founder <- dataset_name
  return(founder)
}


#' Save arguments for integration/subsetting
#'
#' @param args Arguments to save
#' @param dataset_name Name of analysis
#' @param sc_dir Directory with single-cell analyses
#'
#' @return NULL
#' @keywords internal
#'
save_scseq_args <- function(args, dataset_name, sc_dir) {
  jsonlite::write_json(args,
                       path = file.path(sc_dir, dataset_name, 'args.json'),
                       auto_unbox = TRUE,
                       null = 'null',
                       pretty = TRUE)
}

#' Add QC metrics to integrated SingleCellExperiment
#'
#' @param combined Integrated \code{SingleCellExperiment}
#' @param scseqs list of \code{SingleCellExperiment}s used to create \code{combined}
#'
#' @return \code{combined} with QC metrics from \code{scseqs} added
#'
#' @keywords internal
add_combined_metrics <- function(combined, scseqs) {
  # add original QC metrics
  metrics <- c('log10_sum', 'log10_detected', 'mito_percent', 'ribo_percent', 'doublet_score')
  for (metric in metrics) combined[[metric]] <- unlist(sapply(scseqs, `[[`, metric), use.names = FALSE)

  return(combined)
}

add_combined_ambience <- function(combined, scseqs) {
  samples <- unique(combined$batch)
  genes <- row.names(combined)

  for (sample in samples) {
    col <- paste0(sample, '_ambience')
    sce <- scseqs[[sample]]
    rdata <- SingleCellExperiment::rowData(sce[genes, ])
    # legacy calculated pct_ambient from droplets < 10
    scol <- ifelse('ambience' %in% colnames(rdata), 'ambience', 'pct_ambient')
    SummarizedExperiment::rowData(combined)[[col]] <- rdata[[scol]]
  }

  return(combined)
}


#' Load subsets of scseqs for integration/subsetting
#'
#' Sets orig.ident to \code{ident}.
#'
#' @param dataset_names Character vector of single cell analysis names to load.
#' @param sc_dir The directory with single-cell datasets
#' @param ident Either \code{'test'} \code{'ctrl'} for integration. The from dataset for subsetting.
#'
#' @return List of \code{SingleCellExperiment} objects.
#'
#' @keywords internal
load_scseq_subsets <- function(dataset_names, sc_dir, subset_metrics = NULL, is_include = NULL, ...) {

  scseqs <- list()

  # load each scseq and exclude based on metrics
  for (dataset_name in dataset_names) {
    dataset_dir <- file.path(sc_dir, dataset_name)
    scseq <- load_scseq_qs(dataset_dir, ...)

    # set orig.resoln to track resolution of origin datasets
    resoln_name <- load_resoln(dataset_dir)
    scseq$orig.resoln <- resoln_name

    # set orig.ident to original dataset name (subset)
    scseq$orig.ident <- factor(dataset_name)

    if (length(subset_metrics)) {
      cdata <- scseq@colData
      metrics <- lapply(subset_metrics, evaluate_custom_metric, scseq)
      metrics <- do.call(cbind, metrics)
      if (!is.null(metrics)) cdata <- cbind.safe(cdata, metrics)

      exclude <- cdata[, subset_metrics, drop = FALSE]

      for (i in seq_len(ncol(exclude)))
        if (is_include) exclude[,i] <- !exclude[,i]

      exclude <- apply(exclude, 1, any)
      scseq <- scseq[, !exclude]
      gc()

      # remove subset_metrics hardcoded in scseq
      # no/all cells will meet them by definition
      scseq@colData <- scseq@colData[, !colnames(scseq@colData) %in% subset_metrics, drop = FALSE]
    }


    # require that have cells left
    if (ncol(scseq)>0) scseqs[[dataset_name]] <- scseq
  }

  if (!length(scseqs)) return(NULL)
  return(scseqs)
}

get_exclude_cells <- function(scseq, exclude_clusters) {
  # subset scseq to excluded cells
  exclude_num <- gsub('^.+?_(\\d+)$', '\\1', exclude_clusters)
  is.exclude <- as.numeric(scseq$cluster) %in% exclude_num
  scseq <- scseq[, is.exclude]

  # remove de-dup postfix from integration
  cells <- gsub('[.]\\d+$', '', colnames(scseq))

  # group excluded by dataset name
  exclude <- split(cells, scseq$batch)

  return(exclude)
}


#' Save Single Cell RNA-seq data for app
#'
#' @param scseq_data Named list with \code{scseq}, \code{markers}, and/or \code{annot}
#' @param dataset_name The analysis name.
#' @param sc_dir Path to directory with single-cell datasets.
#' @param overwrite overwrite \code{dataset_name} sub-directory?
#'
#' @return NULL
#' @keywords internal
save_scseq_data <- function(scseq_data, dataset_name, sc_dir, overwrite = TRUE) {
  dataset_dir <- file.path(sc_dir, dataset_name)

  # remove all previous data in case overwriting
  if (overwrite) unlink(dataset_dir, recursive = TRUE)

  dir.create(dataset_dir, showWarnings = FALSE)
  for (type in names(scseq_data)) {
    item <- scseq_data[[type]]

    if (type == 'scseq') {
      dataset_dir <- file.path(sc_dir, dataset_name)
      split_save_scseq(item, dataset_dir)

    } else {
      qs::qsave(item, scseq_part_path(sc_dir, dataset_name, type), preset = 'fast')
    }
  }

  return(NULL)
}

#' Load SingleCellExperiment from qs file
#'
#' Also attaches clusters from last applied leiden resolution and stores resolution.
#'
#' @param dataset_dir Path to folder with scseq.qs file
#' @param meta data.frame with column \code{group} and \code{row.names} as sample
#' names corresponding to \code{scseq$batch}. Default (\code{NULL}) loads previous
#' specification from file.
#' @param groups character vector of length two. First value is test group name
#' and second in control group name.
#' @param with_logs should logcounts be loaded? Default is \code{FALSE} to increase
#' speed and reduce memory usage.
#' @param with_counts should counts be loaded? Default is \code{FALSE} to increase
#' speed and reduce memory usage.
#'
#' @return SingleCellExperiment
#' @export
#'
load_scseq_qs <- function(dataset_dir, meta = NULL, groups = NULL, with_logs = FALSE, with_counts = FALSE) {
  shell_path <- file.path(dataset_dir, 'shell.qs')
  scseq <- qs::qread(shell_path)

  if (with_logs) {
    dgclogs <- qs::qread(file.path(dataset_dir, 'dgclogs.qs'))
    SingleCellExperiment::logcounts(scseq) <- dgclogs
  }

  if (with_counts) {
    counts <- qs::qread(file.path(dataset_dir, 'counts.qs'))
    SingleCellExperiment::counts(scseq) <- counts
  }

  resoln_name <- load_resoln(dataset_dir)
  scseq <- attach_clusters(scseq, file.path(dataset_dir, resoln_name))
  scseq <- attach_meta(scseq, dataset_dir, meta, groups)
  if (is.null(scseq$batch)) scseq$batch <- scseq$project

  return(scseq)
}

#' Validate dataset selection for integration
#'
#'
#' @return \code{NULL} if valid, otherwise an error message
#'
#' @keywords internal
validate_integration <- function(types, name, ref_name, dataset_names, sc_dir) {
  msg <- NULL

  species <- get_integration_species(dataset_names, sc_dir)


  if ('reference' %in% types & !isTruthy(ref_name)) {
    msg <- 'Select reference'
  } else if (is.null(types)) {
    msg <- 'Select one or more integration types'
  } else if (name == '') {
    msg <- 'Enter name for integrated dataset'
  } else if (grepl('/', name)) {
    msg <- "Remove '/' from dataset name"
  } else if (length(dataset_names) < 2) {
    msg <- 'Select atleast two datasets'
  } else if (length(species) > 1) {
    msg <- 'Datasets must be from the same species'
  }

  return(msg)
}

get_integration_species <- function(integration_datasets, sc_dir) {
  species_paths <- file.path(sc_dir, integration_datasets, 'species.qs')
  species <- lapply(species_paths, qread.safe)
  unique(unlist(species))
}

validate_subset <- function(from_dataset, subset_name, subset_features, is_include, hvgs) {

  have.dataset <- shiny::isTruthy(from_dataset)
  if (!have.dataset) {
    return('No dataset selected')
  }

  have.name <- shiny::isTruthy(subset_name)
  if (!have.name) {
    return('No name given')
  }


  all.excluded <- is_include & !shiny::isTruthy(subset_features)
  if (all.excluded) {
    return('All excluded')
  }

  return(validate_not_path(subset_name))
}


#' Get path to saved scseq part
#'
#' @param data_dir Path to directory with analyses.
#' @param dataset_name Name of analysis.
#' @param part either \code{'annot'}, \code{'scseq'}, or \code{'markers'}.
#'
#' @return Path to analysis \code{part}.
#'
#' @keywords internal
scseq_part_path <- function(data_dir, dataset_name, part) {
  fname <- paste0(part, '.qs')
  file.path(data_dir, dataset_name, fname)
}

#' Run drug queries
#'
#' Used by run_comparison
#'
#' @param res_paths List of paths with names \code{'cmap'} \code{'l1000'} and \code{'anal'} to
#' saved cmap, l1000, and differential expression analysis results.
#' @param session Shiny session object used for progress bar.
#'
#' @return \code{res} with drug query results added to \code{'cmap'} \code{'l1000'} slots.
#'
#' @keywords internal
run_drug_queries <- function(top_table, drug_paths, es, ngenes = 200) {

  # get dprime effect size values for analysis
  dprimes <- get_dprimes(top_table)

  # get correlations between query and drug signatures
  res <- list(
    cmap = query_drugs(dprimes, es$cmap, ngenes = ngenes),
    l1000_drugs = query_drugs(dprimes, es$l1000_drugs, ngenes = ngenes),
    l1000_genes = query_drugs(dprimes, es$l1000_genes, ngenes = ngenes)
  )

  qs::qsave(res$cmap, drug_paths$cmap)
  qs::qsave(res$l1000_drugs, drug_paths$l1000_drugs)
  qs::qsave(res$l1000_genes, drug_paths$l1000_genes)

  return(res)
}

get_celldex_species <- function(ref_name) {
  ifelse(ref_name %in% c('ImmGenData', 'MouseRNAseqData'), 'Mus musculus', 'Homo sapiens')
}

convert_species <- function(x, tx2gene_dir, species, other_species = 'Homo sapiens') {
  if (methods::is(x, 'SingleCellExperiment')) species <- x@metadata$species
  if (species == other_species) return(x)

  species_tx2gene <- load_tx2gene(species, tx2gene_dir)
  other_tx2gene <- load_tx2gene(other_species, tx2gene_dir)

  other <- species_symbols_to_other(
    symbols = row.names(x),
    species_tx2gene = species_tx2gene,
    other_tx2gene = other_tx2gene)

  na.other <- is.na(other)
  x <- x[!na.other, ]
  row.names(x) <- other[!na.other]

  return(x)
}

species_symbols_to_other <- function(symbols, species_tx2gene, other_tx2gene) {

  # df with species gene name and hgnc homologous ensemble id
  species_tx2gene <-
    species_tx2gene %>%
    dplyr::filter(.data$gene_name %in% symbols) %>%
    dplyr::select(.data$gene_name, .data$hsapiens_homolog_ensembl_gene) %>%
    stats::na.omit() %>%
    dplyr::filter(!duplicated(.data$gene_name))

  # ensure rows have same order as symbols
  species_tx2gene <-
    data.frame(gene_name = symbols) %>%
    dplyr::left_join(species_tx2gene) %>%
    dplyr::rename('species_symbol' = 'gene_name')


  # df with other symbol and hgnc homologous ensemble id
  other_tx2gene <- other_tx2gene %>%
    dplyr::select(.data$gene_name, .data$hsapiens_homolog_ensembl_gene) %>%
    dplyr::filter(!duplicated(.data$hsapiens_homolog_ensembl_gene)) %>%
    dplyr::rename('other_symbol' = 'gene_name') %>%
    stats::na.omit()

  # join the two
  map <- dplyr::left_join(species_tx2gene, other_tx2gene)

  return(map$other_symbol)
}

#' Load drug effect size matrices for drug queries
#'
#' @return list of matrices
#'
#' @keywords internal
load_drug_es <- function() {

  cmap  <- dseqr.data::load_data('cmap_es_ind.qs')
  l1000_drugs  <- dseqr.data::load_data('l1000_drugs_es.qs')
  l1000_genes  <- dseqr.data::load_data('l1000_genes_es.qs')

  return(list(
    cmap = cmap,
    l1000_drugs = l1000_drugs,
    l1000_genes = l1000_genes
  ))
}

#' Used to generate file names for single cell analyses.
#'
#' @param x Character vector of selected cluster numbers
#' @return String with sorted cluster numbers comma collapsed.
#'
#' @keywords internal
collapse_sorted <- function(x, collapse = ',') {
  paste(sort(as.numeric(x)), collapse = ',')
}


#' Search for closest row in a truth matrix for each row in a test matrix.
#'
#' from Stack Overflow: 40668623
#'
#' @param truth matrix to pick nearest from for all rows in \code{test}.
#' @param test matrix to test each row of and find closest rows in \code{truth}.
#'
#' @return vector of rows in \code{truth} that are the nearest to each row in \code{test}.
#' @keywords internal
#'
#' @examples
#' set.seed(123)  ## for reproducibility
#' D <- 2 #amount of dimensions
#' K <- 5
#' events <- 2*K #number of events
#' truth <- matrix(data=runif(events, min = 0, max = 1), nrow=K)
#' E <- 2
#' test <- matrix(data=runif(2*E, min = 0, max = 1), nrow=E)
#'
#' \dontrun{
#' #[1] 4 3
#' get_nearest_row(truth, test)
#' }
#'
get_nearest_row <- function(truth, test) {
  ntest <- nrow(test)
  ntruth <- nrow(truth)
  diffs   <- truth[rep(seq_len(ntruth), ntest),] - test[rep(seq_len(ntest), each=ntruth),]
  eucdiff <- function(x) sqrt(rowSums(x^2))
  max.col(-matrix(eucdiff(diffs), nrow=ntest, byrow=TRUE), "first")
}


#' Get median x-y coordinates for clusters in TSNE plot data
#'
#' @param plot_data data.frame with columns \code{'ident'}, \code{'TSNE_1'}, and \code{'TSNE_2'}.
get_tsne_coords <- function(plot_data) {
  ident <- TSNE_1 <- TSNE_2 <- ident <- NULL

  plot_data %>%
    dplyr::group_by(ident) %>%
    dplyr::summarise(TSNE_1 = stats::median(TSNE_1),
                     TSNE_2 = stats::median(TSNE_2))
}


#' Update Progress from Background Processes
#'
#' @param bgs \code{reactivevalues} of \link[callr]{r_bg}
#' @param progs \code{reactivevalues} of \link[shiny]{Progress}
#' @param new_dataset \code{reactive} that triggers update of available datasets.
#'
#' @return NULL
#' @keywords internal
#'
handle_sc_progress <- function(bgs, progs, new_dataset) {
  bg_names <- names(bgs)
  if (!length(bg_names)) return(NULL)

  todel <- c()
  for (name in bg_names) {
    bg <- bgs[[name]]
    progress <- progs[[name]]
    if(is.null(bg)) next

    msgs <- bg$read_output_lines()
    if (length(msgs)) print(msgs)

    # for some reason this un-stalls bg process for integration
    # also nice to see things printed to stderr
    errs <- bg$read_error_lines()
    if (length(errs)) print(errs)

    if (length(msgs)) {
      last <- utils::tail(msgs, 1)
      progress$set(value = progress$getValue() + length(msgs),
                   detail = gsub('^\\d+ +', '', last))

    }

    if (!bg$is_alive()) {
      res <- bg$get_result()
      progress$close()
      if (res) new_dataset(name)
      todel <- c(todel, name)
    }
  }
  for (name in todel) {
    bgs[[name]] <- NULL
    progs[[name]] <- NULL
  }
}


#' Transfer annotation when change resolution
#'
#' @param resoln new resolution
#' @param prev_resoln previous resolution
#' @param dataset_name name of single-cell dataset
#' @param sc_dir directory with folder named \code{dataset_name}
#'
#' @return saves new annotation in resolution subdirectory
#' @keywords internal
#'
transfer_prev_annot <- function(resoln, prev_resoln, dataset_name, sc_dir) {

  # ref clusters are from previous resolution
  ref_resoln_name <- file.path(dataset_name, get_resoln_dir(prev_resoln))
  ref_cluster <- qs::qread(file.path(sc_dir, ref_resoln_name, 'clusters.qs'))

  # query clusters are new resolution
  query_resoln_name <- file.path(dataset_name, get_resoln_dir(resoln))
  query_cluster <- qs::qread(file.path(sc_dir, query_resoln_name, 'clusters.qs'))

  # transfer labels
  tab <- table(assigned = ref_cluster, cluster = query_cluster)
  pred <- row.names(tab)[apply(tab, 2, which.max)]
  annot <- get_pred_annot(pred, ref_resoln_name, query_resoln_name, sc_dir)
  annot_nums <- as.character(seq_along(annot))

  # keep ordered nums where prediction is numeric
  suppressWarnings(is.num <- !is.na(as.numeric(remove.unique(annot))))
  annot[is.num] <- annot_nums[is.num]

  qs::qsave(annot, file.path(sc_dir, query_resoln_name, 'annot.qs'))
  return(annot)
}

#' Get the applied resolution dataset name
#'
#' e.g. PBMCS_1/snn1
#'
#' @param sc_dir directory with \code{dataset_name} folder
#' @param dataset_name name of single-cell dataset
#'
#' @return the resolution dataset name based on the last applied resolution
#' @keywords internal
#'
get_resoln_name <- function(sc_dir, dataset_name) {
  dataset_dir <- file.path(sc_dir, dataset_name)

  # so that don't add '/snn1' to e.g. 'reset' and 'BlueprintEncodeData'
  if (!dir.exists(dataset_dir)) return(dataset_name)

  resoln <- load_resoln(dataset_dir)
  file.path(dataset_name, resoln)
}

#' Load the applied resolution name
#'
#' e.g. snn1
#'
#' @param dataset_dir directory corresponding to single-cell dataset
#'
#' @return the last applied resolution name (e.g. \code{'snn1'})
#' @keywords internal
#'
load_resoln <- function(dataset_dir) {
  resoln_path <- file.path(dataset_dir, 'resoln.qs')
  resoln <- qread.safe(resoln_path, .nofile = 1)
  get_resoln_dir(resoln)
}


#' Run limma on pseudobulk experiment
#'
#' @param meta data.frame of sample metadata with column \code{group} and \code{row.names}
#'   set with \code{summed$batch}
#' @param species species name used to retrieve feature annotation.
#' @param trend Should limma-trend analysis be run (default is \code{FALSE})? Slightly faster than
#'  limma-voom (default). Used for pseudobulk grid differential expression analyses.
#' @param summed pseudobulk \code{SingleCellExperiment}. If \code{NULL}, \code{summed_path} must be supplied
#'   and the call is assumed to originate from \link[callr]{r_bg}.
#' @param with_fdata if \code{TRUE}, adds Entrez IDs to fit object which are needed to
#'  downstream GO pathway analyses.
#' @param progress progress object or \code{NULL}.
#' @param value Initial value of progress.
#' @param ... additional arguments to \link[edgeR]{filterByExpr}.
#' @inheritParams edgeR::calcNormFactors
#'
#' @return List with a fit object, model matrix, and normalized expression matrix.
#' @keywords internal
#'
run_limma_scseq <- function(summed, meta, species, trend = FALSE, method = 'TMMwsp', with_fdata = FALSE, progress = NULL, value = 0, ...) {


  if (is.null(progress)) {
    progress <- list(set = function(value, message = '', detail = '') {
      cat(value, message, detail, '...\n')
    })
  }

  groups <- meta$group
  names(groups) <- row.names(meta)

  # need at least two groups
  gtab <- table(groups)
  if (length(gtab) < 2) return(NULL)

  # add groups
  idents <- unname(groups[summed$batch])
  summed$orig.ident <- factor(idents)

  progress$set(value = value+1)
  subsets <- construct_pbulk_subsets(summed, method, ...)

  # need entrezids for pathway analyses
  fdata <- NULL
  if (with_fdata) {
    # get feature annotation
    release <- switch(species,
                      'Homo sapiens' = '94',
                      'Mus musculus' = '98',
                      '103')

    fdata <- rkal::setup_fdata(species, release)
    fdata <- unique(fdata, by = c('gene_name'))
    fdata <- fdata[row.names(summed), ]
    fdata <- as.data.frame(fdata)
    row.names(fdata) <- fdata$gene_name
  }

  # fitting models
  type <- ifelse(trend, 'grid:', 'cluster:')
  progress$set(message = paste('Fitting', type), detail = '', value = value+2)
  inc <- 2/length(subsets)

  fit <- list()
  for (i in seq_along(subsets)) {
    clust <- names(subsets)[i]
    progress$set(detail = clust)
    progress$inc(inc)
    subset <- subsets[[i]]
    fit[[clust]] <- fit_lm(subset, fdata, clust, trend)
  }

  progress$set(value = value+4)
  return(fit)
}

fit_lm <- function(subset, fdata, clust, trend = FALSE) {
  pdata <- subset$pdata
  group <- pdata$group
  mod <- stats::model.matrix(~0 + group)
  colnames(mod) <- gsub("^group", "", colnames(mod))

  lib.size <- pdata$lib.size * pdata$norm.factors

  # fit for unfiltered genes
  yik <- subset$yik
  if (trend) {
    v <- t(log2(t(yik + 0.5)/(lib.size + 1) * 1e+06))
    aw <- limma::arrayWeights(v, mod, method = "reml")
    fit <- limma::lmFit(v, mod, weights = aw)
    E <- NULL

  } else {
    v <- quickVoomWithQualityWeights(yik, mod, lib.size)
    fit <- limma::lmFit(v, mod)
    E <- v$E
    colnames(E) <- gsub(paste0('^', clust, '_'), '', colnames(E))
  }

  if (!is.null(fdata)) {
    fit$genes <- data.frame(fdata[row.names(fit), 'ENTREZID', drop = FALSE])
  }
  return(list(fit=fit, mod=mod, E=E))
}

quickVoomWithQualityWeights <- function (y, mod, lib.size) {
  # need log-expressed values for arrayWeights
  v <- t(log2(t(y + 0.5)/(lib.size + 1) * 1e+06))
  aw <- limma::arrayWeights(v, mod, method = "reml")
  suppressWarnings(v <- limma::voom(y, mod, lib.size, weights = aw))

  #	incorporate the array weights into the voom weights
  v$weights <- t(aw * t(v$weights))
  v$targets$sample.weights <- aw

  return(v)
}


get_grid <- function(scseq) {
  reds <- SingleCellExperiment::reducedDimNames(scseq)
  red <- reds[reds %in% c('UMAP', 'TSNE')]

  red.mat <- SingleCellExperiment::reducedDim(scseq, red)
  grid_size <- get_grid_size(red.mat)
  nx <- grid_size[1]
  ny <- grid_size[2]

  dat <- data.frame(x=red.mat[,1], y=red.mat[,2])

  # create grid
  x <- seq(min(dat$x), max(dat$x), length.out = nx)
  y <- seq(min(dat$y), max(dat$y), length.out = ny)

  # in nx*ny grid: get xy bin for each point
  points <- cbind(dat$x, dat$y)
  grid <- data.frame(xi = findInterval(points[,1], x),
                     yi = findInterval(points[,2], y))

  diff.x <- diff(x[c(1, 2)])
  diff.y <- diff(y[c(1, 2)])

  x2 <- x + diff.x
  y2 <- y + diff.y

  grid$x1 <- x[grid$xi]
  grid$x2 <- x2[grid$xi]

  grid$y1 <- y[grid$yi]
  grid$y2 <- y2[grid$yi]

  grid$cluster <- paste(grid$xi, grid$yi, sep='-')
  return(grid)
}

# split up scseq parts so that can load what need (faster)
split_save_scseq <- function(scseq, dataset_dir) {

  # save as seperate parts
  dgc.logs <- SingleCellExperiment::logcounts(scseq)
  counts <- SingleCellExperiment::counts(scseq)
  t.logs <- Matrix::t(dgc.logs)

  SingleCellExperiment::logcounts(scseq) <- NULL
  SingleCellExperiment::counts(scseq) <- NULL

  qs::qsave(scseq, file.path(dataset_dir, 'shell.qs'), preset = 'fast')
  qs::qsave(dgc.logs, file.path(dataset_dir, 'dgclogs.qs'), preset = 'fast')
  HDF5Array::writeTENxMatrix(t.logs, file.path(dataset_dir, 'tlogs.tenx'), group = 'mm10')
  qs::qsave(counts, file.path(dataset_dir, 'counts.qs'), preset = 'fast')
}

# subset scseq to cells in test vs control contrast
subset_contrast <- function(scseq) {
  is.con <- scseq$orig.ident %in% c('test', 'ctrl')
  scseq <- scseq[, is.con]
  scseq$orig.ident <- droplevels(scseq$orig.ident)
  return(scseq)
}


# to allow pasting comma seperated genes into gene search
format_comma_regex <- function(regex) {
  regex <- gsub(', ', '$|^', regex)
  regex <- paste0('^', regex, '$')
  return(regex)
}


# get point with central x value cluster labels
get_label_coords <- function(coords, labels) {
  colnames(coords) <- c('x', 'y')
  coords$label <- labels

  # get point
  coords %>%
    dplyr::group_by(.data$label) %>%
    dplyr::summarise(
      x = stats::median(.data$x),
      y = stats::median(.data$y)) %>%
    dplyr::mutate(label = as.character(.data$label))
}


#' Take currently selected row and keep it in the same position
#'
#' Used to stop dataset change when adding new datasets
#'
#' @param datasets current data.frame of single-cell datasets
#' @param prev previous data.frame of single-cell datasets
#' @param curr name of currently selected row from \code{prev}
#'
#' @return \code{datasets} with \code{curr} at same row as it is in \code{prev}
#'
keep_curr_selected <- function(datasets, prev, curr) {

  # get currently selected row
  curr <- as.numeric(curr)
  curr_name <- prev[curr, 'name']

  # position in new datasets
  new_posn <- utils::tail(which(curr_name == datasets$name), 1)
  ndata <- nrow(datasets)

  # in case delete current that is last
  if (!length(new_posn) || curr > ndata) return(datasets)

  # move so that row at new_posn is at curr
  idx <- seq_len(nrow(datasets))
  idx_new <- replace(idx, c(curr, new_posn), c(new_posn, curr))
  datasets <- datasets[idx_new, ]
  datasets$value <- idx
  return(datasets)

}

h5_format_msg <- sprintf(
  paste0(
    'H5 file must be in feature barcode matrix format. ',
    '<a href="%s" target="_blank"><i class="fas fa-external-link-alt"></i></a>'
  ),
  'https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices')


validate_scseq_import <- function(up_df, samples) {

  msg <- NULL

  if (sum(is.na(samples))) {
    msg <- 'Specify sample for all files'
    return(msg)
  }

  uniq_samples <- unique(samples)
  for (sample in uniq_samples) {

    upi <- up_df[samples %in% sample,, drop = FALSE]
    files <- upi$name

    # matrix files
    mtx.file <- grep('.mtx', files, fixed = TRUE, value = TRUE)
    genes.file <- grep('features.tsv|genes.tsv', files, value = TRUE)
    barcodes.file <-  grep('barcodes.tsv', files, fixed = TRUE, value = TRUE)

    # if have matrix files, have exactly one of each type
    neach <- c(length(mtx.file), length(genes.file), length(barcodes.file))
    if (any(neach)) {
      if (!all(neach == 1) | length(files) != 3) {
        msg <- 'Need exactly one .mtx, genes.tsv, and barcodes.tsv files'
        return(msg)
      }
      next()
    }

    # no more than one H5 file per sample
    is_h5 <- grepl('[.]h5$|[.]hdf5$', files)
    if (sum(is_h5) > 1) {
      msg <- 'Specify only one sample per H5 file'
      return(msg)
    }

    # correct format for H5 file
    expect_h5 <- c('data', 'indices', 'indptr', 'shape', 'barcodes', 'gene_names', 'genes')
    if (sum(is_h5) == 1) {
      infile <- hdf5r::H5File$new(upi$datapath[is_h5], 'r')
      genomes <- names(infile)
      if ('matrix' %in% genomes) next()

      for (genome in genomes) {
        genome.names <- names(infile[[genome]])

        if (!all(expect_h5 %in% genome.names))
          return(h5_format_msg)
      }

      infile$close()
      next()
    }
  }

  return(msg)
}

validate_scseq_add_sample <- function(sample, rows) {
  msg <- NULL

  if (is.null(rows)) {
    msg <- 'No rows selected.'
    return(msg)
  }

  if (!shiny::isTruthy(sample)) {
    msg <- 'No sample name provided'
    return(msg)
  }

  msg <- validate_not_path(sample)
  return(msg)
}


#' Get HTML to make delete buttons in datatable rows
#'
#' @param session Shiny session object used for namespacing.
#' @param len Numeric number of rows
#' @param title Title for buttons. Default is \code{'Delete file'}
#'
#' @return Character vector with length \code{len} with delete buttons. Event
#' can be observed in a reactive by \code{input$delete_row} and the row clicked
#' will have have \code{'input$delete_row_i'} where 'i' is the row.
#' @export
#'
getDeleteRowButtons <- function(session, len, title='Delete file') {

  inputs <- character(len)
  for (i in seq_len(len)) {
    inputs[i] <- sprintf(
      '<span title="%s" class="btn dt-btn" id="%s" onclick="Shiny.onInputChange(\'%s\', this.id, {priority: \'event\'})"><icon class="%s"></icon></span>',
      title,
      paste0(session$ns('delete_'), i),
      session$ns('delete_row'),
      'far fa-trash-alt'
    )
  }
  inputs
}


# modal to integrate datasets
integrationModal <- function(session, choices) {
  ns <- session$ns

  modalDialog(
    withTags({
      div(id = ns('integration-form'),
          shinyWidgets::pickerInput(
            inputId = ns('integration_datasets'),
            label = 'Select datasets to integrate:',
            choices = choices,
            width = '100%',
            options = shinyWidgets::pickerOptions(
              `selected-text-format` = "count > 0",
              actionsBox = TRUE,
              liveSearch = TRUE,
              size=14),
            multiple = TRUE),
          tags$div(style='display: none', id=ns('integration_options_container'),
                   shinyWidgets::checkboxGroupButtons(
                     ns('integration_types'),
                     'Integration types:',
                     choices = c('harmony', 'fastMNN', 'reference'),
                     justified = TRUE,
                     selected = 'harmony',
                     checkIcon = list(
                       yes = icon("ok", lib = "glyphicon"),
                       no = icon("remove", lib = "glyphicon", style="color: transparent;"))
                   ),
                   div(id=ns('ref_name_container'), style='display: none;',
                       selectizeInput(
                         ns('ref_name'),
                         HTML('Select reference:'),
                         choices = refs, width = '100%')
                   ),
                   shinypanel::textInputWithValidation(
                     ns('integration_name'),
                     container_id = ns('name-container'),
                     label = 'Name for integrated dataset:',
                     help_id = ns('error_msg')
                   )
          ),
          div(tags$i(class = 'fas fa-exclamation-triangle', style='color: red;'), ' Need 3 or more samples for p-values and grid plots.', style='color: grey; font-style: italic;')
      )
    }),
    title = 'Integrate Single Cell Datasets',
    size = 'm',
    footer = tagList(
      actionButton(ns("submit_integration"), "Integrate Datasets", class = 'btn-warning'),
      tags$div(class='pull-left', modalButton('Cancel'))
    ),
    easyClose = FALSE
  )
}


# modal to upload single-cell dataset
importSingleCellModal <- function(session, show_init) {

  modalDialog(
    tags$div(
      class='alert alert-warning', role = 'alert',
      tags$div(tags$b("For each sample upload files:")),
      tags$br(),
      tags$div("- ", tags$code("filtered_feature_bc_matrix.h5"), "or"),
      tags$br(),
      tags$div("- ", tags$code("matrix.mtx"), ", ", tags$code("barcodes.tsv"), ", and", tags$code("features.tsv"), 'or'),
      tags$br(),
      tags$div("- ", tags$code("*.rds"), "or", tags$code("*.qs"), "with", tags$code("Seurat"), "or", tags$code("SingleCellExperiment"), "objects",
               tags$a(href = "https://docs.dseqr.com/docs/single-cell/add-dataset/#r-objects", target="_blank", "(requirements)")),
      hr(),
      '\U1F331 Add prefixes e.g.', tags$i(tags$b('sample_matrix.mtx')), ' to auto-name samples:',
      tags$a(href = 'https://dseqr.s3.amazonaws.com/GSM3972011_involved.zip', target = '_blank', 'example files.')
    ),
    div(class='upload-validation dashed-upload',
        attrib_replace(
          fileInput(
            placeholder = 'drag files here',
            session$ns('up_raw'),
            label = '',
            buttonLabel = 'Upload',
            width = '100%',
            accept = c('.h5', '.hdf5', '.tsv', '.fastq.gz', '.mtx', '.rds', '.qs'),
            multiple = TRUE
          ),
          list(id = session$ns("up_raw"), type = "file"),
          onchange = sprintf("checkSingleCellFileName(this, '%s');", session$ns("up_raw_errors"))
        ),
        tags$div(
          id = session$ns('validate-up-filetype'),
          tags$span(class = 'help-block', id = session$ns('error_msg_filetype'))
        )
    ),
    tags$div(
      id = session$ns('sample_name_container'),
      style = ifelse(show_init, '', 'display: none;'),
      hr(),
      shinypanel::textInputWithButtons(
        session$ns('sample_name'),
        'Sample name for selected rows:',
        actionButton(session$ns('add_sample'), '', icon('plus', class='fa-fw')),
        container_id = session$ns('validate-up'),
        help_id = session$ns('error_msg')
      ),
      hr()
    ),
    div(
      id = session$ns('up_table_container'),
      class= ifelse(show_init, 'dt-container', 'invisible-height dt-container'),
      DT::dataTableOutput(session$ns('up_table'), width = '100%'),
    ),
    title = 'Upload Single Cell Datasets',
    size = 'l',
    footer = tagList(
      actionButton(
        inputId = session$ns("import_datasets"),
        label = "Import Datasets",
        class = ifelse(show_init, 'btn-warning', 'btn-warning disabled')
      ),
      tags$div(class='pull-left', modalButton('Cancel'))
    ),
    easyClose = FALSE,
  )
}


# modal to delete dataset
deleteModal <- function(session, choices, type) {

  modalDialog(
    tags$div(class='selectize-fh',
             shinyWidgets::pickerInput(session$ns('remove_datasets'),
                                       label='Select datasets to delete:',
                                       width = '100%',
                                       choices = choices,
                                       options = shinyWidgets::pickerOptions(
                                         `selected-text-format` = "count > 0",
                                         actionsBox = TRUE,
                                         liveSearch = TRUE,
                                         width='fit',
                                         size=14),
                                       multiple = TRUE)
    ),
    tags$div(id=session$ns('confirm_delete_container'), style='display:none;',
             tags$div(class='alert alert-danger', 'This action cannot be undone.'),
             br(),
             textInput(session$ns('confirm_delete'), HTML('<span> Type <i><span style="color: gray">delete</span></i> to confirm:</span>'),
                       placeholder = "delete", width = '100%'
             ),
    ),
    title = paste('Delete', type, 'Datasets'),
    size = 'm',
    footer = tagList(
      actionButton(session$ns("delete_dataset"), "Delete Datasets"),
      tags$div(class='pull-left', modalButton("Cancel"))
    ),
    easyClose = FALSE,
  )
}


# modal to confirm adding single-cell dataset
confirmSubsetModal <- function(session,
                               new_dataset_name,
                               ref_name,
                               subset_metrics,
                               subset_clusters,
                               is_include) {

  metrics_ui <- ref_ui <- clusters_ui <- NULL
  direction <- ifelse(is_include, 'Include', 'Exclude')

  if (length(subset_metrics)) {
    metrics_ui <- tags$div(
      tags$br(),
      tags$div(tags$b(paste(direction, "cells with:"))),
      tags$div(lapply(subset_metrics, function(m) tags$div(tags$code(m))))
    )
  }

  if (length(subset_clusters)) {
    clusters_ui <- tags$div(
      tags$br(),
      tags$div(tags$b(paste(direction, "clusters:"))),
      tags$div(paste(subset_clusters, collapse = ', '))
    )
  }

  if (isTruthy(ref_name)) {
    is.ref <- refs$name == ref_name

    ref_ui <- tags$div(
      tags$br(),
      tags$div(tags$b("Reference dataset:")),
      tags$div(refs$label[is.ref])
    )

    new_dataset_name <- paste0(new_dataset_name, '_', refs$type[is.ref])
  }


  UI <- tags$div(
    class='alert alert-info', role = 'alert',
    tags$div(tags$b("New dataset name:")),
    tags$div(title = new_dataset_name, style = "text-overflow: ellipsis; overflow: hidden", new_dataset_name),
    ref_ui,
    clusters_ui,
    metrics_ui,
    hr(),
    '\U1F331 Click cancel to change settings'
  )

  modalDialog(
    UI,
    title = 'Create new single-cell dataset?',
    size = 'm',
    footer = tagList(
      actionButton(session$ns('confirm_subset'), 'Submit', class = 'btn-warning'),
      tags$div(class='pull-left', modalButton('Cancel'))
    )
  )
}


# modal to export dataset
exportModal <- function(session, choices, selected, options) {

  modalDialog(
    tags$div(
      class='alert alert-warning', role = 'alert',
      tags$div(tags$b("Download format:")),
      tags$br(),
      tags$div(tags$code("file.qs"), "with", tags$code("SingleCellExperiment"), "object."),
      hr(),
      tags$div(tags$b('\U1F331 To load:')),
      br(),
      tags$code("install.packages('qs')"),
      br(),
      tags$code("scseq <- qs::qread('file.qs')")
    ),
    selectizeInput(
      session$ns('export_dataset'),
      label = 'Select a single-cell dataset:',
      choices = c('', choices),
      options = options,
      selected = selected,
      width = '100%'),
    title = "Export Single Cell Dataset",
    size = 'm',
    footer = tagList(
      actionButton(session$ns("confirm_export"),
                   "Download Dataset",
                   class='btn-success'),
      tags$div(class='pull-left', modalButton("Cancel"))
    ),
    easyClose = TRUE
  )
}

#' shortest distance between cluster centroids
#'
#' TODO: use to order clusters by similarity
#'
#'
#' @param label_coords data.frame with columns label, x, and y giving coordinates
#' of cluster labels
#'
#' @return \code{label_coords} sorted by shortest path
#' @keywords internal
#' @importFrom TSP TSP solve_TSP
#'
sort_clusters <- function(label_coords) {

  m <- stats::dist(label_coords[, c('x', 'y')])

  # solve for shortest path
  tsp <- TSP::TSP(m)
  tour <- TSP::solve_TSP(tsp)

  # permutation vector for shortest tour
  perm_vec <- as.integer(tour)

  label_coords <- label_coords[perm_vec, ]
  return(label_coords)
}


#' Get expression colors for scatterplot
#'
#' @param ft.scaled expression values scaled from 0 to 1
#'
#' @return character vector of colors
#' @export
#'
#' @examples
#'
#' ft <- rnorm(100)
#' ft.scaled <- scales::rescale(ft)
#' colors <- get_expression_colors(ft.scaled)
#'
get_expression_colors <- function(ft.scaled) {
  cols <- const$colors$ft
  colors <- scales::gradient_n_pal(cols)(ft.scaled)

  # zero expression same as boolean off (white)
  colors[ft.scaled == 0] <- const$colors$bool[1]

  return(colors)
}
