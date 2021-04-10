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
    pred_annot <- make.unique(ref_preds, '_')

  } else if (ref_name == 'reset') {
    # reset annotation
    # load query annotation
    query_annot_path <- scseq_part_path(sc_dir, dataset_name, 'annot')
    query_annot <- qs::qread(query_annot_path)
    pred_annot <- as.character(seq_along(query_annot))

  } else {
    ref_annot_path <- scseq_part_path(sc_dir, ref_name, 'annot')
    ref_annot <- qs::qread(ref_annot_path)
    ref_annot <- gsub('_\\d+$', '', ref_annot)

    ref_preds <- ref_annot[as.numeric(ref_preds)]
    pred_annot <- make.unique(ref_preds, '_')
  }
  return(pred_annot)
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
  senv <- loadNamespace('SingleR')

  ref_names <- names(preds)
  dated <- c()

  dated <- sapply(ref_names, function(ref_name) {
    if (ref_name %in% ls(senv)) return(FALSE)
    resoln_name <- get_resoln_name(sc_dir, ref_name)

    ref_path <- scseq_part_path(sc_dir, resoln_name, 'scseq_sample')
    ref_date <- file.info(ref_path)$ctime
    ref_date <- as.character(ref_date)
    res_date <- names(preds[[ref_name]])[1]
    dated <- !identical(res_date, ref_date)

    return(dated)
  })

  return(preds[!dated])
}

#' Run differential abundance analysis
#'
#' @param scseq \code{SingleCellExperiment}
#'
#'
#' @keywords internal
diff_abundance <- function(scseq, annot, pairs = NULL) {

  abundances <- table(scseq$cluster, scseq$batch)
  abundances <- unclass(abundances)
  row.names(abundances) <- annot


  extra.info <- scseq@colData[match(colnames(abundances), scseq$batch),]
  y.ab <- edgeR::DGEList(abundances, samples=extra.info)

  adj <- edgeR::equalizeLibSizes(y.ab)$pseudo.counts
  adj <- round(adj)
  colnames(adj) <- paste0(colnames(adj), '_adjusted_counts')

  group <- y.ab$samples$orig.ident
  group <- stats::relevel(group, 'ctrl')
  y.ab$samples$group <- group

  # if only two samples, sort by logFC and return
  if (ncol(adj) == 2) {
    ord <- match(c('test', 'ctrl'), y.ab$samples$group)
    adj <- as.data.frame(adj[, ord])
    adj$logFC <- log2(adj[[1]]) - log2(adj[[2]])
    adj <- adj[order(abs(adj$logFC), decreasing = TRUE), ]
    return(adj)
  }


  if (!is.null(pairs))
    y.ab$samples$pair <- factor(pairs[colnames(y.ab),])


  keep <- edgeR::filterByExpr(y.ab, group=y.ab$samples$orig.ident)
  y.ab <- y.ab[keep,]

  # setup eset so that it will work with run_limma
  eset <- Biobase::ExpressionSet(y.ab$counts, Biobase::AnnotatedDataFrame(y.ab$samples))
  Biobase::assayDataElement(eset, 'vsd') <- Biobase::exprs(eset)
  Biobase::fData(eset)[, c('SYMBOL', 'ENTREZID')] <- row.names(eset)

  prev <- list(pdata = Biobase::pData(eset))
  eset <- crossmeta::run_limma_setup(eset, prev)
  lm_fit <- crossmeta::run_limma(eset)

  tt <- crossmeta::get_top_table(lm_fit, with.es = FALSE)
  tt <- cbind(tt, adj[row.names(tt), ])

  tt$ENTREZID <- NULL
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
get_label_transfer_choices <- function(anal_options, selected_anal, preds, species) {

  # omit previous and current
  is.recent <- anal_options$type == 'Previous Session'
  is.sel <- anal_options$name == selected_anal
  type   <- utils::tail(anal_options$type[is.sel], 1)
  type   <-  ifelse(type %in% c('Integrated', 'Individual'), selected_anal, type)

  anal_options <- anal_options[!is.sel & !is.recent, ]


  if (species == 'Homo sapiens') external <- 'Blueprint Encode Data'
  else if (species == 'Mus musculus') external <- 'Mouse RNAseq Data'

  choices <- data.frame(
    label = c('Reset Labels', external, anal_options$name),
    value = c('reset', gsub(' ', '', external), anal_options$name),
    optionLabel = c('Reset Labels', external, anal_options$optionLabel),
    type = factor(c(type, 'External Reference', anal_options$type), ordered = TRUE),
    stringsAsFactors = FALSE
  )


  choices$preds <- choices$value %in% names(preds)
  choices <- rbind(NA, choices)

  return(choices)
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

#' Get choices for included cluster in integration
#' @param dataset_names Names of analyses selected for integration
#' @param anal_colors Character vector of colors to indicate analysis
#' @param data_dir Directory with single cell analyses.
#'
#' @return data.frame with columns for rendering selectizeInput include choices
#'
#' @keywords internal
get_exclude_choices <- function(dataset_names, data_dir, anal_colors = NA) {

  if (is.null(dataset_names)) return(NULL)

  # load markers and annotation for each
  resolns <- load_resoln(file.path(data_dir, dataset_names))
  resoln_names <- file.path(dataset_names, resolns)

  annot_paths <- scseq_part_path(data_dir, resoln_names, 'annot')
  marker_paths <- scseq_part_path(data_dir, resoln_names, 'markers')

  annots <- lapply(annot_paths, qs::qread)
  clusters <- lapply(annots, seq_along)

  exclude_choices <- lapply(seq_along(dataset_names), function(i) {
    data.frame(
      name = stringr::str_trunc(annots[[i]], 27),
      value = paste(dataset_names[i], clusters[[i]], sep = '_'),
      anal = dataset_names[i],
      label = annots[[i]],
      color = anal_colors[i], stringsAsFactors = FALSE
    )
  })

  do.call(rbind, exclude_choices)
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

    choices <- rbind(tail(choices, 1), head(choices, nrow(choices)-1))

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

check_has_scseq <- function(dataset_names, sc_dir) {
  sapply(dataset_names, function(dataset_name) {
    fnames <- list.files(file.path(sc_dir, dataset_name))
    any(fnames %in% c('scseq.qs', 'scseq.qs'))
  })
}


#' Get data.frame of single-cell dataset choices
#'
#' @param sc_dir Directory containing single-cell datasets
#'
#' @return data.frame of single-cell dataset choices for selectizeInput
#'
#' @keywords internal
get_sc_dataset_choices <- function(sc_dir) {


  # exclude missing from integrated (e.g. manual delete)
  integrated <- qread.safe(file.path(sc_dir, 'integrated.qs'))
  has.scseq <- check_has_scseq(integrated, sc_dir)
  integrated <- integrated[has.scseq]

  int_type <- lapply(integrated, function(int) qread.safe(file.path(sc_dir, int, 'founder.qs'),
                                                          .nofile = 'Integrated',
                                                          .nullfile = 'Integrated'))
  int_type <- unlist(int_type)

  sub <- duplicated(int_type) | duplicated(int_type, fromLast = TRUE)
  int_type[!sub] <- 'Integrated'
  int_opt <- integrated
  int_opt[sub] <- stringr::str_replace(int_opt[sub], paste0(int_type[sub], '_'), '')

  individual <- setdiff(list.files(sc_dir), c(integrated, 'integrated.qs'))

  # exclude individual without scseq (e.g. folder with fastq.gz files only)
  has.scseq <- check_has_scseq(individual, sc_dir)
  individual <- individual[unlist(has.scseq)]

  # founder for subsets as type
  ind_type <- sapply(individual, function(ind) qread.safe(file.path(sc_dir, ind, 'founder.qs'),
                                                          .nofile = 'Individual',
                                                          .nullfile = 'Individual'), USE.NAMES = FALSE)


  # add founders that are currently individual
  # from saving without QC0 QC1
  self.founder <- individual %in% ind_type
  ind_type[self.founder] <- individual[self.founder]

  # exclude founder name from option label
  sub <- ind_type != 'Individual'
  ind_opt <- individual
  ind_opt[sub] <- stringr::str_replace(ind_opt[sub], paste0(ind_type[sub], '_?'), '')

  label <- c(integrated, individual)
  type <- c(int_type, ind_type)
  opt <- c(int_opt, ind_opt)

  # get previously selected
  if (length(individual)) {
    prev <- qread.safe(file.path(sc_dir, 'prev_dataset.qs'), .nullfile = individual[1])
    if (!prev %in% label) prev <- individual[1]

    prev_type <- type[label == prev]
    founder <- get_founder(sc_dir, prev)
    if (prev_type == founder)
      prev <- label[type %in% founder]

    type <- c(rep('Previous Session', length(prev)), type)
    label <- c(prev, label)
    opt <- c(prev, opt)
  }

  choices <- data.frame(value = seq_along(label),
                        name = label,
                        label = label,
                        type = type,
                        optionLabel = opt,
                        stringsAsFactors = FALSE)

  return(choices)
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
  res <- lapply(unique(types), function(type) c('', res[types==type]))
  names(res) <- unique(types)
  return(res)

}


#' Get single-cell directories without fastq or cellranger files
#'
#' Used so that all directories shown for new single-cell dataset contain
#' raw files.
#'
#' @param sc_dir Path to single-cell data directory
#'
#' @return Character vector of directories/files to hide
#' @keywords internal
#'
get_exclude_dirs <- function(sc_dir) {
  dirs <- list.dirs(sc_dir, full.names = FALSE, recursive = FALSE)
  exclude <- c()
  for (dir in dirs) {
    full_dir <- file.path(sc_dir, dir)
    files <- list.files(full_dir)
    if (!any(grepl('fastq.gz$|.mtx$|.h5$', files))) {
      exclude <- c(exclude, dir)
    }
  }

  exclude <- c(exclude, 'integrated.qs', 'prev_dataset.qs')
  return(exclude)
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

  ft <- strsplit(metric, '[|><=\\)\\(&!*]')[[1]]
  ft <- gsub(' ', '', ft)
  ft <- ft[ft != '']
  not.num <- is.na(suppressWarnings(as.numeric(ft)))
  ft[not.num]
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
  expr <- expr[row.names(expr) %in% ft,, drop = FALSE]
  expr <- t(as.matrix(expr))

  qcs <- scseq@colData
  qcs <- qcs[, colnames(qcs) %in% ft, drop = FALSE]

  dat <- cbind(expr, qcs)
  dat <- as.data.frame(cbind(expr, qcs))
  colnames(dat) <- c(colnames(expr), colnames(qcs))
  dat <- dat[, match(ft, colnames(dat)), drop = FALSE]

  # convert features to generic (in case e.g. minus)
  seq <- 1:ncol(dat)
  ft.num <- paste0('ft', seq)
  colnames(dat) <- ft.num

  metric.num <- metric
  for (i in seq) metric.num <- gsub(paste0('\\b', ft[i], '\\b'), ft.num[i], metric.num)

  dat <- within(dat, metric <- eval(rlang::parse_expr(metric.num)))
  dat <- dat[, 'metric', drop = FALSE]
  colnames(dat) <- metric
  return(dat)
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
  if (is.null(x)) return(NULL)
  x <- format(as.character(x), justify = justify)
  gsub(' ', '&nbsp;&nbsp;', x)
}

#' Get/Save cluster stats for single-cell related selectizeInputs
#'
#' @param resoln_dir Sub directory with single cell dataset info specific to resolution.
#' @param scseq \code{SingleCellExperiment} object to get/save stats for.
#'   if \code{NULL} (Default), will be loaded.
#' @param top_tables List of \code{limma::topTable} results used by
#'   scSampleComparison for integrated datasets.
#' @param has_replicates Boolean indicating if integrated dataset has multiple
#'   samples in at least one group. Used by scSampleComparison.
#' @param use_disk Should results by saved to disk? used by scSampleComparison
#'   to persist results to disk.
#' @inheritParams get_cluster_choices
#'
#' @return List with cluster stats
get_cluster_stats <- function(resoln_dir = NULL, scseq = NULL, top_tables = NULL, has_replicates = FALSE, use_disk = FALSE, sample_comparison = FALSE) {

  stats_path <- file.path(resoln_dir, 'cluster_stats.qs')
  if (file.exists(stats_path) && use_disk) return(qs::qread(stats_path))

  if (is.null(scseq)) {
    dataset_dir <- dirname(resoln_dir)
    scseq <- load_scseq(dataset_dir, default_clusters = FALSE)
    scseq <- attach_clusters(scseq, resoln_dir)
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
  if (sample_comparison & has_replicates) {
    nsig <- rep(0, nbins)
    names(nsig) <- seq_len(nbins)

    test_clusters <- names(top_tables)
    nsig[test_clusters] <- sapply(top_tables, function(tt) {sum(tt$adj.P.Val.Amb < 0.05 & !tt$ambient)})
    stats$nsig <- nsig
  }

  # show number of non-ambient with logFC > 1
  if (sample_comparison) {
    nbig <- rep(0, nbins)
    names(nbig) <- seq_len(nbins)

    test_clusters <- names(top_tables)
    nbig[test_clusters] <- sapply(top_tables, function(tt) {sum(abs(tt$logFC) > 1 & !tt$ambient)})
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
get_contrast_markers <- function(con, dataset_dir) {
  con <- strsplit(con, '-vs-')[[1]]
  tests_dir <- file.path(dataset_dir, 'tests')
  pairs_path <- file.path(tests_dir, 'pairs.qs')
  pairs <- qs::qread(pairs_path)

  keep <- apply(pairs, 1, function(row) all(con %in% row))
  keep <- which(keep)

  stat_paths <- file.path(tests_dir, paste0('statistics_pair', keep, '.qs'))
  tests <- list(pairs = pairs[keep, ],
                statistics = lapply(stat_paths, qs::qread))

  # returns both directions
  con_markers <- get_scseq_markers(tests)
  names(con_markers) <- paste0(names(con_markers), '-vs-', rev(names(con_markers)))

  return(con_markers)
}


#' Add cell percents to gene choices for single cell
#'
#' @param scseq \code{SingleCellExperiment} object
#' @param markers data.frame of marker genes.
#'
#' @return data.frame of all genes, with markers on top and cell percent columns
#'
#' @keywords internal
get_gene_choices <- function(markers,
                             qc_metrics = NULL,
                             type = NULL,
                             qc_first = FALSE,
                             species = 'Homo sapiens',
                             tx2gene = NULL) {
  markers <- row.names(markers)

  qc_type <- ifelse(names(qc_metrics) == 'numeric', 'QC Score', 'Boolean Features')
  gene_type <- rep('Gene', length(markers))

  if (qc_first) {
    choices <- c(qc_metrics, markers)
    type <- c(qc_type, gene_type)

  } else {
    choices <- c(markers, qc_metrics)
    type <- c(gene_type, qc_type)
  }

  if (!grepl('sapiens|musculus', species))
    stop('Only Homo sapiens and Mus musculus supported')

  if (is.null(tx2gene)) tx2gene <- dseqr.data::load_tx2gene(species)

  idx <- match(choices, tx2gene$gene_name)
  desc <- tx2gene$description[idx]
  desc[is.na(desc)] <- choices[is.na(desc)]

  choices <- data.table::data.table(label = choices,
                                    value = choices,
                                    description = desc,
                                    type = type,
                                    stringsAsFactors = FALSE)

  choices




  return(choices)
}


#' Integrate previously saved SingleCellExperiments
#'
#' Performs integration and saves as a new analysis.
#' Used by \code{explore_scseq_clusters} shiny app.
#'
#' @param sc_dir Directory with saved single-cell datasets.
#' @param test Character vector of test analysis names.
#' @param ctrl Character vector of control analysis names.
#' @param exclude_clusters Charactor vector of clusters for excluding cells.
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
  test,
  ctrl,
  integration_name,
  exclude_clusters,
  exclude_cells = NULL,
  integration_type = c('harmony', 'fastMNN'),
  subset_metrics = NULL,
  is_include = NULL,
  founder = integration_name,
  pairs = NULL,
  hvgs = NULL,
  azimuth_ref = NULL,
  npcs = 30,
  cluster_alg = 'leiden',
  resoln = 1,
  progress = NULL,
  value = 0) {

  # for save_scseq_args
  args <- c(as.list(environment()))
  args$progress <- args$sc_dir <- NULL
  args$date <- Sys.time()
  args$value <- NULL

  if (is.null(progress)) {
    progress <- list(set = function(value, message = '', detail = '') {
      cat(value, message, detail, '...\n')
    })
  }

  dataset_name <- paste(integration_name, integration_type, sep = '_')

  # save dummy data if testing shiny
  if (isTRUE(getOption('shiny.testmode'))) {
    scseq_data <- list(scseq = NULL, markers = NULL, annot = NULL)
    save_scseq_data(scseq_data, dataset_name, sc_dir, add_integrated = TRUE)
    return(TRUE)
  }

  progress$set(value+1, detail = 'loading')
  test_scseqs <- load_scseq_subsets(test,
                                    sc_dir,
                                    exclude_clusters,
                                    subset_metrics,
                                    is_include,
                                    'test',
                                    exclude_cells = exclude_cells)

  if (!length(test_scseqs)) {
    progress$set(value+1, detail = 'error: no test cells')
    Sys.sleep(3)
    return(FALSE)
  }

  ctrl_scseqs <- load_scseq_subsets(ctrl,
                                    sc_dir,
                                    exclude_clusters,
                                    subset_metrics,
                                    is_include,
                                    'ctrl',
                                    exclude_cells = exclude_cells)

  if (!length(ctrl_scseqs)) {
    progress$set(value+1, detail = 'error: no control cells')
    Sys.sleep(3)
    return(FALSE)
  }

  # preserve identity of original samples and integrate
  scseqs <- c(test_scseqs, ctrl_scseqs)
  ambient <- get_integrated_ambient(scseqs)

  # make sure all the same species
  species <- unique(sapply(scseqs, function(x) x@metadata$species))
  if(length(species) > 1) stop('Multi-species integration not supported.')

  progress$set(value+2, detail = 'integrating')
  combined <- integrate_scseqs(scseqs, type = integration_type, pairs = pairs, hvgs = hvgs, azimuth_ref = azimuth_ref)
  combined$project <- dataset_name

  # retain original QC metrics
  combined <- add_combined_metrics(combined, scseqs)

  rm(scseqs, test_scseqs, ctrl_scseqs); gc()

  # add ambient outlier info
  combined <- add_integrated_ambient(combined, ambient)

  combined@metadata$species <- species
  combined@metadata$npcs <- npcs

  # save what is stable with resolution change
  scseq_data <- list(scseq = combined,
                     species = species,
                     ambient = ambient,
                     pairs = pairs,
                     founder = founder,
                     resoln = resoln)

  is_azimuth <- integration_type == 'Azimuth'
  if (is_azimuth) {
    resoln <- switch(azimuth_ref, 'human_pbmc' = 2)
    scseq_data$resoln <- resoln
    scseq_data$azimuth_ref <- azimuth_ref

    save_scseq_data(scseq_data, dataset_name, sc_dir, add_integrated = TRUE)
    save_azimuth_clusters(combined@colData, dataset_name, sc_dir)

  } else {
    # TSNE on corrected reducedDim
    progress$set(value+3, detail = 'reducing')
    combined <- run_reduction(combined, dimred = 'corrected')

    # get graph
    scseq_data$snn_graph <- get_snn_graph(combined, npcs = npcs)
    combined$cluster <- get_clusters(scseq_data$snn_graph, cluster_alg, resoln)

    save_scseq_data(scseq_data, dataset_name, sc_dir, add_integrated = TRUE)
  }

  progress$set(value+4, detail = 'saving loom')
  save_scle(combined, file.path(sc_dir, dataset_name))

  # run things that can change with resolution change
  run_post_cluster(combined,
                   dataset_name,
                   sc_dir,
                   resoln,
                   progress,
                   value+4,
                   reset_annot = !is_azimuth)

  save_scseq_args(args, dataset_name, sc_dir)

  return(TRUE)
}

#' Post-Cluster Processing of SingleCellExperiment
#'
#' @param scseq SingleCellExperiment
#' @param resoln resolution parameter
#' @param integrated is \code{scseq} integrated?
#' @param progress progress
#' @param value value
#'
#' @return NULL
#'
run_post_cluster <- function(scseq, dataset_name, sc_dir, resoln, progress = NULL, value = 0, reset_annot = TRUE) {

  integrated <- 'corrected' %in% SingleCellExperiment::reducedDimNames(scseq)
  species <- scseq@metadata$species

  if (is.null(progress)) {
    progress <- list(set = function(value, message = '', detail = '') {
      cat(value, message, detail, '...\n')
    })
  }

  progress$set(value, detail = 'getting markers')
  markers <- get_presto_markers(scseq)


  # used for label transfer
  scseq_sample <- downsample_clusters(scseq)


  anal <- list(scseq_sample = scseq_sample,
               markers = markers,
               clusters = scseq$cluster,
               applied = TRUE)

  if (reset_annot) anal$annot <- names(markers)

  if (integrated) {

    progress$set(value+1, detail = 'pseudobulking')
    summed <- aggregate_across_cells(scseq)

    obj <- scseq
    pbulk_esets <- NULL
    has_replicates <- length(unique(scseq$batch)) > 2

    if (has_replicates) {
      pairs_path <- scseq_part_path(sc_dir, dataset_name, 'pairs')
      pairs <- qread.safe(pairs_path)

      if (species == 'Homo sapiens') release <- '94'
      else if (species == 'Mus musculus') release <- '98'

      pbulk_esets <- obj <- construct_pbulk_esets(summed, pairs, species, release)
    }

    progress$set(value+2, detail = 'fitting models')
    lm_fit <- run_limma_scseq(obj)
    anal_int <- list(summed = summed,
                     lm_fit_0svs = lm_fit,
                     pbulk_esets = pbulk_esets,
                     has_replicates = has_replicates)

    anal <- c(anal, anal_int)
  }

  # save in subdirectory e.g. snn0.8
  progress$set(value+3, detail = 'saving')
  dataset_subname <- file.path(dataset_name, paste0('snn', resoln))
  save_scseq_data(anal, dataset_subname, sc_dir, overwrite = FALSE)

}

aggregate_across_cells <- function(scseq) {
  counts_mat <- SingleCellExperiment::counts(scseq)
  batch <- scseq$batch
  cluster <- scseq$cluster
  orig.ident <- scseq$orig.ident
  y <- paste(cluster, batch, sep='_')
  levs <- unique(y)

  agg <- t(presto::sumGroups(counts_mat, factor(y, levs), 1))
  row.names(agg) <- row.names(counts_mat)
  colnames(agg) <- levs

  # setup colData
  dup.y <- duplicated(y)
  cdata <- S4Vectors::DataFrame(
    cluster = cluster[!dup.y],
    batch = batch[!dup.y],
    orig.ident = orig.ident[!dup.y])

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = agg),
    rowData = SingleCellExperiment::rowData(scseq),
    colData = cdata
  )

  return(sce)
}




#' Integration Utility Function for Background Process
#'
#' Used as background process via \link{callr::r_bg}
#'
#' @inheritParams integrate_saved_scseqs
#' @param integration_types Character vector of integration types to run.
#'
#' @return TRUE is successful, otherwise FALSE
#' @keywords internal
#'
run_integrate_saved_scseqs <- function(sc_dir, test, ctrl, integration_name,
                                       integration_types, exclude_clusters, pairs, azimuth_ref) {

  for (i in seq_along(integration_types)) {
    type <- integration_types[i]
    ref <- NULL
    if (type == 'Azimuth') ref <- azimuth_ref

    # run integration
    res <- integrate_saved_scseqs(sc_dir,
                                  test = test,
                                  ctrl = ctrl,
                                  integration_name = integration_name,
                                  integration_type = type,
                                  exclude_clusters = exclude_clusters,
                                  pairs = pairs,
                                  azimuth_ref = ref,
                                  value = i*8-8)

    # stop subsequent integration types if error
    if (!res) return(FALSE)
  }
  return(TRUE)
}


#' Subset previously saved SingleCellExperiment
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
                               exclude_clusters,
                               exclude_cells,
                               subset_metrics,
                               is_integrated,
                               is_include,
                               progress = NULL,
                               hvgs = NULL,
                               azimuth_ref = NULL) {

  if (is.null(progress)) {
    progress <- list(set = function(value, message = '', detail = '') {
      cat(value, message, detail, '...\n')
    })
  }
  progress$set(1, detail = 'loading')

  is_include <- rep(is_include, length(subset_metrics))

  if (is_integrated) {

    # add to previous subsets
    args <- load_args(sc_dir, from_dataset)
    args$is_include <- c(args$is_include, is_include)
    args$subset_metrics <- c(args$subset_metrics, subset_metrics)

    for (ds in names(args$exclude_cells))
      args$exclude_cells[[ds]] <- c(args$exclude_cells[[ds]], exclude_cells[[ds]])

    res <- integrate_saved_scseqs(sc_dir,
                                  test = args$test,
                                  ctrl = args$ctrl,
                                  integration_name = dataset_name,
                                  integration_type = args$integration_type,
                                  exclude_clusters = args$exclude_clusters,
                                  exclude_cells =  args$exclude_cells,
                                  subset_metrics = args$subset_metrics,
                                  is_include = args$is_include,
                                  founder = founder,
                                  pairs = args$pairs,
                                  hvgs = hvgs,
                                  azimuth_ref = azimuth_ref,
                                  progress = progress,
                                  value = 1)
    return(res)

  } else {

    # for save_scseq_args
    args <- c(as.list(environment()))
    args$progress <- args$sc_dir <- NULL
    args$date <- Sys.time()

    scseq <- load_scseq_subsets(from_dataset,
                                sc_dir,
                                exclude_clusters,
                                subset_metrics,
                                is_include,
                                ident=from_dataset)

    if (!length(scseq)) {
      progress$set(1, detail = "error: no cells")
      Sys.sleep(3)
      return(FALSE)
    }

    process_raw_scseq(scseq[[1]],
                      dataset_name,
                      sc_dir,
                      hvgs = hvgs,
                      progress = progress,
                      value = 1,
                      founder = founder,
                      azimuth_ref = azimuth_ref)

    save_scseq_args(args, dataset_name, sc_dir)
    return(TRUE)
  }
}

load_args <- function(sc_dir, dataset_name) {
  jsonlite::read_json(file.path(sc_dir, dataset_name, 'args.json'), simplifyVector = TRUE)
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
  founder <- qs::qread(fpath)
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

  # add mixing (Control, Test, sample origin) related metrics
  combined$is_test <- combined$orig.ident == 'test'
  combined$is_ctrl <- !combined$is_test

  samples <- unique(combined$batch)
  if (length(samples > 2)) {
    for (sample in samples)
      combined[[sample]] <- combined$batch == sample
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
load_scseq_subsets <- function(dataset_names, sc_dir, exclude_clusters, subset_metrics = NULL, is_include = NULL, ident = 'test', exclude_cells = NULL) {

  exclude_datasets <- gsub('^(.+?)_\\d+$', '\\1', exclude_clusters)
  exclude_clusters <- gsub('^.+?_(\\d+)$', '\\1', exclude_clusters)

  scseqs <- list()

  # load each scseq and exclude based on metrics
  for (dataset_name in dataset_names) {
    dataset_dir <- file.path(sc_dir, dataset_name)
    scseq <- load_scseq_qs(dataset_dir)

    # set orig.resoln to track resolution of origin datasets
    resoln_name <- load_resoln(dataset_dir)
    scseq$orig.resoln <- resoln_name

    # set orig.ident to ctrl/test (integration) or original dataset name (subset)
    scseq$orig.ident <- factor(ident)

    if (length(subset_metrics)) {
      cdata <- scseq@colData
      metrics <- lapply(subset_metrics, evaluate_custom_metric, scseq)
      metrics <- do.call(cbind, metrics)
      if (!is.null(metrics)) cdata <- cbind(cdata, metrics)

      exclude <- cdata[, subset_metrics, drop = FALSE]

      for (i in seq_len(ncol(exclude)))
        if (is_include[i]) exclude[,i] <- !exclude[,i]

      exclude <- apply(exclude, 1, any)
      scseq <- scseq[, !exclude]
      gc()

      # remove subset_metrics hardcoded in scseq
      # no/all cells will meet them by definition
      scseq@colData <- scseq@colData[, !colnames(scseq@colData) %in% subset_metrics, drop = FALSE]
    }

    # remove excluded clusters
    is.exclude <- exclude_datasets == dataset_name
    if (any(is.exclude)) {
      exclude <- exclude_clusters[is.exclude]
      scseq <- scseq[, !scseq$cluster %in% exclude]
      gc()
    }

    # remove excluded cells
    cells <- exclude_cells[[dataset_name]]
    if (!is.null(cells)) {
      exclude <- colnames(scseq) %in% cells
      scseq <- scseq[, !exclude]

      if (sum(exclude) != length(cells))
        warning("Exclude cells didn't match perfectly.")
    }

    # require that have cells left
    if (ncol(scseq)>0) scseqs[[dataset_name]] <- scseq
  }

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
#' @param add_integrated Add analysis to integrated.qs file? Default is \code{FALSE}
#'
#' @return NULL
#' @keywords internal
save_scseq_data <- function(scseq_data, dataset_name, sc_dir, add_integrated = FALSE, overwrite = TRUE) {
  dataset_dir <- file.path(sc_dir, dataset_name)

  if (add_integrated) {
    # add to integrated if new
    int_path <- file.path(sc_dir, 'integrated.qs')
    int_options <- c(qs::qread(int_path), dataset_name)
    qs::qsave(unique(int_options), int_path)
  }

  # remove all previous data in case overwriting
  if (overwrite) unlink(dataset_dir, recursive = TRUE)

  dir.create(dataset_dir, showWarnings = FALSE)
  for (type in names(scseq_data)) {
    item <- scseq_data[[type]]

    if (type == 'markers') {
      # save marker data.frames individually for fast loading

      for (i in names(item))
        qs::qsave(item[[i]], scseq_part_path(sc_dir, dataset_name, paste0('markers_', i)))

    } else if (type == 'tests') {
      # save pairwise test statistics for fast single group comparisons

      tests_dir <- file.path(dataset_name, 'tests')
      dir.create(file.path(sc_dir, tests_dir))

      qs::qsave(item$pairs, scseq_part_path(sc_dir, tests_dir, 'pairs'))

      for (i in seq_along(item$statistics))
        qs::qsave(item$statistics[[i]], scseq_part_path(sc_dir, tests_dir, paste0('statistics_pair', i)))

    } else {
      qs::qsave(item, scseq_part_path(sc_dir, dataset_name, type))
    }

  }

  return(NULL)
}

#' Load SingleCellExperiment from qs file
#'
#' Also attaches clusters from last applied leiden resolution and stores resolution.
#'
#' @param dataset_dir Path to folder with scseq.qs file
#'
#' @return SingleCellExperiment
#' @keywords internal
#'
load_scseq_qs <- function(dataset_dir) {
  scseq_path <- file.path(dataset_dir, 'scseq.qs')
  resoln_name <- load_resoln(dataset_dir)

  scseq <- qs::qread(scseq_path)
  scseq <- attach_clusters(scseq, file.path(dataset_dir, resoln_name))

  return(scseq)
}

#' Validate dataset selection for integration
#'
#' @param test Character vector of test dataset names
#' @param ctrl Character vector of control dataset names
#'
#' @return \code{NULL} if valid, otherwise an error message
#'
#' @keywords internal
validate_integration <- function(test, ctrl, pairs) {
  msg <- NULL

  if (!is.null(pairs)) {

    if (!all(pairs$sample %in% c(test, ctrl)))
      msg <- 'Samples missing from pairs csv'
  }

  return(msg)
}

validate_subset <- function(from_dataset, subset_name, subset_clusters, is_include, hvgs) {
  msg <- NULL

  have.dataset <- isTruthy(from_dataset)
  have.name <- isTruthy(subset_name)
  all.excluded <- is_include & !isTruthy(subset_clusters)
  all.included <- !isTruthy(subset_clusters) & !is_include
  no.hvgs <- !isTruthy(hvgs)

  if (!have.dataset) {
    msg <- 'No dataset selected'

  } else if (!have.name) {
    msg <- 'No name given'

  } else if (all.excluded) {
    msg <- 'All excluded'

  } else if (all.included & no.hvgs) {
    msg <- 'All included without custom genes'
  }

  return(msg)
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
#' @param ambient Character vector of ambient genes to exclude from drug queries.
#'
#' @return \code{res} with drug query results added to \code{'cmap'} \code{'l1000'} slots.
#'
#' @keywords internal
run_drug_queries <- function(top_table, drug_paths, es, ambient = NULL, species = NULL, ngenes = 200) {

  # get dprime effect size values for analysis
  dprimes <- get_dprimes(top_table)

  if (isTRUE(species == 'Mus musculus')) {
    names(dprimes) <- toupper(names(dprimes))
    ambient <- toupper(ambient)
  }

  # exclude ambient (for single cell only)
  dprimes <- dprimes[!names(dprimes) %in% ambient]

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
#' dseqr:::get_nearest_row(truth, test)
#' #[1] 4 3
get_nearest_row <- function(truth, test) {
  diffs   <- truth[rep(1:nrow(truth), nrow(test)),] -test[rep(1:nrow(test), each=nrow(truth)),]
  eucdiff <- function(x) sqrt(rowSums(x^2))
  max.col(-matrix(eucdiff(diffs), nrow=nrow(test), byrow=TRUE), "first")
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
      last <- tail(msgs, 1)
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
  ref_resoln_name <- file.path(dataset_name, paste0('snn', prev_resoln))
  ref_cluster <- qs::qread(file.path(sc_dir, ref_resoln_name, 'clusters.qs'))

  # query clusters are new resolution
  query_resoln_name <- file.path(dataset_name, paste0('snn', resoln))
  query_cluster <- qs::qread(file.path(sc_dir, query_resoln_name, 'clusters.qs'))

  # transfer labels
  tab <- table(assigned = ref_cluster, cluster = query_cluster)
  pred <- row.names(tab)[apply(tab, 2, which.max)]
  annot <- get_pred_annot(pred, ref_resoln_name, query_resoln_name, sc_dir)
  annot_nums <- as.character(seq_along(annot))

  # keep ordered nums where prediction is numeric
  suppressWarnings(is.num <- !is.na(as.numeric(gsub('_\\d+$', '', annot))))
  annot[is.num] <- annot_nums[is.num]

  qs::qsave(annot, file.path(sc_dir, query_resoln_name, 'annot.qs'))
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
  paste0('snn', sapply(resoln_path, qread.safe, .nofile = 1))
}
