#' Get predicted annotation for label transfer
#'
#' Clusters with an average prediction scores below \code{min.score} retain their original labels.
#'
#' @param ref_preds data.frame generated in \code{\link{labelTransferForm}} on event \code{submit_transfer}
#' @param ref_name Name of reference analysis that labels are transfered from.
#' @param dataset_name Name of analysis that labels are transfered to.
#' @param sc_dir Directory containing folders with analyses for \code{ref_name} and \code{dataset_name}.

#' @return Character vector of predicted labels from \code{ref_name}.
#' @export
#' @keywords internal
get_pred_annot <- function(ref_preds, ref_name, dataset_name, sc_dir) {


  # load query annotation
  query_annot_path <- scseq_part_path(sc_dir, dataset_name, 'annot')
  query_annot <- readRDS(query_annot_path)

  senv <- loadNamespace('SingleR')

  if (ref_name %in% ls(senv)) {
    pred_annot <- make.unique(ref_preds, '_')

  } else if (ref_name == 'reset') {
    # reset annotation
    pred_annot <- as.character(seq_along(query_annot))

  } else {
    ref_annot_path <- scseq_part_path(sc_dir, ref_name, 'annot')
    ref_annot <- readRDS(ref_annot_path)

    ref_preds <- ref_annot[as.numeric(ref_preds)]
    pred_annot <- make.unique(ref_preds, '_')
  }
  return(pred_annot)
}

#' Run differential abundance analysis
#'
#' @param scseq \code{SingleCellExperiment}
#'
#' @export
#' @keywords internal
diff_abundance <- function(scseq, annot, pairs = NULL) {

  abundances <- table(scseq$cluster, scseq$batch)
  abundances <- unclass(abundances)
  row.names(abundances) <- annot

  if (ncol(abundances) == 2) return(abundances)

  extra.info <- scseq@colData[match(colnames(abundances), scseq$batch),]
  y.ab <- edgeR::DGEList(abundances, samples=extra.info)

  if (!is.null(pairs))
    y.ab$samples$pair <- factor(pairs[colnames(y.ab),])

  group <- y.ab$samples$orig.ident
  group <- relevel(group, 'ctrl')
  y.ab$samples$group <- group

  keep <- edgeR::filterByExpr(y.ab, group=y.ab$samples$orig.ident)
  y.ab <- y.ab[keep,]

  # setup eset so that it will work with run_limma
  eset <- Biobase::ExpressionSet(y.ab$counts, Biobase::AnnotatedDataFrame(y.ab$samples))
  Biobase::assayDataElement(eset, 'vsd') <- Biobase::exprs(eset)
  Biobase::fData(eset)[, c('SYMBOL', 'ENTREZID')] <- row.names(eset)

  lm_fit <- run_limma(eset, prev_anal = list(pdata = Biobase::pData(eset)))

  tt <- get_top_table(lm_fit, with.es = FALSE)
  return(tt)
}

run_varPart <- function(eset) {

  pdata <- Biobase::pData(eset)
  pdata$pair <- factor(pdata$pair)
  lib.size <- pdata$lib.size * pdata$norm.factors

  y <- Biobase::exprs(eset)
  v <- limma::voom(y, lm_fit$mod, lib.size = lib.size)
  form <- ~ (1|group) + (1|pair)
  param <- BiocParallel::SerialParam()
  varPart <- variancePartition::fitExtractVarPartModel(v, form, info, BPPARAM = param)
  return(varPart)
}

run_dream <- function(eset) {

  pdata <- Biobase::pData(eset)
  pdata$pair <- factor(pdata$pair)
  lib.size <- pdata$lib.size * pdata$norm.factors

  form <- ~group + (1|pair)
  y <- Biobase::exprs(eset)
  param <- BiocParallel::SerialParam()
  v <- variancePartition::voomWithDreamWeights(y, form, pdata, lib.size = lib.size, BPPARAM = param)
  fitmm <- variancePartition::dream(v, form, pdata, BPPARAM = param)
  tt <- limma::topTable(fitmm, coef = 'grouptest', number = Inf, sort.by = 'P')
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
#' @export
#' @keywords internal
get_label_transfer_choices <- function(anal_options, selected_anal, preds, species) {

  # determine if is subset
  is.recent <- anal_options$type == 'Previous Session'
  is.sel <- anal_options$name == selected_anal
  type   <- tail(anal_options$type[is.sel], 1)
  type   <-  ifelse(type %in% c('Integrated', 'Individual'), selected_anal, type)

  anal_options <- anal_options[!is.sel & !is.recent, ]


  if (species == 'Homo sapiens') external <- 'Blueprint Encode Data'
  else if (species == 'Mus musculus') external <- 'Mouse RNAseq Data'

  choices <- data.frame(
    value = c('reset', gsub(' ', '', external), anal_options$name),
    itemLabel = c('Reset Labels', external, anal_options$itemLabel),
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
#' @export
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
#' @export
#' @keywords internal
get_exclude_choices <- function(dataset_names, data_dir, anal_colors = NA) {

  if (is.null(dataset_names)) return(NULL)

  # load markers and annotation for each
  annot_paths <- scseq_part_path(data_dir, dataset_names, 'annot')
  marker_paths <- scseq_part_path(data_dir, dataset_names, 'markers')

  annots <- lapply(annot_paths, readRDS)
  clusters <- lapply(annots, function(x) seq(0, length(x)-1))

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
#' @export
#' @keywords internal
get_cluster_choices <- function(clusters, sample_comparison = FALSE, ...) {

  testColor <- get_palette(clusters)

  inds <- seq_along(clusters)
  not.inds <- clusters != inds
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
    # non-formatted for item/hover
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

  } else {
    # show the cell numbers/percentages
    choices$ncells <- cluster_stats$ncells
    choices$pcells <- round(cluster_stats$pcells)
    choices$pspace <- strrep('&nbsp;&nbsp;', max(0, 2 - nchar(choices$pcells)))
  }

  return(choices)
}


#' Get choices data.frame for custom metrics
#'
#' @param scseq \code{SingleCellExperiment}
#'
#' @return data.frame
#' @export
get_metric_choices <- function(scseq) {

  metrics <- scseq@colData
  names <- colnames(metrics)
  names <- names[sapply(metrics, is.logical)]
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

html_space <- function(x) {
  x <- format(as.character(x))
  gsub(' ', '&nbsp;&nbsp;', x)
}

#' Get/Save cluster stats for single-cell related selectizeInputs
#'
#' @param dataset_dir Directory with single cell dataset.
#' @param scseq \code{SingleCellExperiment} object to get/save stats for. if \code{NULL} (Default), will be loaded.
#'
#' @return List with cluster stats
#' @export
get_cluster_stats <- function(dataset_dir = NULL, scseq = NULL, top_tables = NULL, has_replicates = FALSE, use_disk = FALSE, sample_comparison = FALSE) {

  stats_path <- file.path(dataset_dir, 'cluster_stats.rds')
  if (file.exists(stats_path) && use_disk) return(readRDS(stats_path))

  if (is.null(scseq)) scseq <- load_scseq(dataset_dir)

  ncells <- tabulate(scseq$cluster)
  pcells <- ncells / sum(ncells) * 100
  stats <- list(ncells = ncells, pcells = pcells)

  if (sample_comparison) {

    # number of total test and ctrl cells (shown)
    nbins <- length(levels(scseq$cluster))
    is.test <- scseq$orig.ident == 'test'
    stats$ntest <- tabulate(scseq$cluster[is.test], nbins = nbins)
    stats$nctrl <- tabulate(scseq$cluster[!is.test], nbins = nbins)

    # number of test and ctrl cells in each sample (title)
    test <- unique(scseq$batch[is.test])
    ctrl <- unique(scseq$batch[!is.test])

    neach <- tapply(scseq$cluster, list(scseq$batch, scseq$cluster), length)
    neach[is.na(neach)] <- 0

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

  if (use_disk) saveRDS(stats, stats_path)
  return(stats)
}



#' Get contrast choices data.frame for selectize dropdown
#'
#' @param clusters Character vector of cluster names
#' @param test Name of test contrast
#'
#' @return data.frame with columns for rendering selectizeInput contrast choices
#' @export
#' @keywords internal
get_contrast_choices <- function(clusters, test) {

  # group choices are as compared to other clusters
  test_name <- clusters[as.numeric(test)]
  ctrl_names <- clusters[clusters != test_name]
  ctrls <- setdiff(seq_along(clusters), test)

  colours <- get_palette(clusters)
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


#' Add cell percents to gene choices for single cell
#'
#' @param scseq \code{SingleCellExperiment} object
#' @param markers data.frame of marker genes.
#'
#' @return data.frame of all genes, with markers on top and cell percent columns
#' @export
#' @keywords internal
get_gene_choices <- function(markers, qc_metrics = NULL, type = NULL, qc_first = FALSE) {
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

  idx <- match(choices, tx2gene$gene_name)
  choices <- data.table::data.table(label = choices,
                                    value = choices,
                                    description = tx2gene$description[idx],
                                    type = type,
                                    stringsAsFactors = FALSE)




  return(choices)
}




#' Integrate previously saved scseqs
#'
#' Performs integration and saves as a new analysis.
#' Used by \code{explore_scseq_clusters} shiny app.
#'
#' @param sc_dir Directory with saved single-cell datasets.
#' @param test Character vector of test analysis names.
#' @param ctrl Character vector of control analysis names.
#' @param dataset_name Name for new integrated analysis.
#' @param progress optional Shiny \code{Progress} object.
#'
#' @return NULL
#' @export
#' @keywords internal
integrate_saved_scseqs <- function(sc_dir, test, ctrl, exclude_clusters, integration_name, integration_type = c('harmony', 'liger', 'fastMNN'), pairs = NULL, progress = NULL) {
  # for save_scseq_args
  args <- c(as.list(environment()))
  args$progress <- args$sc_dir <- NULL
  args$date <- Sys.time()

  if (is.null(progress)) {
    progress <- list(set = function(value, message = '', detail = '') {
      cat(value, message, detail, '...\n')
    })
  }

  dataset_name <- paste(integration_name, integration_type, sep = '_')

  # save dummy data if testing shiny
  if (isTRUE(getOption('shiny.testmode'))) {
    scseq_data <- list(scseq = NULL, markers = NULL, annot = NULL)
    save_scseq_data(scseq_data, dataset_name, sc_dir, integrated = TRUE)
    return(NULL)
  }

  progress$set(1, detail = 'loading')
  test_scseqs <- load_scseq_subsets(test, sc_dir, exclude_clusters, ident = 'test')
  ctrl_scseqs <- load_scseq_subsets(ctrl, sc_dir, exclude_clusters, ident = 'ctrl')

  # preserve identity of original samples and integrate
  scseqs <- c(test_scseqs, ctrl_scseqs)
  ambient <- get_integrated_ambient(scseqs)

  # make sure all the same species
  species <- unique(sapply(scseqs, function(x) x@metadata$species))
  if(length(species) > 1) stop('Multi-species integration not supported.')

  if (species == 'Homo sapiens') release <- '94'
  else if (species == 'Mus musculus') release <- '98'

  progress$set(2, detail = 'integrating')
  combined <- integrate_scseqs(scseqs, type = integration_type)
  combined$project <- dataset_name

  # retain original QC metrics
  combined <- add_combined_metrics(combined, scseqs)

  rm(scseqs, test_scseqs, ctrl_scseqs); gc()

  # add clusters
  progress$set(3, detail = 'clustering')
  choices <- get_npc_choices(combined, type = 'corrected')

  combined@metadata$species <- species
  combined@metadata$npcs <- choices$npcs
  combined$cluster <- choices$cluster

  # TSNE on corrected reducedDim
  progress$set(4, detail = 'reducing')
  combined <- run_tsne(combined, dimred = 'corrected')

  # add ambient outlier info
  combined <- add_integrated_ambient(combined, ambient)

  progress$set(5, detail = 'getting markers')
  tests <- pairwise_wilcox(combined, block = combined$batch, groups = combined$cluster)
  markers <- get_scseq_markers(tests)

  # top markers for SingleR
  top_markers <- scran::getTopMarkers(tests$statistics, tests$pairs)

  # generate pseudo-bulk so that can exclude counts (large)
  summed <- scater::aggregateAcrossCells(combined,
                                         id = S4Vectors::DataFrame(
                                           cluster = combined$cluster,
                                           batch = combined$batch))


  progress$set(6, detail = 'fitting linear models')
  obj <- combined
  pbulk_esets <- NULL
  has_replicates <- length(unique(combined$batch)) > 2
  if (has_replicates) pbulk_esets <- obj <- construct_pbulk_esets(summed, pairs, species, release)
  lm_fit <- run_limma_scseq(obj)

  progress$set(7, detail = 'saving')
  scseq_data <- list(scseq = combined,
                     species = species,
                     summed = summed,
                     markers = markers,
                     ambient = ambient,
                     tests = tests,
                     pairs = pairs,
                     top_markers = top_markers,
                     has_replicates = has_replicates,
                     founder = integration_name,
                     lm_fit_0svs = lm_fit,
                     pbulk_esets = pbulk_esets,
                     annot = names(markers))

  save_scseq_data(scseq_data, dataset_name, sc_dir, integrated = TRUE)
  save_scseq_args(args, dataset_name, sc_dir)
  rm(scseq_data, summed, markers, ambient, tests, pairs, top_markers, lm_fit, pbulk_esets); gc()


  # don't save raw counts for loom (saved in non-sparse format)
  SummarizedExperiment::assay(combined, 'counts') <- NULL; gc()

  progress$set(8, detail = 'saving loom')
  save_scle(combined, file.path(sc_dir, dataset_name))

  return(NULL)
}

#' Save arguments for integration/subsetting
#'
#' @param args Arguments to save
#' @param dataset_name Name of analysis
#' @param sc_dir Directory with single-cell analyses
#'
#' @return NULL
#' @export
#'
save_scseq_args <- function(args, dataset_name, sc_dir) {
  jsonlite::write_json(args,
                       path = paste0(file.path(sc_dir, dataset_name), '.json'),
                       auto_unbox = TRUE,
                       null = 'null',
                       pretty = TRUE)
}

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
#' @param ident Either \code{'test'} or \code{'ctrl'}
#'
#' @return List of \code{SingleCellExperiment} objects.
#' @export
#' @keywords internal
load_scseq_subsets <- function(dataset_names, sc_dir, exclude_clusters, exclude_metrics = NULL, ident = scseq$project) {

  exclude_datasets <- gsub('^(.+?)_\\d+$', '\\1', exclude_clusters)
  exclude_clusters <- gsub('^.+?_(\\d+)$', '\\1', exclude_clusters)

  scseqs <- list()

  # load each scseq and exclude based on metrics
  for (dataset_name in dataset_names) {
    scseq <- readRDS(scseq_part_path(sc_dir, dataset_name, 'scseq'))

    if (!is.null(exclude_metrics)) {
      cdata <- scseq@colData
      metrics <- readRDS.safe(file.path(sc_dir, dataset_name, 'saved_metrics.rds'))
      if (!is.null(metrics)) cdata <- cbind(cdata, metrics)

      exclude <- cdata[, exclude_metrics, drop = FALSE]
      exclude <- apply(exclude, 1, any)
      scseq <- scseq[, !exclude]
    }
    scseqs[[dataset_name]] <- scseq
  }


  # remove excluded clusters
  for (i in seq_along(scseqs)) {
    dataset_name <- names(scseqs)[i]
    scseq <- scseqs[[dataset_name]]

    # set orig.ident to ctrl/test (integration) or original dataset name (subset)
    scseq$orig.ident <- factor(ident)

    # only remove excluded clusters if present
    is.exclude <- exclude_datasets == dataset_name
    if (any(is.exclude)) {
      exclude <- exclude_clusters[is.exclude]
      scseq <- scseq[, !scseq$cluster %in% exclude]
    }
    scseqs[[dataset_name]] <- scseq
  }

  return(scseqs)
}


#' Save Single Cell RNA-seq data for app
#'
#' @param scseq_data Named list with \code{scseq}, \code{markers}, and/or \code{annot}
#' @param dataset_name The analysis name.
#' @param sc_dir Path to directory with single-cell datasets.
#' @param integrated is the analysis integration. Default is \code{FALSE}
#'
#' @return NULL
#' @export
save_scseq_data <- function(scseq_data, dataset_name, sc_dir, integrated = FALSE) {
  dataset_dir <- file.path(sc_dir, dataset_name)

  if (integrated) {
    # add to integrated if new
    int_path <- file.path(sc_dir, 'integrated.rds')
    int_options <- c(readRDS(int_path), dataset_name)
    saveRDS(unique(int_options), int_path)
  }

  # remove all previous data in case overwriting
  unlink(dataset_dir, recursive = TRUE)

  dir.create(dataset_dir)
  for (type in names(scseq_data)) {

    if (type == 'markers') {
      # save marker data.frames individually for fast loading

      markers <- scseq_data[[type]]
      for (i in names(markers))
        saveRDS(markers[[i]], scseq_part_path(sc_dir, dataset_name, paste0('markers_', i)))

    } else if (type == 'tests') {
      # save pairwise test statistics for fast single group comparisons

      tests <- scseq_data[[type]]
      tests_dir <- file.path(dataset_name, 'tests')
      dir.create(file.path(sc_dir, tests_dir))

      saveRDS(tests$pairs, scseq_part_path(sc_dir, tests_dir, 'pairs'))

      for (i in seq_along(tests$statistics))
        saveRDS(tests$statistics[[i]], scseq_part_path(sc_dir, tests_dir, paste0('statistics_pair', i)))

    } else {
      saveRDS(scseq_data[[type]], scseq_part_path(sc_dir, dataset_name, type))
    }

  }

  return(NULL)
}

#' Validate dataset selection for integration
#'
#' @param test Character vector of test dataset names
#' @param ctrl Character vector of control dataset names
#'
#' @return \code{NULL} if valid, otherwise an error message
#' @export
#' @keywords internal
validate_integration <- function(test, ctrl, pairs) {
  msg <- NULL

  if (!is.null(pairs)) {

    if (!all(pairs$sample %in% c(test, ctrl)))
      msg <- 'Samples missing from pairs csv'
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
#' @export
#' @keywords internal
scseq_part_path <- function(data_dir, dataset_name, part) {
  fname <- paste0(part, '.rds')
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
#' @export
#' @keywords internal
run_drug_queries <- function(top_table, drug_paths, es, ambient = NULL, species = NULL) {

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
    cmap = query_drugs(dprimes, es$cmap),
    l1000_drugs = query_drugs(dprimes, es$l1000_drugs),
    l1000_genes = query_drugs(dprimes, es$l1000_genes)
  )

  saveRDS(res$cmap, drug_paths$cmap)
  saveRDS(res$l1000_drugs, drug_paths$l1000_drugs)
  saveRDS(res$l1000_genes, drug_paths$l1000_genes)

  return(res)
}

#' Load drug effect size matrices for drug queries
#'
#' @return list of matrices
#' @export
#' @keywords internal
load_drug_es <- function() {

  cmap_path  <- system.file('extdata', 'cmap_es_ind.rds', package = 'drugseqr.data', mustWork = TRUE)
  l1000_drugs_path <- system.file('extdata', 'l1000_drugs_es.rds', package = 'drugseqr.data', mustWork = TRUE)
  l1000_genes_path <- system.file('extdata', 'l1000_genes_es.rds', package = 'drugseqr.data', mustWork = TRUE)

  cmap  <- readRDS(cmap_path)
  l1000_drugs <- readRDS(l1000_drugs_path)
  l1000_genes <- readRDS(l1000_genes_path)

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
#' @export
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
#' @export
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
#' get_nearest_row(truth, test)
#' #[1] 4 3
get_nearest_row <- function(truth, test) {
  diffs   <- truth[rep(1:nrow(truth), nrow(test)),] -test[rep(1:nrow(test), each=nrow(truth)),]
  eucdiff <- function(x) sqrt(rowSums(x^2))
  max.col(-matrix(eucdiff(diffs), nrow=nrow(test), byrow=TRUE), "first")
}

#' Get plot to compare original labels for integrated single cell dataset.
#'
#' @param anal Name of original analysis in \code{scseq$project} to plot.
#' @param scseq Integrated \code{SingleCellExperiment} object.
#' @param annot Character vector of original labels for \code{anal}.
#' @param plot \code{ggplot} object from integrated \code{scseq}.
#'
#' @return \code{ggplot} object showing cells in \code{anal} with original labels but coordinates from \code{plot}
#'  of integrated \code{scseq}.
#' @export
#' @keywords internal
#'
get_label_plot <- function(anal, scseq, annot, plot) {

  # current colors used in plot with associated cluster
  cols <- get_palette(levels(scseq$cluster))
  names(cols) <- levels(scseq$cluster)


  # median coordinates for clusters for integrated dataset
  plot_data <- plot$data
  int_coords <- get_tsne_coords(plot_data)

  # use same limits as integrated
  xlims <- range(plot_data$TSNE_1)
  ylims <- range(plot_data$TSNE_2)

  # median coordinates for clusters for individual dataset
  in_anal <- scseq$batch %in% anal
  scseq <- scseq[, in_anal]
  plot_data <- plot_data[colnames(scseq), ]
  plot_data$ident <- scseq$orig.cluster
  anal_coords <- get_tsne_coords(plot_data)

  # closest match for each anal cluster
  match <- get_nearest_row(int_coords[ ,-1], anal_coords[, -1])
  anal_coords <- anal_coords %>%
    mutate(match = int_coords$ident[match]) %>%
    arrange(as.numeric(ident))

  # change identity to previous labels
  cl <- scseq$orig.cluster
  lv <- as.character(sort(unique(as.numeric(cl))))
  names(annot) <- seq_along(annot)

  scseq$cluster <- factor(cl, levels = lv)
  levels(scseq$cluster) <- annot[lv]

  # expand to best match
  anal_coords$cols <- cols[anal_coords$match]

  # replace any duplicates
  is.dup <- which(duplicated(anal_coords$cols))
  if (length(is.dup)) {
    new_cols <- setdiff(get_palette(lv), anal_coords$cols)
    anal_coords[is.dup, 'cols'] <- new_cols[seq_along(is.dup)]
  }

  plot_tsne_cluster(scseq, cols=anal_coords$cols) +
    ggplot2::xlim(xlims) +
    ggplot2::ylim(ylims)
}

#' Get median x-y coordinates for clusters in TSNE plot data
#'
#' @param plot_data data.frame with columns \code{'ident'}, \code{'TSNE_1'}, and \code{'TSNE_2'}.
get_tsne_coords <- function(plot_data) {
  plot_data %>%
    group_by(ident) %>%
    summarise(TSNE_1 = median(TSNE_1),
              TSNE_2 = median(TSNE_2))
}

