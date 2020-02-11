#' Get predicted annotation for label transfer
#'
#' Clusters with an average prediction scores below \code{min.score} retain their original labels.
#'
#' @param ref_preds data.frame generated in \code{\link{labelTransferForm}} on event \code{submit_transfer}
#' @param ref_name Name of reference analysis that labels are transfered from.
#' @param anal_name Name of analysis that labels are transfered to.
#' @param sc_dir Directory containing folders with analyses for \code{ref_name} and \code{anal_name}.

#' @return Character vector of predicted labels from \code{ref_name}.
#' @export
#' @keywords internal
get_pred_annot <- function(ref_preds, ref_name, anal_name, sc_dir) {


  # load query annotation
  query_annot_path <- scseq_part_path(sc_dir, anal_name, 'annot')
  query_annot <- readRDS(query_annot_path)

  senv <- loadNamespace('SingleR')

  if (ref_name %in% ls(senv)) {
    pred_annot <- make.unique(ref_preds, '_')

  } else if (ref_name == '') {
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
diff_abundance <- function(scseq, annot) {

  abundances <- table(scseq$cluster, scseq$batch)
  abundances <- unclass(abundances)
  row.names(abundances) <- annot

  if (ncol(abundances) == 2) return(abundances)

  extra.info <- scseq@colData[match(colnames(abundances), scseq$batch),]
  y.ab <- edgeR::DGEList(abundances, samples=extra.info)

  keep <- edgeR::filterByExpr(y.ab, group=y.ab$samples$orig.ident)
  y.ab <- y.ab[keep,]

  group <- y.ab$samples$orig.ident
  group <- relevel(group, 'ctrl')
  design <- model.matrix(~group)
  colnames(design) <- gsub('^group', '', colnames(design))
  y.ab <- edgeR::estimateDisp(y.ab, design, trend="none")

  fit.ab <- edgeR::glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)

  res <- edgeR::glmQLFTest(fit.ab)
  res <- edgeR::topTags(res, n = Inf)
  as.data.frame(res)
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

  anal_options <- anal_options[anal_options$value != selected_anal, ]

  if (species == 'Homo sapiens') external <- 'Blueprint Encode Data'
  else if (species == 'Mus musculus') external <- 'Mouse RNAseq Data'

  choices <- data.frame(
    value = c(gsub(' ', '', external), anal_options$value),
    label = stringr::str_trunc(c(external, anal_options$value), 35),
    type = c('External Reference', anal_options$type),
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
#' @param anal_names Names of analyses selected for integration
#' @param anal_colors Character vector of colors to indicate analysis
#' @param data_dir Directory with single cell analyses.
#'
#' @return data.frame with columns for rendering selectizeInput include choices
#' @export
#' @keywords internal
get_exclude_choices <- function(anal_names, data_dir, anal_colors = NA) {

  if (is.null(anal_names)) return(NULL)

  # load markers and annotation for each
  annot_paths <- scseq_part_path(data_dir, anal_names, 'annot')
  marker_paths <- scseq_part_path(data_dir, anal_names, 'markers')

  annots <- lapply(annot_paths, readRDS)
  clusters <- lapply(annots, function(x) seq(0, length(x)-1))

  exclude_choices <- lapply(seq_along(anal_names), function(i) {
    data.frame(
      name = stringr::str_trunc(annots[[i]], 27),
      value = paste(anal_names[i], clusters[[i]], sep = '_'),
      anal = anal_names[i],
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
#' @param sample_comparison is this for test vs control comparion? Default is \code{FALSE}.
#'
#' @return data.frame with columns for rendering selectizeInput cluster choices
#' @export
#' @keywords internal
get_cluster_choices <- function(clusters, dataset_dir, scseq = NULL, sample_comparison = FALSE, top_tables = NULL, has_replicates = FALSE) {

  testColor <- get_palette(clusters)

  # cluster choices are the clusters themselves
  # value is original cluster number so that saved pathway analysis name
  # isn't affected by updates to cluster annotation
  choices <- data.frame(name = stringr::str_trunc(clusters, 27),
                        value = seq(1, along.with = clusters),
                        label = clusters,
                        testColor,
                        row.names = NULL, stringsAsFactors = FALSE)

  cluster_stats <- get_cluster_stats(dataset_dir, scseq, top_tables = top_tables, has_replicates = has_replicates)

  if (sample_comparison) {
    # non-formatted for item/hover
    choices$ntest <- cluster_stats$ntest
    choices$nctrl <- cluster_stats$nctrl
    choices$nsig  <- cluster_stats$nsig
    choices$nbig  <- cluster_stats$nbig
    choices$ntest_each <- cluster_stats$ntest_each
    choices$nctrl_each <- cluster_stats$nctrl_each

    # formatted for options
    choices$nctrlf <- format(cluster_stats$nctrl)
    choices$nctrlf <- gsub(' ', '&nbsp;&nbsp;', choices$nctrlf)
    choices$nsigf  <- format(cluster_stats$nsig)
    choices$nsigf  <- gsub(' ', '&nbsp;&nbsp;', choices$nsigf)
    choices$nbigf  <- format(cluster_stats$nbig)
    choices$nbigf  <- gsub(' ', '&nbsp;&nbsp;', choices$nbigf)

  } else {
    # show the cell numbers/percentages
    choices$ncells <- cluster_stats$ncells
    choices$pcells <- round(cluster_stats$pcells)
    choices$pspace <- strrep('&nbsp;&nbsp;', 2 - nchar(choices$pcells))
  }

  return(choices)
}

#' Get/Save cluster stats for single-cell related selectizeInputs
#'
#' @param dataset_dir Directory with single cell dataset.
#' @param scseq \code{SingleCellExperiment} object to get/save stats for. if \code{NULL} (Default), will be loaded.
#'
#' @return List with cluster stats
#' @export
get_cluster_stats <- function(dataset_dir, scseq = NULL, top_tables = NULL, has_replicates = FALSE) {

  # return previously saved stats if exists
  stats_path <- file.path(dataset_dir, 'cluster_stats.rds')
  if (file.exists(stats_path)) return(readRDS(stats_path))

  # otherwise generate and save
  if (is.null(scseq)) scseq <- load_scseq(dataset_dir)

  ncells <- tabulate(scseq$cluster)
  pcells <- ncells / sum(ncells) * 100
  stats <- list(ncells = ncells, pcells = pcells)

  is.integrated <- !is.null(top_tables)
  if (is.integrated) {

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
  if (is.integrated & has_replicates) {
    nsig <- rep(0, nbins)
    names(nsig) <- seq_len(nbins)

    test_clusters <- names(top_tables)
    nsig[test_clusters] <- sapply(top_tables, function(tt) {sum(tt$adj.P.Val.Amb < 0.05 & !tt$ambient)})
    stats$nsig <- nsig
  }

  # show number of non-ambient with logFC > 1
  if (is.integrated) {
    nbig <- rep(0, nbins)
    names(nbig) <- seq_len(nbins)

    test_clusters <- names(top_tables)
    nbig[test_clusters] <- sapply(top_tables, function(tt) {sum(abs(tt$logFC) > 1 & !tt$ambient)})
    stats$nbig <- nbig
  }

  saveRDS(stats, stats_path)
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

  contrast_choices <- data.frame(test = stringr::str_trunc(test_name, 11, ellipsis = '..'),
                                 ctrl = stringr::str_trunc(c('all', ctrl_names), 11, ellipsis = '..'),
                                 name = test_name,
                                 value = c(test, paste0(test, '-vs-', ctrls)),
                                 testColor = colours[test_name],
                                 ctrlColor = c('white', colours[ctrl_names]), row.names = NULL, stringsAsFactors = FALSE)

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
get_gene_choices <- function(markers, type = 'samples') {

  markers$label <- markers$value <- row.names(markers)

  # add description for title
  idx <- match(row.names(markers), tx2gene$gene_name)
  markers$description <- tx2gene$description[idx]

  markers <- markers[, c('label', 'value', 'description')]
  markers <- data.table::data.table(markers)
  return(markers)
}




#' Integrate previously saved scseqs
#'
#' Performs integration and saves as a new analysis.
#' Used by \code{explore_scseq_clusters} shiny app.
#'
#' @param sc_dir Directory with saved single-cell datasets.
#' @param test Character vector of test analysis names.
#' @param ctrl Character vector of control analysis names.
#' @param anal_name Name for new integrated analysis.
#' @param progress optional Shiny \code{Progress} object.
#'
#' @return NULL
#' @export
#' @keywords internal
integrate_saved_scseqs <- function(sc_dir, test, ctrl, exclude_clusters, anal_name, pairs = NULL, updateProgress = NULL) {

  # save dummy data if testing shiny
  if (isTRUE(getOption('shiny.testmode'))) {
    scseq_data <- list(scseq = NULL, markers = NULL, annot = NULL)
    save_scseq_data(scseq_data, anal_name, sc_dir, integrated = TRUE)
    return(NULL)
  }

  # default updateProgress and number of steps
  if (is.null(updateProgress)) updateProgress <- function(...) {NULL}
  n = 8

  updateProgress(1/n, 'loading')
  test_scseqs <- load_scseqs_for_integration(test, exclude_clusters = exclude_clusters, sc_dir = sc_dir, ident = 'test')
  ctrl_scseqs <- load_scseqs_for_integration(ctrl, exclude_clusters = exclude_clusters, sc_dir = sc_dir, ident = 'ctrl')

  # preserve identity of original samples and integrate
  scseqs <- c(test_scseqs, ctrl_scseqs)
  ambient <- get_integrated_ambient(scseqs)

  # make sure all the same species
  species <- unique(sapply(scseqs, function(x) x@metadata$species))
  if(length(species) > 1) stop('Multi-species integration not supported.')

  if (species == 'Homo sapiens') release <- '94'
  else if (species == 'Mus musculus') release <- '98'

  updateProgress(2/n, 'integrating')
  combined <- integrate_scseqs(scseqs)

  # retain original doublet scores (needs scaling?)
  combined$doublet_score <- unlist(sapply(scseqs, `[[`, 'doublet_score'))
  rm(scseqs, test_scseqs, ctrl_scseqs); gc()

  # add clusters
  updateProgress(3/n, 'clustering')
  choices <- get_npc_choices(combined, type = 'corrected')

  combined@metadata$species <- species
  combined@metadata$npcs <- choices$npcs
  combined$cluster <- choices$cluster

  # TSNE on corrected reducedDim
  updateProgress(4/n, 'reducing')
  combined <- run_tsne(combined, dimred = 'corrected')

  # add ambient outlier info
  combined <- add_integrated_ambient(combined, ambient)

  updateProgress(5/n, 'getting markers')
  tests <- pairwise_wilcox(combined, block = combined$batch, groups = combined$cluster)
  markers <- get_scseq_markers(tests)

  # top markers for SingleR
  top_markers <- scran::getTopMarkers(tests$statistics, tests$pairs)

  # generate pseudo-bulk so that can exclude counts (large)
  summed <- scater::aggregateAcrossCells(combined,
                                         id = S4Vectors::DataFrame(
                                           cluster = combined$cluster,
                                           batch = combined$batch))

  SummarizedExperiment::assay(combined, 'counts') <- NULL; gc()

  updateProgress(6/n, 'fitting linear models')
  obj <- combined
  pbulk_esets <- NULL
  has_replicates <- length(unique(combined$batch)) > 2
  if (has_replicates) pbulk_esets <- obj <- construct_pbulk_esets(summed, pairs, species, release)
  lm_fit <- run_limma_scseq(obj)


  updateProgress(7/n, 'saving')
  scseq_data <- list(scseq = combined,
                     species = species,
                     summed = summed,
                     markers = markers,
                     ambient = ambient,
                     tests = tests,
                     top_markers = top_markers,
                     has_replicates = has_replicates,
                     lm_fit_0svs = lm_fit,
                     pbulk_esets = pbulk_esets,
                     annot = names(markers))

  save_scseq_data(scseq_data, anal_name, sc_dir, integrated = TRUE)

  updateProgress(8/n, 'saving loom')
  save_scle(scseq, file.path(sc_dir, anal_name))

  return(NULL)
}


#' Load scRNA-Seq datasets for integration
#'
#' Sets orig.ident to \code{ident}.
#'
#' @param anal_names Character vector of single cell analysis names to load.
#' @param sc_dir The directory with single-cell datasets
#' @param ident Either \code{'test'} or \code{'ctrl'}
#'
#' @return List of \code{SingleCellExperiment} objects.
#' @export
#' @keywords internal
load_scseqs_for_integration <- function(anal_names, exclude_clusters, sc_dir, ident) {

  exclude_anals <- gsub('^(.+?)_\\d+$', '\\1', exclude_clusters)
  exclude_clusters <- gsub('^.+?_(\\d+)$', '\\1', exclude_clusters)

  # TODO: make integration work with SingleCellLoomExperiment
  scseqs <- list()
  for (anal in anal_names)
    scseqs[[anal]] <- readRDS(scseq_part_path(sc_dir, anal, 'scseq'))


  for (i in seq_along(scseqs)) {
    anal <- names(scseqs)[i]
    scseq <- scseqs[[anal]]

    # set orig.ident to ctrl/test
    scseq$orig.ident <- factor(ident)

    # only remove excluded clusters if present
    is.exclude <- exclude_anals == anal
    if (any(is.exclude)) {
      exclude <- exclude_clusters[is.exclude]
      scseq <- scseq[, !scseq$cluster %in% exclude]
    }
    scseqs[[anal]] <- scseq
  }

  return(scseqs)
}


#' Save Single Cell RNA-seq data for app
#'
#' @param scseq_data Named list with \code{scseq}, \code{markers}, and/or \code{annot}
#' @param anal_name The analysis name.
#' @param sc_dir Path to directory with single-cell datasets.
#' @param integrated is the analysis integration. Default is \code{FALSE}
#'
#' @return NULL
#' @export
save_scseq_data <- function(scseq_data, anal_name, sc_dir, integrated = FALSE) {
  anal_dir <- file.path(sc_dir, anal_name)

  if (integrated) {
    # add to integrated if new
    int_path <- file.path(sc_dir, 'integrated.rds')
    int_options <- c(readRDS(int_path), anal_name)
    saveRDS(unique(int_options), int_path)

    # remove all previous data in case overwriting
    unlink(anal_dir, recursive = TRUE)
  }

  dir.create(anal_dir)
  for (type in names(scseq_data)) {

    if (type == 'markers') {
      # save marker data.frames individually for fast loading

      markers <- scseq_data[[type]]
      for (i in names(markers))
        saveRDS(markers[[i]], scseq_part_path(sc_dir, anal_name, paste0('markers_', i)))

    } else if (type == 'tests') {
      # save pairwise test statistics for fast single group comparisons

      tests <- scseq_data[[type]]
      tests_dir <- file.path(anal_name, 'tests')
      dir.create(file.path(sc_dir, tests_dir))

      saveRDS(tests$pairs, scseq_part_path(sc_dir, tests_dir, 'pairs'))

      for (i in seq_along(tests$statistics))
        saveRDS(tests$statistics[[i]], scseq_part_path(sc_dir, tests_dir, paste0('statistics_pair', i)))

    } else {
      saveRDS(scseq_data[[type]], scseq_part_path(sc_dir, anal_name, type))
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
validate_integration <- function(test, ctrl, anal_name, anal_options, pairs) {
  msg <- NULL

  if (is.null(anal_name) || anal_name == '') {
    msg <- 'Provide a name for integrated analysis'

  } else if (is.null(test) || is.null(ctrl)) {
    msg <- 'Need control and test datasets'

  } else if (!is.null(pairs)) {

    if (!all(pairs$sample %in% c(test, ctrl)))
      msg <- 'Samples missing from pairs csv'
  }

  return(msg)
}


#' Get path to saved scseq part
#'
#' @param data_dir Path to directory with analyses.
#' @param anal_name Name of analysis.
#' @param part either \code{'annot'}, \code{'scseq'}, or \code{'markers'}.
#'
#' @return Path to analysis \code{part}.
#' @export
#' @keywords internal
scseq_part_path <- function(data_dir, anal_name, part) {
  fname <- paste0(part, '.rds')
  file.path(data_dir, anal_name, fname)
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

