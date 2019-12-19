#' Get predicted annotation for label transfer
#'
#' Clusters with an average prediction scores below \code{min.score} retain their original labels.
#'
#' @param ref_preds data.frame generated in \code{\link{labelTransferForm}} on event \code{submit_transfer}
#' @param ref_name Name of reference analysis that labels are transfered from.
#' @param anal_name Name of analysis that labels are transfered to.
#' @param sc_dir Directory containing folders with analyses for \code{ref_name} and \code{anal_name}.
#' @param min.score Minimum average prediction score used for label transfer.
#'
#' @return Character vector of predicted labels from \code{ref_name}.
#' @export
#' @keywords internal
get_pred_annot <- function(ref_preds, ref_name, anal_name, sc_dir, min.score = 0) {

  # load reference and query annotation
  query_annot_path <- scseq_part_path(sc_dir, anal_name, 'annot')
  query_annot <- readRDS(query_annot_path)

  # for resetting annotation
  if (ref_name == '') {
    pred_annot <- as.character(seq_along(query_annot)-1)

  } else {
    ref_annot_path <- scseq_part_path(sc_dir, ref_name, 'annot')
    ref_annot <- readRDS(ref_annot_path)

    # load query annotation

    # get predicted annotation
    pred.idx <- as.numeric(ref_preds$predicted.id) + 1
    pred_annot <- ref_annot[pred.idx]

    # use original query annotation below min score threshold
    # not sure of utility so default is take all preds
    poor.pred <- ref_preds$mean.score < min.score
    pred_annot[poor.pred] <- query_annot[poor.pred]
    pred_annot <- make.unique(pred_annot, '_')
  }
  return(pred_annot)
}


#' Get percentage of cells expressing each gene
#'
#' @param scseq \code{Seurat} object
#' @param ident.1 Test group level in \code{Idents(scseq)}.
#' @param ident.2 Control group level in \code{Idents(scseq)}.
#'
#' @return data.frame with columns pct.1 and pct.2 indicating fraction of cells that express each gene.
#' @export
#' @keywords internal
get_cell_pcts <- function(scseq, ident.1, ident.2) {
  assay <- get_scseq_assay(scseq)
  data <- scseq[[assay]]@data

  cells <- Seurat::Idents(scseq)
  cells.1 <- cells %in% ident.1
  cells.2 <- cells %in% ident.2

  pct.1 <- round(
    x = Matrix::rowSums(x = data[, cells.1, drop = FALSE] > 0) / sum(cells.1),
    digits = 3
  )

  pct.2 <- round(
    x = Matrix::rowSums(x = data[, cells.2, drop = FALSE] > 0) / sum(cells.2),
    digits = 3
  )

  return(cbind(pct.1, pct.2))
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
get_label_transfer_choices <- function(anal_options, selected_anal, preds) {

  anal_options <- lapply(anal_options, setdiff, selected_anal)

  choices <- data.frame(
    label = c(NA, unlist(anal_options, use.names = FALSE)),
    type = c(NA,
             rep('Individual', length(anal_options$Individual)),
             rep('Integrated', length(anal_options$Integrated)))
  )


  choices$value <- choices$label
  choices$preds <- choices$label %in% names(preds)

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
#' @param scseq \code{Seurat} object
#' @param value Character vector for value column which is returned from \code{selectizeInput}. Default is \code{clusters}.
#' @param sample_comparison is this for test vs control comparion? Default is \code{FALSE}.
#'
#' @return data.frame with columns for rendering selectizeInput cluster choices
#' @export
#' @keywords internal
get_cluster_choices <- function(clusters, anal_name, sc_dir, sample_comparison = FALSE) {

  testColor <- get_palette(clusters)

  # cluster choices are the clusters themselves
  # value is original cluster number so that saved pathway analysis name
  # isn't affected by updates to cluster annotation
  choices <- data.frame(name = stringr::str_trunc(clusters, 27),
                        value = seq(0, along.with = clusters),
                        label = clusters,
                        testColor,
                        row.names = NULL, stringsAsFactors = FALSE)

  cluster_stats <- get_cluster_stats(sc_dir, anal_name)

  if (sample_comparison) {
    choices$ntest <- cluster_stats$ntest
    choices$nctrl <- format(cluster_stats$nctrl)
    choices$nctrl <- gsub(' ', '&nbsp;&nbsp;', choices$nctrl)

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
#' @param sc_dir Path to directory containing single-cell dataset folders
#' @param anal_name Name of single cell analysis to get/save stats for
#' @param scseq \code{Seurat} object to get/save stats for. if \code{NULL} (Default), will be loaded.
#'
#' @return List with cluster stats
#' @export
get_cluster_stats <- function(sc_dir, anal_name, scseq = NULL) {

  # return previously saved stats if exists
  stats_path <- scseq_part_path(sc_dir, anal_name, 'cluster_stats')
  if (file.exists(stats_path)) return(readRDS(stats_path))

  # otherwise generate and save
  if (is.null(scseq)) {
    scseq_path <- scseq_part_path(sc_dir, anal_name, 'scseq')
    scseq <- readRDS(scseq_path)
  }

  ncells <- tabulate(scseq$seurat_clusters)
  pcells <- ncells / sum(ncells) * 100
  stats <- list(ncells = ncells, pcells = pcells)

  is.integrated <- 'integrated' %in% names(scseq@assays)
  if (is.integrated) {
    stats$ntest <- tabulate(scseq$seurat_clusters[scseq$orig.ident == 'test'])
    stats$nctrl <- tabulate(scseq$seurat_clusters[scseq$orig.ident == 'ctrl'])
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
  test_name <- clusters[as.numeric(test)+1]
  ctrl_names <- clusters[clusters != test_name]
  ctrls <- setdiff(seq(0, along.with = clusters), test)

  colours <- get_palette(clusters)
  names(colours) <- clusters

  contrast_choices <- data.frame(test = stringr::str_trunc(test_name, 11, ellipsis = '..'),
                                 ctrl = stringr::str_trunc(c('all', ctrl_names), 11, ellipsis = '..'),
                                 value = c(test, paste0(test, ' vs ', ctrls)),
                                 testColor = colours[test_name],
                                 ctrlColor = c('white', colours[ctrl_names]), row.names = NULL, stringsAsFactors = FALSE)

  return(contrast_choices)

}


#' Add cell percents to gene choices for single cell
#'
#' @param scseq Seurat object
#' @param markers data.frame of marker genes.
#' @param selected_cluster Character vector indicating selected cluster.
#' @param comparison_type Either \code{'samples'} or \code{'clusters'}
#'
#' @return data.frame of all genes, with markers on top and cell percent columns
#' @export
#' @keywords internal
get_gene_choices <- function(scseq, markers, selected_cluster, comparison_type) {
  clusters <- Seurat::Idents(scseq)

  # get cell.pcts for all genes (previously just had for markers)
  # allows selecting non-marker genes (at bottom of list)
  if (comparison_type == 'clusters') {

    markers <- markers[, c('pct.1', 'pct.2')]
    idents <- strsplit(selected_cluster, ' vs ')[[1]]

    if (length(idents) == 2) {
      cell_pcts <- get_cell_pcts(scseq, idents[1], idents[2])
    } else {
      cell_pcts <- get_cell_pcts(scseq, idents, setdiff(levels(clusters), idents))
    }

    cell_pcts <- cell_pcts[!row.names(cell_pcts) %in% row.names(markers), ]
    markers <- rbind(markers, cell_pcts)

  } else {
    Seurat::Idents(scseq) <- scseq$orig.ident
    scseq <- scseq[, clusters %in% selected_cluster]
    cell_pcts <- get_cell_pcts(scseq, 'test', 'ctrl')
    # most positive on top
    markers <- markers[order(markers$t, decreasing = TRUE), ]
    markers <- as.data.frame(cell_pcts)[row.names(markers), ]
  }

  markers$pct.1 <- round(markers$pct.1 * 100)
  markers$pct.2 <- round(markers$pct.2 * 100)

  pspace <- 2 - nchar(markers$pct.2)
  pspace[pspace < 0] <- 0
  markers$pspace <- strrep('&nbsp;&nbsp;', pspace)

  markers$label <- markers$value <- row.names(markers)

  # add description for title
  markers$description <- tx2gene$description[match(row.names(markers), tx2gene$gene_name)]
  return(markers)
}




#' Integrate previously saved scseqs
#'
#' Performs integration and saves as a new analysis.
#' Used by \code{explore_scseq_clusters} shiny app.
#'
#' @param data_dir Directory with saved analyses.
#' @param test Character vector of test analysis names.
#' @param ctrl Character vector of control analysis names.
#' @param anal_name Name for new integrated analysis.
#' @param progress optional Shiny \code{Progress} object.
#'
#' @return NULL
#' @export
#' @keywords internal
integrate_saved_scseqs <- function(data_dir, test, ctrl, exclude_clusters, anal_name, updateProgress = NULL, use_scalign = FALSE) {

  reduction <- ifelse(use_scalign, 'embed', 'pca')
  dims <- if (use_scalign) 1:32 else 1:30

  # save dummy data if testing shiny
  if (isTRUE(getOption('shiny.testmode'))) {
    scseq_data <- list(scseq = NULL, markers = NULL, annot = NULL)
    save_scseq_data(scseq_data, anal_name, data_dir, integrated = TRUE)
    return(NULL)
  }

  # default updateProgress and number of steps
  if (is.null(updateProgress)) updateProgress <- function(...) {NULL}
  n = 6

  updateProgress(1/n, 'loading')
  test_scseqs <- load_scseqs_for_integration(test, exclude_clusters = exclude_clusters, data_dir = data_dir, ident = 'test')
  ctrl_scseqs <- load_scseqs_for_integration(ctrl, exclude_clusters = exclude_clusters, data_dir = data_dir, ident = 'ctrl')

  # preserve identity of original samples and integrate
  scseqs <- c(test_scseqs, ctrl_scseqs)
  scseqs <- add_project_scseqs(scseqs)

  updateProgress(2/n, 'integrating')
  combined <- integrate_scseqs(scseqs, use_scalign = use_scalign)
  rm(scseqs); gc()

  updateProgress(3/n, 'clustering')
  combined <- add_scseq_clusters(combined, reduction = reduction, dims = dims)

  updateProgress(4/n, 'reducing')
  combined <- run_umap(combined)

  updateProgress(5/n, 'getting markers')
  markers <- get_scseq_markers(combined)

  updateProgress(6/n, 'saving')
  scseq_data <- list(scseq = combined, markers = markers, annot = names(markers))
  save_scseq_data(scseq_data, anal_name, data_dir, integrated = TRUE)

  # get and save cluster stats
  get_cluster_stats(data_dir, anal_name, combined)

  return(NULL)
}


#' Load scRNA-Seq datasets for integration
#'
#' Will restore SCT assay if was saved seperately. Downsamples very large datasets.
#' Also sets orig.ident to \code{ident}.
#'
#' @param anal_names Character vector of single cell analysis names to load.
#' @param data_dir The directory with single cell RNA seq datasets
#' @param ident Either \code{'test'} or \code{'ctrl'}
#'
#' @return List of \code{Seurat} objects.
#' @export
#' @keywords internal
load_scseqs_for_integration <- function(anal_names, exclude_clusters, data_dir, ident) {

  exclude_anals <- gsub('^(.+?)_\\d+$', '\\1', exclude_clusters)
  exclude_clusters <- gsub('^.+?_(\\d+)$', '\\1', exclude_clusters)

  scseqs <- list()
  for (anal in anal_names) {
    scseqs[[anal]] <- load_saved_scseq(anal, data_dir)
  }

  for (i in seq_along(scseqs)) {
    anal <- names(scseqs)[i]
    scseq <- scseqs[[anal]]

    # set orig.ident to ctrl/test
    scseq$orig.ident <- factor(ident)

    # only remove excluded clusters if present
    is.exclude <- exclude_anals == anal
    if (any(is.exclude)) {
      exclude <- exclude_clusters[is.exclude]
      scseq <- scseq[, !scseq$seurat_clusters %in% exclude]
    }
    scseqs[[anal]] <- scseq
  }



  return(scseqs)
}

#' Load saved single cell RNA-seq dataset
#'
#' Used for both dataset integration and label transfer.
#'
#' @param anal Name of single cell analysis and containing folder.
#' @param data_dir Path to directory containing \code{anal} folder.
#' @param downsample Should the loaded Seurat object be downsampled? For reducing speed/memory burden. Default is FALSE.
#'
#' @return \code{Seurat} object
#' @export
#' @keywords internal
load_saved_scseq <- function(anal, data_dir, downsample = FALSE) {
  # load scseq
  scseq_path <- scseq_part_path(data_dir, anal, 'scseq')
  scseq <- readRDS(scseq_path)

  # add original clusters
  scseq$orig_clusters <- scseq$seurat_clusters

  # restore SCT assay
  sct_path <- scseq_part_path(data_dir, anal, 'sct')
  if (file.exists(sct_path)) {
    sct <- readRDS(sct_path)
    scseq[['SCT']] <- sct
  }

  # restore counts slot
  scseq[['RNA']]@counts <- scseq[['RNA']]@data

  # downsample very large datasets
  if (downsample) scseq <- downsample_scseq(scseq)

  # bug in Seurat, need for integration
  scseq[['SCT']]@misc$vst.out$cell_attr <- scseq[['SCT']]@misc$vst.out$cell_attr[colnames(scseq), ]

  return(scseq)
}


#' Downsample very large scseq objects for integration
#'
#' Used by \code{load_scseqs_for_integration}.
#'
#' @param scseq \code{Seurat} object
#' @param max.cells Maximum number of cells to keep. Default is 10000.
#' @param seed Integer used for reproducibility.
#'
#' @return \code{scseq} with maximum \code{max.cells} cells.
#' @export
#' @keywords internal
downsample_scseq <- function(scseq, max.cells = 1000, seed = 0L) {
  if (ncol(scseq) > max.cells) {
    set.seed(seed)
    scseq <- subset(scseq, cells = sample(Seurat::Cells(scseq), max.cells))
    gc()
  }

  return(scseq)
}



#' Adds project name to meta.data of scseq
#'
#' Used to keep track of which sample is which for integrated datasets. This is used for generating pseudo-bulk
#' counts for each sample.
#'
#' @param scseqs List of \code{Seurat} objects.
#'
#' @return \code{scseqs} with \code{meta.data$project} column equal to \code{project.name} slot.
#' @export
#' @keywords internal
add_project_scseqs <- function(scseqs) {

  # preserve identities of original samples
  projects <- sapply(scseqs, function(x) x@project.name)
  projects <- make.unique(projects)

  for (i in seq_along(projects))
    scseqs[[i]]@meta.data$project <- projects[i]

  return(scseqs)

}


#' Save Single Cell RNA-seq data for app
#'
#' @param scseq_data Named list with \code{scseq}, \code{markers}, and/or \code{annot}
#' @param anal_name The analysis name.
#' @param data_dir Path to directory to save in
#' @param integrated is the analysis integration. Default is \code{FALSE}
#'
#' @return NULL
#' @export
save_scseq_data <- function(scseq_data, anal_name, data_dir, integrated = FALSE, reduce_size = FALSE) {
  if (integrated) {
    int_path <- file.path(data_dir, 'integrated.rds')
    int_options <- readRDS(int_path)
    saveRDS(c(int_options, anal_name), int_path)
  }


  if (reduce_size) {
    # optionally seperate larger parts off
    scseq <- scseq_data$scseq
    sct <- scseq[['SCT']]

    # keep @data which has corrected log counts for visualization
    scseq[['SCT']]@scale.data <- matrix(nrow = 0, ncol = 0)
    scseq[['SCT']]@counts <- matrix(nrow = 0, ncol = 0)
    scseq[['SCT']]@misc <- NULL

    # redundant
    stopifnot(all.equal(scseq[['RNA']]@counts, scseq[['RNA']]@data))
    scseq[['RNA']]@counts <- matrix(nrow = 0, ncol = 0)

    scseq_data$scseq <- scseq
    scseq_data$sct <- sct
  }


  dir.create(file.path(data_dir, anal_name))
  for (type in names(scseq_data)) {
    saveRDS(scseq_data[[type]], scseq_part_path(data_dir, anal_name, type))
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
validate_integration <- function(test, ctrl, anal_name, anal_options) {
  msg <- NULL

  if (is.null(anal_name) || anal_name == '') {
    msg <- 'Provide a name for integrated analysis'

  } else if (anal_name %in% unlist(anal_options)) {
    msg <- 'Analysis name already exists'

  } else if (is.null(test) || is.null(ctrl)) {
    msg <- 'Need control and test datasets'
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


#' Get genes that are ambient in at least one test and control sample
#'
#' @param scseqs List of Seurat objects.
#'
#' @return List with test and control ambient genes
#' @export
#' @keywords interal
get_integrated_ambient <- function(scseqs) {

  # datasets that are test samples
  is.test <- sapply(scseqs, function(x) levels(x$orig.ident) == 'test')

  # genes that are ambient in at least one test sample
  ambient.test <- lapply(scseqs[is.test], function(x) x[['RNA']]@meta.features)
  ambient.test <- lapply(ambient.test, function(x) row.names(x)[x$out_ambient])
  ambient.test <- unique(unlist(ambient.test))

  # genes that are ambient in at least one ctrl sample
  ambient.ctrl <- lapply(scseqs[!is.test], function(x) x[['RNA']]@meta.features)
  ambient.ctrl <- lapply(ambient.ctrl, function(x) row.names(x)[x$out_ambient])
  ambient.ctrl <- unique(unlist(ambient.ctrl))

  return(list(test = ambient.test, ctrl = ambient.ctrl))
}


#' Integrate Seurat objects
#'
#' @param scseqs List of Seurat objects
#' @param use_scalign whether or not to use scAlign. Default is \code{FALSE}
#'
#' @return Integrated Seurat object.
#' @export
#' @keywords internal
integrate_scseqs <- function(scseqs, use_scalign = FALSE) {

  genes  <- Seurat::SelectIntegrationFeatures(object.list = scseqs, nfeatures = 3000)
  scseqs <- Seurat::PrepSCTIntegration(object.list = scseqs, anchor.features = genes)

  ambient <- get_integrated_ambient(scseqs)

  if (use_scalign) {
    combined <- scalign_integrate(scseqs, genes)

  } else {
    combined <- cca_integrate(scseqs, genes)
  }


  # add ambient outlier info
  combined <- add_integrated_ambient(combined, ambient)

  return(combined)
}

#' Integrate with CCA from Seurat
#'
#' Uses SCT method
#'
#' @param scseqs Seurat objects
#' @param genes Character vector of genes to integrate on
#'
#' @return Seurat object
#' @export
#' @keywords internal
cca_integrate <- function(scseqs, genes) {

  k.filter <- min(200, min(sapply(scseqs, ncol)))
  k.score <- min(c(sapply(scseqs, ncol)-1, 30))
  anchors <- Seurat::FindIntegrationAnchors(scseqs, k.filter = k.filter, normalization.method = "SCT",
                                            anchor.features = genes, dims = 1:k.score, k.score = k.score)

  combined <- Seurat::IntegrateData(anchors, normalization.method = "SCT")
  combined$orig.ident <- factor(combined$orig.ident)

  return(combined)

}

#' Integrate Single Cell datasets using scAlign
#'
#' @param scseqs List of \code{Seurat} objects
#' @param genes Highly variable genes to integrate with
#' @importFrom SingleCellExperiment colData
#'
#' @export
#' @keywords internal
#'
scalign_integrate <- function(scseqs, genes) {

  common_meta <- lapply(scseqs, function(x) colnames(x@meta.data))
  common_meta <- Reduce(intersect, common_meta)

  scalign.list <- lapply(scseqs, function(x) {
    SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = x[['SCT']]@data[genes, ],
                                                             scale.data = x[['SCT']]@scale.data[genes, ]),
                                               colData = x@meta.data[, common_meta])
  })

  scalign <- scAlign::scAlignCreateObject(sce.objects = scalign.list, project.name = "sjia")

  scalign <- scAlign::scAlignMulti(scalign,
                                   encoder.data="scale.data",
                                   decoder.data="logcounts",
                                   supervised='none',
                                   run.encoder=TRUE,
                                   run.decoder=TRUE,
                                   log.results=TRUE,
                                   log.dir=file.path('./tmp'),
                                   device="GPU")

  scseq <- scalign_to_scseq(scalign)

  return(scseq)


}



#' Convert scAlign integrated dataset into Seurat object
#'
#' @param scalign SingleCellExperiment returned from \code{scalign_integrate}
#'
#' @return Seurat object
#' @export
#' @keywords internal
scalign_to_scseq <- function(scalign) {

  # empty RNA assay
  empty_counts <- matrix(nrow = 0, ncol = ncol(scalign),dimnames = list(NULL, colnames(scalign)))

  scseq <- Seurat::CreateSeuratObject(
    empty_counts,
    meta.data = as.data.frame(scalign@colData)
  )

  # add SCT assay
  scseq[['SCT']] <- Seurat::CreateAssayObject(data = assay(scalign, 'logcounts'))
  scseq[['SCT']]@scale.data <- assay(scalign, 'scale.data')

  # embedding as reduced dim
  embed <- SingleCellExperiment::reducedDim(scalign, 'ALIGNED-GENE')
  row.names(embed) <- colnames(scalign)
  colnames(embed) <- paste0('EMBED_', seq_len(ncol(embed)))
  scseq@reductions$embed <- Seurat::CreateDimReducObject(embed, assay = 'SCT', key = 'EMBED_')

  # each decoder output logcount matrix as seperate assay
  decodes <- setdiff(SingleCellExperiment::reducedDimNames(scalign), 'ALIGNED-GENE')
  for (decode in decodes) {
    logcounts <- t(SingleCellExperiment::reducedDim(scalign, decode))
    dimnames(logcounts) <- dimnames(scalign)
    scseq[[decode]] <- Seurat::CreateAssayObject(data = logcounts)
  }

  return(scseq)
}


#' Get assay to use from Seurat object.
#'
#' Default is 'SCT' but if not present then 'RNA'
#'
#' @param scseq Seurat object
#'
#' @return either \code{'SCT'} or \code{'RNA'}
#' @export
#' @keywords internal
get_scseq_assay <- function(scseq) {
  ifelse('SCT' %in% names(scseq@assays), 'SCT', 'RNA')
}


#' Run drug queries
#'
#' Used by run_comparison
#'
#' @param res_paths List of paths with names \code{'cmap'} \code{'l1000'} and \code{'anal'} to
#' saved cmap, l1000, and differential expression analysis results.
#' @param session Shiny session object used for progress bar.
#' @param res List to store drug query results in.
#' @param ambient Character vector of ambient genes to exclude from drug queries.
#'
#' @return \code{res} with drug query results added to \code{'cmap'} \code{'l1000'} slots.
#' @export
#' @keywords internal
run_drug_queries <- function(top_table, res_paths, session, res = list(), ambient = NULL) {

  progress <- Progress$new(session, min = 0, max = 4)
  progress$set(message = "Querying drugs", value = 1)
  on.exit(progress$close())

  cmap_path  <- system.file('extdata', 'cmap_es_ind.rds', package = 'drugseqr.data', mustWork = TRUE)
  l1000_drugs_path <- system.file('extdata', 'l1000_drugs_es.rds', package = 'drugseqr.data', mustWork = TRUE)
  l1000_genes_path <- system.file('extdata', 'l1000_genes_es.rds', package = 'drugseqr.data', mustWork = TRUE)

  cmap_es  <- readRDS(cmap_path)
  progress$inc(1)
  l1000_drugs_es <- readRDS(l1000_drugs_path)
  l1000_genes_es <- readRDS(l1000_genes_path)
  progress$inc(1)

  # get dprime effect size values for analysis
  dprimes <- get_dprimes(top_table)

  # exclude ambient (for single cell only)
  dprimes <- dprimes[!names(dprimes) %in% ambient]

  # get correlations between query and drug signatures
  res$cmap <- query_drugs(dprimes, cmap_es)
  res$l1000_drugs <- query_drugs(dprimes, l1000_drugs_es)
  res$l1000_genes <- query_drugs(dprimes, l1000_genes_es)
  progress$inc(1)

  saveRDS(res$cmap, res_paths$cmap)
  saveRDS(res$l1000_drugs, res_paths$l1000_drugs)
  saveRDS(res$l1000_genes, res_paths$l1000_genes)

  return(res)
}

#' Run single cell RNA-Seq test vs ctrl cluster comparison
#'
#' Used by \code{run_comparison} for differential expression analysis.
#'
#' @param scseq \code{Suerat} object
#' @param selected_clusters Character vector of selected clusters to run comparison for.
#' @param sc_dir Path to folder with single cell analyses.
#' @param anal_name Folder name in \code{sc_dir} that contains single cell analysis.
#' @param res_paths List of paths to saved results for markers of cells in \code{selected_clusters} and
#'  differential expression analysis comparing test to control cells in \code{selected_clusters}.
#'  Paths are stored in \code{'markers'} and \code{'anal'} slots.
#'
#' @return List with slots \code{'markers'} and \code{'anal'} containing results for markers of
#'  cells in \code{selected_clusters} and differential expression analysis comparing test to control
#'  cells in \code{selected_clusters}.
run_cluster_comparison <- function(scseq, selected_clusters, sc_dir, anal_name) {

  clusters_name <- collapse_sorted(selected_clusters)
  clusters <- as.character(Seurat::Idents(scseq))
  in.sel <- clusters %in% selected_clusters

  new.idents <- clusters
  new.idents[in.sel] <- 'ident.1'
  Seurat::Idents(scseq) <- factor(new.idents)

  cluster_markers <- get_scseq_markers(scseq, ident.1 = 'ident.1')

  fname <- paste0("markers_", clusters_name, '.rds')
  fpath <- file.path(sc_dir, anal_name, fname)
  saveRDS(cluster_markers, fpath)

  Seurat::Idents(scseq) <- scseq$orig.ident
  scseq <- scseq[, in.sel]
  fit <- fit_lm_scseq(scseq = scseq,
                          data_dir = sc_dir,
                          anal_name = anal_name,
                          clusters_name = clusters_name)

  res <- list(
    cluster_markers = cluster_markers,
    fit = fit
  )
  return(res)
}

#' Run single cell RNA-Seq test vs ctrl pathway comparison
#'
#' Used by \code{run_comparison} for pathway analysis.
#'
#' @param scseq \code{Suerat} object
#' @param selected_clusters Character vector of selected clusters to run comparison for.
#' @param sc_dir Path to folder with single cell analyses.
#' @param anal_name Folder name in \code{sc_dir} that contains single cell analysis.
#' @param res List with previous differential expression analysis results in slot \code{'anal'}.
#' @param ambient Character vector of ambient genes to exclude from pathway analysis.
#'
#' @return \code{res} with slot \code{'path'} containing pathway analysis results. The \code{'anal'}
#'  slot is also subsetted such that ambient genes are excluded.
run_path_comparison <- function(scseq, selected_clusters, sc_dir, anal_name, res, ambient) {

  clusters_name <- collapse_sorted(selected_clusters)
  clusters <- as.character(Seurat::Idents(scseq))
  in.sel <- clusters %in% selected_clusters

  # run pathway analysis (will load if exists)


  Seurat::Idents(scseq) <- scseq$orig.ident
  scseq <- scseq[, in.sel]
  res$path <- diff_path_scseq(scseq,
                              prev_anal = anal,
                              ambient = ambient,
                              data_dir = sc_dir,
                              anal_name = anal_name,
                              clusters_name = clusters_name,
                              NI = 24)

  # remove ambient from markers
  is.ambient <-row.names(res$anal$top_table) %in% ambient

  res$anal$top_table <- res$anal$top_table[!is.ambient, ]
  res$anal$ebayes_sv$df.residual <- res$anal$ebayes_sv$df.residual[!is.ambient]

  return(res)
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
#' @param scseq Integrated \code{Seurat} object.
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
  cols <- get_palette(levels(scseq$seurat_clusters))
  names(cols) <- levels(scseq$seurat_clusters)

  # median coordinates for clusters for integrated dataset
  plot_data <- plot$data
  int_coords <- get_umap_coords(plot_data)

  # median coordinates for clusters for individual dataset
  in_anal <- scseq$project %in% anal
  scseq <- scseq[, in_anal]
  plot_data <- plot_data[colnames(scseq), ]
  plot_data$ident <- scseq$orig_clusters
  anal_coords <- get_umap_coords(plot_data)

  # closest match for each anal cluster
  match <- get_nearest_row(int_coords[ ,-1], anal_coords[, -1])
  anal_coords <- anal_coords %>%
    mutate(match = int_coords$ident[match]) %>%
    arrange(as.numeric(ident))

  # change identity to previous labels
  cl <- scseq$orig_clusters
  lv <- as.character(sort(unique(as.numeric(cl))))
  names(annot) <- 0:(length(annot)-1)

  scseq$seurat_clusters <- factor(cl, levels = lv)
  levels(scseq$seurat_clusters) <- annot[lv]
  Seurat::Idents(scseq) <- scseq$seurat_clusters

  # expand to best match
  anal_coords$cols <- cols[anal_coords$match]

  # replace any duplicates
  is.dup <- which(duplicated(anal_coords$cols))
  if (length(is.dup)) {
    new_cols <- setdiff(get_palette(lv), anal_coords$cols)
    anal_coords[is.dup, 'cols'] <- new_cols[seq_along(is.dup)]
  }

  plot_umap_cluster(scseq, cols=anal_coords$cols)
}

#' Get median x-y coordinates for clusters in UMAP plot data
#'
#' @param plot_data data.frame with columns \code{'ident'}, \code{'UMAP_1'}, and \code{'UMAP_2'}.
get_umap_coords <- function(plot_data) {
  plot_data %>%
    group_by(ident) %>%
    summarise(UMAP_1 = median(UMAP_1),
              UMAP_2 = median(UMAP_2))
}

