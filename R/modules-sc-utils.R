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
get_exclude_choices <- function(anal_names, anal_colors, data_dir) {

  if (is.null(anal_names)) return(NULL)

  # load markers and annotation for each
  annot_paths <- scseq_part_path(data_dir, anal_names, 'annot')
  marker_paths <- scseq_part_path(data_dir, anal_names, 'markers')

  annots <- lapply(annot_paths, readRDS)
  clusters <- lapply(annots, function(x) seq(0, length(x)-1))
  colors <- lapply(annots, get_palette)

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
#'
#' @return data.frame with columns for rendering selectizeInput cluster choices
#' @export
#' @keywords internal
get_cluster_choices <- function(clusters, scseq, value = clusters) {

  # show the cell numbers/percentages
  ncells <- tabulate(scseq$seurat_clusters)
  pcells <- round(ncells / sum(ncells) * 100)
  pspace <- strrep('&nbsp;&nbsp;', 2 - nchar(pcells))

  # cluster choices are the clusters themselves
  testColor <- get_palette(clusters)
  cluster_choices <- data.frame(name = stringr::str_trunc(clusters, 27),
                                value = seq(0, along.with = clusters),
                                label = clusters,
                                testColor,
                                ncells, pcells, pspace, row.names = NULL, stringsAsFactors = FALSE)

  return(cluster_choices)
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
  sct_paths <- scseq_part_path(data_dir, anal_names, 'sct')
  scseq_paths <- scseq_part_path(data_dir, anal_names, 'sct')

  exclude_anals <- gsub('^(.+?)_\\d+$', '\\1', exclude_clusters)
  exclude_clusters <- gsub('^.+?_(\\d+)$', '\\1', exclude_clusters)

  scseqs <- list()
  for (anal in anal_names) {

    # load scseq
    scseq_path <- scseq_part_path(data_dir, anal, 'scseq')
    scseq <- readRDS(scseq_path)

    # add annotation
    annot_path <- scseq_part_path(data_dir, anal, 'annot')
    annot <- readRDS(annot_path)

    scseq$annot_clusters <- scseq$seurat_clusters
    levels(scseq$annot_clusters) <- annot

    # restore SCT assay
    sct_path <- scseq_part_path(data_dir, anal, 'sct')
    if (file.exists(sct_path)) {
      sct <- readRDS(sct_path)
      scseq[['SCT']] <- sct
    }

    # restore counts slot
    scseq[['RNA']]@counts <- scseq[['RNA']]@data

    # set orig.ident to ctrl/test
    scseq$orig.ident <- factor(ident)

    # only remove excluded clusters if present
    is.exclude <- exclude_anals == anal
    if (any(is.exclude)) {
      exclude <- exclude_clusters[is.exclude]
      scseq <- scseq[, !scseq$seurat_clusters %in% exclude]

    }

    # downsample very large datasets
    # scseq <- downsample_scseq(scseq)

    # bug in Seurat, need for integration
    scseq[['SCT']]@misc$vst.out$cell_attr <- scseq[['SCT']]@misc$vst.out$cell_attr[colnames(scseq), ]

    # add to scseqs
    scseqs[[anal]] <- scseq
  }

  return(scseqs)
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
#' @return \code{NULL} is valid, otherwise an error message
#' @export
#' @keywords internal
validate_integration <- function(test, ctrl, anal_name, anal_options) {
  msg <- NULL
  # make sure both control and test analyses provided
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
#' @return
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


#' Run genes and pathways test vs ctrl scseq comparisons
#'
#' Used in Single Cell tab (no pathway analysis) and Pathways tab (with pathway analysis).
#'
#' @param scseq \code{Suerat} object
#' @param selected_clusters the selected clusters in \code{Seurat::Idents(scseq)}
#' @param sc_dir Path to directory with single cell analysis folders.
#' @param anal_name Name of analysis. A directory in \code{sc_dir}.
#' @param with_path Boolean to include pathway analysis or not.
#'
#' @return Named list with names. \code{anal} \code{cluster_markers} and optionally \code{path}.
#' @export
#' @keywords internal
run_comparison <- function(scseq, selected_clusters, sc_dir, anal_name, with_path = FALSE) {


  clusters_name <- collapse_sorted(selected_clusters)

  markers_path <- scseq_part_path(sc_dir, anal_name, paste0('markers_', clusters_name))
  anal_path <- scseq_part_path(sc_dir, anal_name, paste0('diff_expr_symbol_scseq_', clusters_name))

  # get markers for selected cluster(s)
  # so that don't exclude marker genes as ambient
  clusters <- as.character(Seurat::Idents(scseq))
  in.sel <- clusters %in% selected_clusters

  if (file.exists(markers_path)) {
    cluster_markers <- readRDS(markers_path)

  } else {

    new.idents <- clusters
    new.idents[in.sel] <- 'ident.1'
    Seurat::Idents(scseq) <- factor(new.idents)

    cluster_markers <- get_scseq_markers(scseq, ident.1 = 'ident.1')
    saveRDS(cluster_markers, markers_path)
  }

  # get markers for test group
  if (file.exists(anal_path)) {
    anal <- readRDS(anal_path)

  } else {

    Seurat::Idents(scseq) <- scseq$orig.ident
    scseq <- scseq[, in.sel]
    anal <- diff_expr_scseq(scseq = scseq,
                            data_dir = sc_dir,
                            anal_name = anal_name,
                            clusters_name = clusters_name)
  }

  res <- list(
    cluster_markers = cluster_markers,
    anal = anal)

  if (with_path) {
    # run pathway analysis (will load if exists)
    ambient <- get_ambient(scseq, markers = anal$top_table, cluster_markers = cluster_markers)

    Seurat::Idents(scseq) <- scseq$orig.ident
    scseq <- scseq[, in.sel]
    res$path <- diff_path_scseq(scseq,
                                prev_anal = anal,
                                ambient = ambient,
                                data_dir = sc_dir,
                                anal_name = anal_name,
                                clusters_name = clusters_name)

    # remove ambient from markers
    is.ambient <-row.names(anal$top_table) %in% ambient

    res$anal$top_table <- anal$top_table[!is.ambient, ]
    res$anal$ebayes_sv$df.residual <- anal$ebayes_sv$df.residual[!is.ambient]

  }

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

