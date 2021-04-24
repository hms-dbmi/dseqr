#' Process raw single cell fastq or cellranger files for app
#'
#' @param dataset_name Name of dataset
#' @param fastq_dir Directory with fastq or cellranger files
#' @param sc_dir Single cell directory for app. Will store results in \code{dataset_name} subdirectory
#' @param progress Optional shiny \code{Progress} object. Default will print progress.
#' @param value Integer indicating step of pipeline.
#' @param founder Name of dataset that \code{dataset_name} originates from.
#' @inheritParams run_dseqr
#' @inheritParams run_kallisto_scseq
#'
#' @return NULL
#' @keywords internal
#'
load_raw_scseq <- function(dataset_name,
                           fastq_dir,
                           sc_dir,
                           indices_dir,
                           progress = NULL,
                           recount = FALSE,
                           value = 0,
                           founder = dataset_name,
                           npcs = 30,
                           cluster_alg = 'leiden',
                           resoln = 1,
                           azimuth_ref = NULL,
                           metrics = c('low_lib_size',
                                       'low_n_features',
                                       'high_subsets_mito_percent',
                                       'low_subsets_ribo_percent',
                                       'high_doublet_score')) {
  if (is.null(progress)) {
    progress <- list(set = function(value, message = '', detail = '') {
      cat(value, message, detail, '...\n')
    })
  }

  # standardize cellranger files
  is.cellranger <- check_is_cellranger(fastq_dir)

  progress$set(message = "quantifying files", value = value + 1)
  if (!is.cellranger) run_kallisto_scseq(indices_dir, fastq_dir, recount = recount)

  progress$set(message = "loading", value + 2)
  type <- ifelse(is.cellranger, 'cellranger', 'kallisto')
  scseq <- create_scseq(fastq_dir, project = dataset_name, type = type)
  gc()

  progress$set(message = "running QC", value = value + 3)
  scseq <- normalize_scseq(scseq)
  scseq <- add_hvgs(scseq)
  gc()

  scseq <- add_doublet_score(scseq)
  scseq <- add_scseq_qc_metrics(scseq, scseq@metadata$species, for_qcplots = TRUE)

  scseq <- run_scseq_qc(scseq, metrics)

  process_raw_scseq(scseq,
                    dataset_name,
                    sc_dir,
                    cluster_alg = cluster_alg,
                    npcs = npcs,
                    resoln = resoln,
                    azimuth_ref = azimuth_ref,
                    founder = founder,
                    progress = progress,
                    value = value + 3)
}


#' Convenience utility to run load_raw_scseq in background

#' @keywords internal
#' @noRd
run_load_raw_scseq <- function(opts, fastq_dir, sc_dir, indices_dir, azimuth_ref = NULL) {

  for (opt in opts) {
    load_raw_scseq(opt$dataset_name,
                   fastq_dir,
                   sc_dir,
                   indices_dir,
                   metrics = opt$metrics,
                   founder = opt$founder,
                   azimuth_ref = azimuth_ref)

  }

  return(TRUE)
}


#' Process Count Data for App
#'
#' @param scseq \code{SingleCellExperiment}
#' @param dataset_name Name of dataset to save
#' @param sc_dir Directory to save dataset to
#' @param progress Shiny progress object. Default (\code{NULL}) prints to stdout.
#' @param value Initial value of progress.
#' @inheritParams subset_saved_scseq
#'
#' @return NULL
#' @export
process_raw_scseq <- function(scseq,
                              dataset_name,
                              sc_dir,
                              cluster_alg = 'leiden',
                              npcs = 30,
                              resoln = 1,
                              hvgs = NULL,
                              azimuth_ref = NULL,
                              founder = NULL,
                              progress = NULL,
                              value = 0) {

  if (is.null(progress)) {
    progress <- list(set = function(value, message = '', detail = '') {
      cat(value, message, detail, '...\n')
    })
  }


  scseq@metadata$npcs <- npcs
  # TODO: can this be skipped when not subset?
  scseq <- normalize_scseq(scseq)
  scseq <- add_hvgs(scseq, hvgs = hvgs)

  is.azimuth <- !is.null(azimuth_ref)
  if (!is.azimuth) {
    progress$set(message = "reducing", value = value + 1)
    scseq <- run_pca(scseq)
    scseq <- run_reduction(scseq)
    gc()

    progress$set(message = "clustering", detail = '', value = value + 2)
    snn_graph <- get_snn_graph(scseq, npcs)
    scseq$cluster <- get_clusters(snn_graph, cluster_alg, resoln)
    gc()

    # save independent of resolution
    anal <- list(scseq = scseq, snn_graph = snn_graph, founder = founder, resoln = resoln)
    save_scseq_data(anal, dataset_name, sc_dir)

  } else {
    progress$set(message = "running Azimuth", detail = '', value = value + 2)
    azres <- run_azimuth(list(one = scseq), azimuth_ref)
    scseq <- transfer_azimuth(azres, scseq)
    rm(azres); gc()

    # because e.g. default resoln is celltype.l2 for human pbmc
    resoln <- switch(azimuth_ref, 'human_pbmc' = 2)

    anal <- list(scseq = scseq, founder = founder, resoln = resoln, azimuth_ref = azimuth_ref)
    save_scseq_data(anal, dataset_name, sc_dir)
    save_azimuth_clusters(scseq@colData, dataset_name, sc_dir)
  }

  progress$set(message = "saving loom", value = value + 3)
  save_scle(scseq, file.path(sc_dir, dataset_name))

  # run what depends on resolution
  run_post_cluster(scseq,
                   dataset_name,
                   sc_dir,
                   resoln,
                   progress,
                   value + 3,
                   reset_annot = !is.azimuth)

  progress$set(value = value + 7)
}

transfer_azimuth <- function(azres, scseq) {
  # add new data back to scseq
  projs <- lapply(azres, function(x) {
    proj <- x@reductions$proj.umap@cell.embeddings
    colnames(proj) <- gsub('_', '', colnames(proj))
    return(proj)
  })

  proj <- do.call(rbind, projs)
  type <- gsub('1$', '', colnames(proj)[1])
  SingleCellExperiment::reducedDim(scseq, type) <- proj

  get_meta <- function(l, m) unlist(lapply(l, function(x) x@meta.data[[m]]), use.names = FALSE)

  scseq$mapping_score <- get_meta(azres, 'mapping.score')
  scseq$cluster_l1_score <- get_meta(azres, 'predicted.celltype.l1.score')
  scseq$cluster_l2_score <- get_meta(azres, 'predicted.celltype.l2.score')
  scseq$cluster_l3_score <- get_meta(azres, 'predicted.celltype.l3.score')

  scseq$cluster_l1 <- get_meta(azres, 'predicted.celltype.l1')
  scseq$cluster_l2 <- get_meta(azres, 'predicted.celltype.l2')
  scseq$cluster_l3 <- get_meta(azres, 'predicted.celltype.l3')

  clus <- factor(scseq$cluster_l2)
  scseq$cluster <- factor(as.numeric(clus))
  return(scseq)
}


run_azimuth <- function(scseqs, azimuth_ref) {

  reference <- dseqr.data::load_data(paste0(azimuth_ref, '.qs'))

  refdata <- lapply(c("celltype.l1", "celltype.l2", "celltype.l3"), function(x) {
    reference$map[[x, drop = TRUE]]
  })

  names(refdata) <- c("celltype.l1", "celltype.l2", "celltype.l3")

  refdata[["impADT"]] <- Seurat::GetAssayData(
    object = reference$map[['ADT']],
    slot = 'data'
  )

  queries <- list()
  for (ds in names(scseqs)) {
    scseq <- scseqs[[ds]]
    counts <- SingleCellExperiment::counts(scseq)
    query <- Seurat::CreateSeuratObject(counts = counts,
                                        min.cells = 1, min.features = 1)

    # Preprocess with SCTransform
    query <- Seurat::SCTransform(
      object = query,
      assay = "RNA",
      new.assay.name = "refAssay",
      residual.features = rownames(reference$map),
      reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
      method = 'glmGamPoi',
      ncells = 2000,
      n_genes = 2000,
      do.correct.umi = FALSE,
      do.scale = FALSE,
      do.center = TRUE,
      verbose = FALSE
    )

    # Find anchors between query and reference
    anchors <- Seurat::FindTransferAnchors(
      reference = reference$map,
      query = query,
      k.filter = NA,
      reference.neighbors = "refdr.annoy.neighbors",
      reference.assay = "refAssay",
      query.assay = "refAssay",
      reference.reduction = "refDR",
      normalization.method = "SCT",
      features = intersect(rownames(reference$map), Seurat::VariableFeatures(query)),
      dims = 1:50,
      n.trees = 20,
      mapping.score.k = 100,
      verbose = FALSE
    )

    # Transfer cell type labels and impute protein expression
    #
    # Transferred labels are in metadata columns named "predicted.*"
    # The maximum prediction score is in a metadata column named "predicted.*.score"
    # The prediction scores for each class are in an assay named "prediction.score.*"
    # The imputed assay is named "impADT" if computed
    query <- Seurat::TransferData(
      reference = reference$map,
      query = query,
      dims = 1:50,
      anchorset = anchors,
      refdata = refdata,
      n.trees = 20,
      store.weights = TRUE,
      verbose = FALSE
    )

    # Calculate the embeddings of the query data on the reference SPCA
    query <- Seurat::IntegrateEmbeddings(
      anchorset = anchors,
      reference = reference$map,
      query = query,
      reductions = "pcaproject",
      reuse.weights.matrix = TRUE,
      verbose = FALSE
    )

    # Calculate the query neighbors in the reference
    # with respect to the integrated embeddings
    query[["query_ref.nn"]] <- Seurat::FindNeighbors(
      object = Seurat::Embeddings(reference$map[["refDR"]]),
      query = Seurat::Embeddings(query[["integrated_dr"]]),
      return.neighbor = TRUE,
      l2.norm = TRUE,
      verbose = FALSE
    )


    # The reference used in the app is downsampled compared to the reference on which
    # the UMAP model was computed. This step, using the helper function NNTransform,
    # corrects the Neighbors to account for the downsampling.
    query <- NNTransform(
      object = query,
      meta.data = reference$map[[]]
    )

    # Project the query to the reference UMAP.
    query[["proj.umap"]] <- Seurat::RunUMAP(
      object = query[["query_ref.nn"]],
      reduction.model = reference$map[["refUMAP"]],
      reduction.key = 'UMAP_',
      verbose = FALSE
    )

    # Calculate mapping score and add to metadata
    query <- Seurat::AddMetaData(
      object = query,
      metadata = Seurat::MappingScore(anchors = anchors),
      col.name = "mapping.score"
    )


    queries[[ds]] <- query
  }

  return(queries)
}

save_azimuth_clusters <- function(meta, dataset_name, sc_dir) {

  # cluster columns
  cols <- grep('^cluster_l\\d$', colnames(meta), value = TRUE)

  # save annotation and clusters for each
  for (i in seq_along(cols)) {
    dataset_subname <- file.path(dataset_name, paste0('snn', i))
    col <- cols[i]
    clusters <- factor(meta[[col]])

    scseq_data <- list(
      annot = levels(clusters),
      clusters = factor(as.numeric(clusters))
    )

    save_scseq_data(scseq_data, dataset_subname, sc_dir)
  }
}





#' Load kallisto/bustools quantification into a SingleCellExperiment object.
#'
#' @param data_dir Directory with raw and kallisto/bustools or CellRanger quantified single-cell RNA-Seq files.
#' @param project String identifying sample.
#' @param type Quantification file type. One of either \code{'kallisto'} or \code{'cellranger'}.
#'
#' @return \code{SingleCellExperiment} object with empty droplets removed and ambient outliers recorded.
#' @export
#'
create_scseq <- function(data_dir, project, type = c('kallisto', 'cellranger')) {

  # load counts
  #TODO suport mouse for kallisto
  if (type[1] == 'kallisto') {
    data_dir <- file.path(data_dir, 'bus_output')
    counts <- load_kallisto_counts(data_dir)
    species <- 'Homo sapiens'

  } else if (type[1] == 'cellranger') {
    counts <- load_cellranger_counts(data_dir)
    species <- get_species(counts)

    counts <- process_cellranger_counts(counts, species)
  }

  # get ambience expression profile/determine outlier genes
  # if pre-filtered cellranger, can't determine outliers/empty droplets
  ncount <- Matrix::colSums(counts)
  if (min(ncount) > 10) {
    ambience <- rep(0, nrow(counts))
    keep_cells <- seq_len(ncol(counts))

  } else {
    ambience <- DropletUtils::estimateAmbience(counts, good.turing = FALSE, round = FALSE)
    keep_cells <- detect_cells(counts, species = species)
  }

  # add ambience metadata for genes
  rowData <- S4Vectors::DataFrame(ambience)
  colData <- S4Vectors::DataFrame(project = rep(project, length(keep_cells)))

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts[, keep_cells]),
    rowData = rowData,
    colData = colData
  )

  # flag to know if need to convert for drug queries
  sce@metadata$species <- species
  return(sce)
}



#' Load SingleCellExperiment
#'
#' Loads scle.loom if exists otherwise scseq.qs
#'
#' @param dataset_dir Directory with scle.loom or scseq.qs
#'
#' @return \code{SingleCellExperiment} or \code{SingleCellLoomExperiment}
#' @export
#' @keywords internal
load_scseq <- function(dataset_dir, default_clusters = TRUE) {
  scle_path <- file.path(dataset_dir, 'scle.loom')

  transition_efs(scle_path)

  # load loom if available (faster and less memory)
  scseq <- tryCatch(LoomExperiment::import(scle_path, type = 'SingleCellLoomExperiment'),
                    error = function(e) {
                      unlink(scle_path)
                      return(load_scseq_qs(dataset_dir))
                    })


  if (default_clusters) {
    resoln_dir <- file.path(dataset_dir, load_resoln(dataset_dir))
    scseq <- attach_clusters(scseq, resoln_dir)
  }

  scseq <- attach_meta(scseq, dataset_dir)

  return(scseq)
}


attach_clusters <- function(scseq, resoln_dir) {
  clusters_path <- file.path(resoln_dir, 'clusters.qs')
  scseq$cluster <- qs::qread(clusters_path)
  return(scseq)
}

attach_meta <- function(scseq, dataset_dir = NULL, meta = NULL, groups = NULL) {

  if (is.null(meta)) {
    meta_path <- file.path(dataset_dir, 'meta.qs')
    meta <- qread.safe(meta_path)
  }

  if (is.null(groups)) {
    group_path <- file.path(dataset_dir, 'prev_groups.qs')
    groups <- qread.safe(group_path)
  }

  if (is.null(meta)) return(scseq)

  meta <- meta[meta$group %in% groups, ]

  groups <- ifelse(meta$group == groups[1], 'test', 'ctrl')
  names(groups) <- row.names(meta)
  groups <- unname(groups[scseq$batch])

  scseq$orig.ident <- factor(groups, levels = c('test', 'ctrl'))
  return(scseq)
}

#' Utility to move files on EFS from IA to SA
#'
#' If time since last access time is greater than \code{Sys.getenv('EFS_LIFECYCLE')},
#' \code{fpath} to copied in order to move it out of EFS infrequent access.
#' Ignored if EFS_LIFECYCLE not set.
#'
#' @param fpath Path of file
#'
#' @return NULL
#'
transition_efs <- function(fpath) {

  # get efs lifecycle from ENV
  efs_diff <- Sys.getenv('EFS_LIFECYCLE')
  efs_diff <- as.difftime(as.numeric(efs_diff), units='days')
  if (is.na(efs_diff)) return(NULL)


  # time since last read
  read_diff <- Sys.time() - file.info(fpath)$atime
  if (read_diff < efs_diff) return(NULL)

  # move out of IA by copying
  message('EFS_LIFECYLCE is: ', efs_diff, ' days.')
  message('moving: ', fpath, ' to standard access.')
  tmp <- file.path(dirname(fpath), 'tmp.loom')
  file.copy(fpath, tmp)
  unlink(fpath)
  file.move(tmp, fpath)
}


#' Save SingleCellExperiment as loom file
#'
#' @param scseq \code{SingleCellExperiment}
#' @param dataset_dir Directory to save in
#' @param overwrite Overwrite existing? Default is \code{TRUE}
#'
#' @return NULL
#' @export
#' @keywords internal
#'
save_scle <- function(scseq, dataset_dir, overwrite = TRUE) {
  # don't save raw counts for loom (saved in non-sparse format)
  SummarizedExperiment::assay(scseq, 'counts') <- NULL

  scle_path <- file.path(dataset_dir, 'scle.loom')


  if (!file.exists(scle_path) | overwrite) {
    unlink(scle_path)

    if ('corrected' %in% SingleCellExperiment::reducedDimNames(scseq)) {

      SingleCellExperiment::reducedDim(scseq, 'corrected') <-
        as.matrix(SingleCellExperiment::reducedDim(scseq, 'corrected'))
    }

    # these objects from fastMNN won't save
    scseq@metadata <- scseq@metadata[!names(scseq@metadata) %in% c('cluster', 'merge.info')]

    scseq <- LoomExperiment::SingleCellLoomExperiment(scseq)
    LoomExperiment::export(scseq, scle_path)
  }
}


#' Read kallisto/bustools market matrix and annotations
#'
#'
#' @param data_dir Path to folder with 'genecount' directory which contains
#'   'genes.mtx', 'genes.genes.txt', and 'genes.barcodes.txt' files generated
#'   by kallisto/bustools.
#'
#' @return sparse dgTMatrix with barcodes in columns and genes in rows.
#' @keywords internal
load_kallisto_counts <- function(data_dir) {

  # read sparse matrix
  counts <- Matrix::readMM(file.path(data_dir, 'genecount', 'genes.mtx'))
  counts <- Matrix::t(counts)
  counts <- methods::as(counts, 'dgCMatrix')

  # read annotations
  row.names(counts) <- readLines(file.path(data_dir, 'genecount', 'genes.genes.txt'))
  colnames(counts) <- readLines(file.path(data_dir, 'genecount', 'genes.barcodes.txt'))

  # remove non-expressed genes
  counts <- counts[Matrix::rowSums(counts) > 0, ]

  return(counts)
}

#' Load cell ranger counts
#'
#' Mainly to avoid having to download massive datasets that have already been
#' quantified.
#'
#' @param data_dir Path to folder with cell ranger files.
#' @importFrom magrittr "%>%"
#' @keywords internal
#'
#' @return dgCMatrix
load_cellranger_counts <- function(data_dir) {
  # read the data in using ENSG features
  h5file <- list.files(data_dir, '.h5$', full.names = TRUE)
  if (length(h5file)) {
    counts <- Read10X_h5(h5file, use.names = FALSE)
  } else {
    counts <- Read10X(data_dir, gene.column = 1)
  }

  if (methods::is(counts, 'list')) counts <- counts$`Gene Expression`
  return(counts)
}

process_cellranger_counts <- function(counts, species) {
  gene_id <- gene_name <- NULL

  if (!grepl('sapiens|musculus', species)) stop('Species not supported')
  tx2gene <- dseqr.data::load_tx2gene(species)

  counts <- counts[row.names(counts) %in% tx2gene$gene_id, ]

  # name genes by tx2gene so that best match with cmap/l1000 data
  map <- tx2gene %>%
    dplyr::filter(gene_id %in% row.names(counts)) %>%
    dplyr::select(gene_id, gene_name) %>%
    dplyr::distinct() %>%
    dplyr::arrange(match(gene_id, row.names(counts)))

  stopifnot(setequal(row.names(counts), map$gene_id))

  row.names(counts) <- map$gene_name

  # remove non-expressed genes
  counts <- counts[Matrix::rowSums(counts) > 0, ]

  # sum counts in rows with same gene
  counts <- Matrix.utils::aggregate.Matrix(counts, row.names(counts), fun = 'sum')

  return(counts)
}

#' Get tx2gene for human or mouse depending on intersections
#'
#' @param counts dgCMatrix where rownames are gene ids
#'
#' @return tx2gene for human or mouse
#' @export
#' @keywords internal
#'
get_species <- function(counts) {

  tx2gene <- dseqr.data::load_tx2gene('Homo sapiens')
  tx2gene_mouse <- dseqr.data::load_tx2gene('Mus musculus')

  nhuman <- sum(row.names(counts) %in% tx2gene$gene_id)
  nmouse <- sum(row.names(counts) %in% tx2gene_mouse$gene_id)

  if (nhuman > nmouse) return('Homo sapiens')
  else if (nmouse > nhuman) return('Mus musculus')
}

#' Load in data from 10X
#'
#' Enables easy loading of sparse data matrices provided by 10X genomics.
#'
#' @param data.dir Directory containing the matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv
#' files provided by 10X. A vector or named vector can be given in order to load
#' several data directories. If a named vector is given, the cell barcode names
#' will be prefixed with the name.
#' @param gene.column Specify which column of genes.tsv or features.tsv to use for gene names; default is 2
#' @param unique.features Make feature names unique (default TRUE)
#'
#' @return If features.csv indicates the data has multiple data types, a list
#'   containing a sparse matrix of the data from each type will be returned.
#'   Otherwise a sparse matrix containing the expression data will be returned.
#'
#' @keywords internal
#'
Read10X <- function(data.dir = NULL, gene.column = 2, unique.features = TRUE) {
  full.data <- list()
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, 'barcodes.tsv')
    gene.loc <- file.path(run, 'genes.tsv')
    features.loc <- file.path(run, 'features.tsv.gz')
    matrix.loc <- file.path(run, 'matrix.mtx')
    # Flag to indicate if this data is from CellRanger >= 3.0
    pre_ver_3 <- file.exists(gene.loc)
    if (!pre_ver_3) {
      addgz <- function(s) {
        return(paste0(s, ".gz"))
      }
      barcode.loc <- addgz(s = barcode.loc)
      matrix.loc <- addgz(s = matrix.loc)
    }
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing")
    }
    if (!pre_ver_3 && !file.exists(features.loc) ) {
      stop("Gene name or features file missing")
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing")
    }
    data <- Matrix::readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names))) {
      cell.names <- as.vector(x = as.character(x = sapply(
        X = cell.names,
        FUN = ExtractField,
        field = 1,
        delim = "-"
      )))
    }
    if (is.null(x = names(x = data.dir))) {
      if (i < 2) {
        colnames(x = data) <- cell.names
      } else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    } else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], "_", cell.names)
    }
    feature.names <- utils::read.delim(
      file = ifelse(test = pre_ver_3, yes = gene.loc, no = features.loc),
      header = FALSE,
      stringsAsFactors = FALSE
    )
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning(
        'Some features names are NA. Replacing NA names with ID from the opposite column requested',
        call. = FALSE,
        immediate. = TRUE
      )
      na.features <- which(x = is.na(x = feature.names[, gene.column]))
      replacement.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[na.features, replacement.column]
    }
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ", gene.column,
                    " but feature.tsv.gz (or genes.tsv) only has ", fcols, " columns.",
                    " Try setting the gene.column argument to a value <= to ", fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, gene.column])
    }
    # In cell ranger 3.0, a third column specifying the type of data was added
    # and we will return each type of data as a separate matrix
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) { # Return Gene Expression first
        lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
      }
      data <- lapply(
        X = lvls,
        FUN = function(l) {
          return(data[data_types == l, ])
        }
      )
      names(x = data) <- lvls
    } else{
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  # Combine all the data from different directories into one big matrix, note this
  # assumes that all data directories essentially have the same features files
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, FUN = `[[`, j))
    # Fix for Issue #913
    list_of_data[[j]] <- methods::as(object = list_of_data[[j]], Class = "dgCMatrix")
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  # If multiple features, will return a list, otherwise
  # a matrix.
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  } else {
    return(list_of_data)
  }
}

# Extract delimiter information from a string.
#
# Parses a string (usually a cell name) and extracts fields based on a delimiter
#
# @param string String to parse.
# @param field Integer(s) indicating which field(s) to extract. Can be a vector multiple numbers.
# @param delim Delimiter to use, set to underscore by default.
#
# @return A new string, that parses out the requested fields, and (if multiple), rejoins them with the same delimiter
#
# @export
#
# @examples
# ExtractField(string = 'Hello World', field = 1, delim = '_')
#
ExtractField <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), split = ",")))
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
}

#' Read 10X hdf5 file
#'
#' Read count matrix from 10X CellRanger hdf5 file.
#' This can be used to read both scATAC-seq and scRNA-seq matrices.
#'
#' @param filename Path to h5 file
#' @param use.names Label row names with feature names rather than ID numbers.
#' @param unique.features Make feature names unique (default TRUE)
#'
#' @return Returns a sparse matrix with rows and columns labeled. If multiple
#' genomes are present, returns a list of sparse matrices (one per genome).
#'
#' @keywords internal
#'
Read10X_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
  genomes <- names(x = infile)
  output <- list()
  if (!infile$attr_exists("PYTABLES_FORMAT_VERSION")) {
    # cellranger version 3
    if (use.names) {
      feature_slot <- 'features/name'
    } else {
      feature_slot <- 'features/id'
    }
  } else {
    if (use.names) {
      feature_slot <- 'gene_names'
    } else {
      feature_slot <- 'genes'
    }
  }
  for (genome in genomes) {
    counts <- infile[[paste0(genome, '/data')]]
    indices <- infile[[paste0(genome, '/indices')]]
    indptr <- infile[[paste0(genome, '/indptr')]]
    shp <- infile[[paste0(genome, '/shape')]]
    features <- infile[[paste0(genome, '/', feature_slot)]][]
    barcodes <- infile[[paste0(genome, '/barcodes')]]
    sparse.mat <- Matrix::sparseMatrix(
      i = indices[] + 1,
      p = indptr[],
      x = as.numeric(x = counts[]),
      dims = shp[],
      repr = "T"
    )
    if (unique.features) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- methods::as(object = sparse.mat, Class = 'dgCMatrix')
    # Split v3 multimodal
    if (infile$exists(name = paste0(genome, '/features'))) {
      types <- infile[[paste0(genome, '/features/feature_type')]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
        sparse.mat <- sapply(
          X = types.unique,
          FUN = function(x) {
            return(sparse.mat[which(x = types == x), ])
          },
          simplify = FALSE,
          USE.NAMES = TRUE
        )
      }
    }
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  if (length(x = output) == 1) {
    return(output[[genome]])
  } else{
    return(output)
  }
}

#' Determine if the selected folder has CellRanger files
#'
#' @param data_dir path to directory to check.
#'
#' @return \code{TRUE} if CellRanger files detected, otherwise \code{FALSE}.
#' @keywords internal
check_is_cellranger <- function(data_dir) {

  # cellranger file names
  files <- list.files(data_dir)
  mtx.file <- grep('.mtx', files, fixed = TRUE, value = TRUE)
  genes.file <- grep('features.tsv|genes.tsv', files, value = TRUE)
  barcodes.file <-  grep('barcodes.tsv', files, fixed = TRUE, value = TRUE)

  # h5 file
  h5.file <- grepl('.h5', files, fixed = TRUE)

  if (length(mtx.file) & length(genes.file) & length(barcodes.file)) {
    standardize_cellranger(data_dir)
    return(TRUE)
  }

  if (h5.file) return(TRUE)
  return(FALSE)
}

#' Rename CellRanger files for loading by Read10X
#'
#' @param data_dir Path to folder with CellRanger files.
#'
#' @return NULL
#' @keywords internal
standardize_cellranger <- function(data_dir) {

  # cellranger file names
  files <- list.files(data_dir)
  mtx.file <- grep('.mtx', files, fixed = TRUE, value = TRUE)
  genes.file <- grep('features.tsv|genes.tsv', files, value = TRUE)
  barcodes.file <-  grep('barcodes.tsv', files, fixed = TRUE, value = TRUE)

  # rename for ?Read10X
  file.rename(file.path(data_dir, c(mtx.file, genes.file, barcodes.file)),
              file.path(data_dir, c('matrix.mtx', 'genes.tsv', 'barcodes.tsv')))
}






#' Load mitochondrial and ribsomal gene names
#'
#' @inheritParams add_scseq_qc_metrics
#'
#' @return Named list with \code{rrna} and \code{mrna} character vectors.
#' @keywords internal
load_scseq_qcgenes <- function(species = 'Homo sapiens') {

  # load mito and ribo genes
  if (species == 'Homo sapiens') {
    rrna <- readLines(system.file('extdata', 'rrna.csv', package = 'dseqr', mustWork = TRUE))
    mrna <- readLines(system.file('extdata', 'mrna.csv', package = 'dseqr', mustWork = TRUE))

  } else if (species == 'Mus musculus') {
    rrna <- readLines(system.file('extdata', 'rrna_mouse.csv', package = 'dseqr', mustWork = TRUE))
    mrna <- readLines(system.file('extdata', 'mrna_mouse.csv', package = 'dseqr', mustWork = TRUE))

  } else {
    stop("Only 'Homo sapiens' and 'Mus musculus' supported")
  }

  return(list(rrna=rrna, mrna=mrna))
}

#' Add Mitochondrial and Ribosomal RNA gene Names to SingelCellExperiment
#'
#' @inheritParams add_scseq_qc_metrics
#'
#'
add_qc_genes <- function(sce, species) {
  # add qc genes as metadata
  qcgenes <- load_scseq_qcgenes(species)
  sce@metadata$mrna <- qcgenes$mrna
  sce@metadata$rrna <- qcgenes$rrna

  return(sce)
}


#' Utility wrapper to run normalization and log transformation
#'
#'
#' @param scseq \code{SingleCellExperiment} object
#'
#' @return Normalized and log transformed \code{scseq}.
#' @export
normalize_scseq <- function(scseq) {

  ncores <- min(parallel::detectCores(), 7)
  doParallel::registerDoParallel(ncores)

  set.seed(100)
  preclusters <- scran::quickCluster(scseq, BPPARAM = BiocParallel::DoparParam())
  scseq <- scran::computeSumFactors(scseq, cluster=preclusters)
  scseq <- scater::logNormCounts(scseq)

  return(scseq)
}

#' Get top highly variable genes for subsequent dimensionality reduction
#'
#' Runs after \code{preprocess_scseq}
#'
#' @param sce \code{SingleCellExperiment} object
#'
#' @export
add_hvgs <- function(sce, hvgs = NULL) {

  if (is.null(hvgs)) {
    dec <- scran::modelGeneVar(sce)
    hvg <- row.names(sce) %in% scran::getTopHVGs(dec, prop=0.1)
  } else {
    hvg <- row.names(sce) %in% hvgs
  }

  SummarizedExperiment::rowData(sce)$hvg <- hvg
  return(sce)
}


#' Run TSNE or UMAP
#'
#' @param sce \code{SingleCellExperiment}
#' @param dimred reducedDim to run TSNE or UMAP on
#'
#' @return \code{sce} with \code{'TSNE'} or \code{'UMAP'} \code{reducedDim}
#' @export
run_reduction <- function(sce, type = c('auto', 'TSNE', 'UMAP'), dimred = 'PCA') {

  if(type[1] == 'auto') {
    ncells <- ncol(sce)
    type <- ifelse(ncells > 5000, 'UMAP', 'TSNE')
  }

  set.seed(1100101001)
  if (type[1] == 'TSNE') {
    sce <- scater::runTSNE(sce, dimred = dimred, n_dimred = sce@metadata$npcs)

  } else if (type[1] == 'UMAP') {
    sce <- scater::runUMAP(sce, dimred = dimred, n_dimred = sce@metadata$npcs, min_dist=0.3)

  }
  colnames(SingleCellExperiment::reducedDim(sce, type)) <- paste0(type, 1:2)
  return(sce)
}





#' Get number of clusters different number of PCs
#'
#' Used to pick number of PCs to retain
#'
#' @param sce \code{SingleCellExperiement}
#' @inheritParams SingleCellExperiment::reducedDim
#'
#' @return result of \code{scran::getClusteredPCs}
#' @keywords internal
#'
get_npc_choices <- function(sce, type = 'PCA') {

  # walktrap very slow if too many cells
  cluster_fun <- ifelse(ncol(sce) > 10000,
                        igraph::cluster_louvain,
                        igraph::cluster_walktrap)

  FUN <- function(x, ...) {
    g <- scran::buildSNNGraph(x, ..., transposed = TRUE)
    cluster_fun(g)$membership
  }

  pcs <- SingleCellExperiment::reducedDim(sce, type = type)
  pcs <- pcs[, utils::head(seq_len(ncol(pcs)), 50)]

  choices <- scran::getClusteredPCs(pcs, FUN = FUN, by=2)
  names(choices$clusters) <- choices$n.pcs

  npcs <- S4Vectors::metadata(choices)$chosen
  cluster <- factor(choices$clusters[[as.character(npcs)]])


  return(list(npcs = npcs, cluster = cluster))
}


get_snn_graph <- function(sce, npcs = 30) {
  types <- SingleCellExperiment::reducedDimNames(sce)
  type <- ifelse('corrected' %in% types, 'corrected', 'PCA')

  # walktrap very slow if too many cells
  pcs <- SingleCellExperiment::reducedDim(sce, type = type)
  pcs <- pcs[, utils::head(seq_len(ncol(pcs)), npcs)]
  g <- scran::buildSNNGraph(pcs, transposed = TRUE)
  return(g)
}

get_clusters <- function(snn_graph, type = c('leiden', 'walktrap'), resolution = 1) {
  resolution <- as.numeric(resolution)

  if (type[1] == 'leiden') {
    cluster <- igraph::cluster_leiden(snn_graph, 'modularity', resolution_parameter = resolution)$membership

  } else if (type[1] == 'walktrap') {
    cluster <- igraph::cluster_walktrap(snn_graph)$membership
  }

  return(factor(cluster))
}

#' Cluster SingleCellExperiment
#'
#' @param sce \code{SingleCellExperiment}
#'
#' @return \code{sce} with column \code{cluster} in colData and \code{'npcs'} in metadata
#' @export
run_pca <- function(sce) {

  # run PCA on HVGs
  rdata <- SummarizedExperiment::rowData(sce)
  subset_row <- row.names(rdata[rdata$hvg, ])

  set.seed(100)
  sce <- scater::runPCA(sce, subset_row = subset_row)
  return(sce)
}


#' Calculate doublet score for each cell
#'
#' @param scseq \code{SingleCellExperiment} object with \code{hvg} column in \code{rowData}
#'   and \code{npcs} in \code{metadata} slot
#'
#' @return \code{scseq} with column \code{doublet_score} added to \code{colData}.
#' @export
#' @keywords internal
add_doublet_score <- function(scseq) {

  hvgs <- SingleCellExperiment::rowData(scseq)$hvg
  hvgs <- row.names(scseq)[hvgs]

  dbl.dens <- scDblFinder::computeDoubletDensity(scseq, subset.row=hvgs)
  scseq$doublet_score <- log10(dbl.dens+1)

  return(scseq)
}



#' Calculate QC metrics for SingleCellExperiment
#'
#' @param sce \code{SingleCellExperiment}
#' @param species Character indicating species. Either \code{'Homo sapiens'},
#' or \code{'Mus musculus'}.
#' @param for_qcplots Are the QC metrics being added for QC plots? Used by
#' \link{load_raw_scseq}.
#'
#' @return \code{sce} with qc metrics added by \code{\link[scater]{addPerCellQC}}
add_scseq_qc_metrics <- function(sce, species, for_qcplots = FALSE) {

  sce <- add_qc_genes(sce, species)
  sce <- scater::addPerCellQC(sce,
                              subsets=list(mito = which(row.names(sce) %in% sce@metadata$mrna),
                                           ribo = which(row.names(sce) %in% sce@metadata$rrna)))

  if (for_qcplots) sce <- add_scseq_qcplot_metrics(sce)

  return(sce)
}

#' Add QC metrics to SingleCellExperiment for plotting
#'
#' @param sce \code{SingleCellExperiment}
#'
#' @return \code{sce} with formated qc metrics.
add_scseq_qcplot_metrics <- function(sce) {

  sce$mito_percent <- sce$subsets_mito_percent
  sce$ribo_percent <- sce$subsets_ribo_percent
  sce$log10_sum <- log10(sce$sum)
  sce$log10_detected <- log10(sce$detected)

  return(sce)
}




#' Run pairwise wilcox tests between single cell clusters
#'
#' @param scseq \code{SingleCellExperiment} object.
#' @inheritParams scran::pairwiseWilcox
#'
#' @return List of \code{data.frame}s, one for each cluster.
#' @keywords internal
pairwise_wilcox <- function(scseq, groups = scseq$cluster, direction = 'up', block = NULL, restrict = NULL) {

  # dont get markers if no clusters
  if (!exist_clusters(scseq)) return(NULL)
  groups <- as.character(groups)

  # only upregulated as more useful for positive id of cell type
  wilcox_tests <- scran::pairwiseWilcox(SingleCellExperiment::logcounts(scseq),
                                        groups = groups,
                                        direction = direction,
                                        block = block,
                                        restrict = restrict)
  return(wilcox_tests)
}

#' Combine pairwise wilcox tests between single-cell clusters
#'
#' @param tests Result of \code{pairwise_wilcox}
#' @inheritParams scran::combineMarkers
#'
#' @return List of data.frames
#' @keywords internal
get_scseq_markers <- function(tests, pval.type = 'some', effect.field = 'AUC', keep = NULL) {
  if (is.null(keep)) keep <- rep(TRUE, nrow(tests$pairs))

  markers <- scran::combineMarkers(tests$statistics[keep],
                                   tests$pairs[keep, ],
                                   pval.type = pval.type,
                                   effect.field = effect.field)

  markers <- lapply(markers, as.data.frame)
  ord <- order(as.numeric(names(markers)))
  markers[ord]
}

get_presto_markers <- function(scseq) {
  markers <- presto::wilcoxauc(scseq, group_by = 'cluster', assay = 'logcounts')
  markers <- markers %>%
    dplyr::filter(logFC > 0) %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(pval -auc) %>%
    dplyr::group_split() %>%
    as.list()

  res <- list()
  for (df in markers) {
    df <- as.data.frame(df)
    row.names(df) <- df$feature
    df$feature <- NULL
    group <- as.character(df$group[1])
    res[[group]] <- df
  }

  ord <- order(as.numeric(names(res)))
  res <- res[ord]
  return(res)
}


#' Run fastMNN integration
#'
#' @param logcounts dgCMatrix of logcounts.
#' @inheritParams  batchelor::fastMNN
#' @param scseqs list of \code{SingleCellExperiment} objects to integrate.
#'
#' @return Integrated \code{SingleCellExperiment} with logcounts assay, corrected reducedDim, and batch annotation.
#' @keywords internal
run_fastmnn <- function(logcounts, subset.row, scseqs) {

  mnn.fun <- function(...) batchelor::fastMNN(
    ...,
    subset.row = subset.row,
    auto.merge = TRUE,
    correct.all = TRUE,
    cos.norm = FALSE,
    prop.k = 0.05)


  set.seed(1000101001)
  cor.out <- do.call(mnn.fun, scseqs) ; gc()

  # store merged (batch normalized) for DE
  SummarizedExperiment::assay(cor.out, 'logcounts') <- logcounts

  # remove things that dont use (bloats loom object)
  SummarizedExperiment::assay(cor.out, 'reconstructed') <- NULL
  SummarizedExperiment::rowData(cor.out)$rotation <- NULL

  return(cor.out)
}


#' Run harmony integration
#'
#' @inheritParams scater::calculatePCA
#' @inheritParams run_fastmnn
#' @param pairs
#'
#' @inherit run_fastmnn return
#' @keywords internal
run_harmony <- function(logcounts, subset_row, batch, pairs = NULL) {

  meta_data <- data.frame(batch = batch)
  vars_use <- 'batch'

  if (!is.null(pairs)) {
    vars_use <- c('batch', 'pairs')
    meta_data$pairs <- factor(pairs[batch, 'pair'])
  }


  set.seed(100)
  pcs <- scater::calculatePCA(logcounts, subset_row = subset_row)
  emb <- harmony::HarmonyMatrix(pcs, meta_data, vars_use, do_pca = FALSE, max.iter.harmony = 40)

  cor.out <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = logcounts),
    reducedDims = list(corrected = emb),
    colData =  S4Vectors::DataFrame(batch = batch))

  return(cor.out)
}



#' Integrate multiple scRNA-seq samples
#'
#' @param scseqs List of \code{SingleCellExperiment} objects
#' @param type One of \code{'harmony'} (default) or \code{'fastMNN'} specifying integration to use.
#' @param pairs data.frame with columns \code{'sample'} with sample
#'   names of \code{scseqs} and \code{'pair'} with integers indicating paired
#'   samples. If not \code{NULL} (default), then harmony includes pairings
#'   as a covariate.
#'
#' @return Integrated \code{SingleCellExperiment} object.
#' @keywords internal
integrate_scseqs <- function(scseqs, type = c('harmony', 'fastMNN'), pairs = NULL, hvgs = NULL, azimuth_ref = NULL) {

  # all common genes
  universe <- Reduce(intersect, lapply(scseqs, row.names))

  # variance modelling results
  decs <- lapply(scseqs, scran::modelGeneVar)

  # subset scseqs and decs
  decs <- lapply(decs, function(x) x[universe, ])
  scseqs <- lapply(scseqs, `[`, universe)

  # rescale each batch for depths
  scseqs <- do.call('multiBatchNorm', scseqs, envir = loadNamespace('batchelor'))

  # feature selection
  decs <- do.call('combineVar', decs, envir = loadNamespace('scran'))
  if (!is.null(hvgs)) hvgs <- hvgs %in% universe
  else hvgs <- decs$bio > 0

  # integration
  no_correct <- function(assay.type) function(...) batchelor::noCorrect(..., assay.type = assay.type)
  combined <- do.call(no_correct('logcounts'), scseqs)
  logcounts <- SummarizedExperiment::assay(combined, 'merged')
  gc()

  if (type[1] == 'harmony') {
    cor.out <- run_harmony(logcounts, hvgs, combined$batch, pairs = pairs)

  } else if (type[1] == 'fastMNN') {
    cor.out <- run_fastmnn(logcounts, hvgs, scseqs)

  } else if (type[1] == 'Azimuth') {
    azres <- run_azimuth(scseqs, azimuth_ref)
    cor.out <- transfer_azimuth(azres, combined)
    rm(azres); gc()

    # re-name merged to logcounts
    SummarizedExperiment::assayNames(cor.out) <- 'logcounts'

    # corrected is used to check if integrated
    SingleCellExperiment::reducedDim(cor.out, 'corrected') <- matrix(nrow = ncol(cor.out))
  }

  cor.out$orig.ident <- unlist(lapply(scseqs, `[[`, 'orig.ident'), use.names = FALSE)
  cor.out$orig.cluster <- unlist(lapply(scseqs, `[[`, 'cluster'), use.names = FALSE)
  cor.out$orig.resoln <- unlist(lapply(scseqs, `[[`, 'orig.resoln'), use.names = FALSE)

  rm(combined); gc()

  # get counts for pseudobulk
  counts <- do.call(no_correct('counts'), scseqs)
  SummarizedExperiment::assay(cor.out, 'counts') <- SummarizedExperiment::assay(counts, 'merged')
  rm(counts, scseqs); gc()

  return(cor.out)
}



#' Get genes that are ambient in at least one test and control sample
#'
#' @param scseqs List of \code{SingleCellExperiment} objects.
#'
#' @return List with test and control ambient genes
#' @keywords internal
get_ambience <- function(combined) {

  # ambient call for each gene in each cluster
  ambience <- SummarizedExperiment::rowData(combined)
  ambience <- ambience[, grepl('ambience$', colnames(ambience)), drop=FALSE]
  ambience <- as.matrix(ambience)
  return(ambience)

}

calc_cluster_ambience <- function(summed, ambience, clus) {
  counts <- SingleCellExperiment::counts(summed)
  counts <- counts[, summed$cluster == clus]
  amb <- DropletUtils::maximumAmbience(counts, ambience, mode = 'proportion')
  amb <- rowMeans(amb, na.rm=TRUE)
  amb <- amb[!is.na(amb)]
  amb <- names(amb)[amb > 0.1]
  return(amb)
}


#' Mark ambient outliers in combined dataset
#'
#' A gene is marked as an ambient outlier if it is an ambient outlier in at least one of the datasets.
#'
#' @param scseqs the original scseqs
#' @param combined the combined scseqs
#'
#' @return \code{combined} with \code{out_ambient} column added to \code{meta.features} slot of \code{SCT} assay.
#' @keywords internal
add_integrated_ambient <- function(combined, ambient) {

  genes <- row.names(combined)
  SummarizedExperiment::rowData(combined)$test_ambient <- genes %in% ambient$test
  SummarizedExperiment::rowData(combined)$ctrl_ambient <- genes %in% ambient$ctrl

  return(combined)
}


#' Test is there is at lest two clusters
#'
#' Used by \code{\link{get_scseq_markers}} to prevent getting markers if there are no clusters to compare
#'
#' @param scseq
#'
#' @return TRUE if more than one cluster exists
#' @keywords internal
exist_clusters <- function(scseq) {
  length(unique(scseq$cluster)) > 1
}

