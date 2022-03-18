#' Import raw single cell fastq or cellranger files for app
#'
#' @param dataset_name Name of dataset
#' @param uploaded_data_dir Directory with fastq or cellranger files
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
import_scseq <- function(dataset_name,
                         uploaded_data_dir,
                         sc_dir,
                         tx2gene_dir,
                         indices_dir = NULL,
                         progress = NULL,
                         recount = FALSE,
                         value = 0,
                         founder = dataset_name,
                         npcs = 30,
                         cluster_alg = 'leiden',
                         resoln = 1,
                         species = NULL,
                         ref_name = NULL,
                         metrics = c('low_lib_size',
                                     'low_n_features',
                                     'high_subsets_mito_percent',
                                     'low_subsets_ribo_percent',
                                     'high_doublet_score')) {

  if (!is.null(metrics) && metrics[1] == 'none') metrics <- NULL
  if (is.null(progress)) {
    progress <- list(set = function(value, message = '', detail = '') {
      cat(value, message, detail, '...\n')
    })
  }
  # check if saved R object
  robject <- find_robject(uploaded_data_dir)

  if (length(robject)) {

    import_robject(
      dataset_name,
      uploaded_data_dir,
      sc_dir,
      species,
      tx2gene_dir,
      metrics,
      progress = progress,
      value = value)

    return(NULL)
  }

  # check if cellranger and if so standardize file names
  is.cellranger <- check_is_cellranger(uploaded_data_dir)

  if (!is.cellranger) {
    progress$set(message = "running pseudoalignment", value = value + 1)
    run_kallisto_scseq(indices_dir, uploaded_data_dir, recount = recount)
  }

  progress$set(message = "loading", value + 2)
  type <- ifelse(is.cellranger, 'cellranger', 'kallisto')
  scseq <- create_scseq(uploaded_data_dir, tx2gene_dir, dataset_name, type)
  gc()

  progress$set(message = "running QC", value = value + 3)
  scseq <- add_doublet_score(scseq)
  scseq <- add_scseq_qc_metrics(scseq, for_qcplots = TRUE)

  scseq <- run_scseq_qc(scseq, metrics)

  process_raw_scseq(scseq,
                    dataset_name,
                    sc_dir,
                    cluster_alg = cluster_alg,
                    npcs = npcs,
                    resoln = resoln,
                    ref_name = ref_name,
                    founder = founder,
                    progress = progress,
                    value = value + 3,
                    tx2gene_dir = tx2gene_dir)
}

check_is_robject <- function(uploaded_data_dir) {

}

find_robject <- function(uploaded_data_dir, load = FALSE) {

  files <- list.files(uploaded_data_dir, full.names = TRUE)
  rds.file <- grep('[.]rds$', files, value = TRUE)
  qs.file <- grep('[.]qs$', files, value = TRUE)

  if (length(rds.file)) {
    fpath <- rds.file[1]
    load_func <- readRDS

  } else if (length(qs.file)) {
    fpath <- qs.file[1]
    load_func <- qs::qread

  } else {
    return(NULL)
  }

  if (load) return(load_func(fpath))
  else return(fpath)
}



import_robject <- function(dataset_name, uploaded_data_dir, sc_dir, species, tx2gene_dir, metrics, progress = NULL, value = 0) {

  if (is.null(progress)) {
    progress <- list(set = function(value, message = '', detail = '') {
      cat(value, message, detail, '...\n')
    })
  }

  ## load the R object
  progress$set(message = "loading R object", value = value + 1)
  scseq <- find_robject(uploaded_data_dir, load = TRUE)

  if (is(scseq, 'Seurat')) {
    scseq <- seurat_to_sce(scseq, dataset_name)
    gc()
  }

  samples <- unique(scseq$batch)
  multisample <- length(samples) > 1

  ## process the R object
  progress$set(message = "processing R object", value = value + 2)

  scseqs <- list()
  need <- get_need_multisample(scseq)

  if (!multisample | need$process_samples)
    scseqs <- process_robject_samples(scseq, tx2gene_dir, species, metrics)

  if (multisample) {
    scseq <- process_robject_multisample(scseq, scseqs, need)
  } else {
    scseq <- scseqs[[1]]
  }

  red.names <- SingleCellExperiment::reducedDimNames(scseq)

  # add species
  scseq$project <- dataset_name
  scseq@metadata$species <- species

  # get snn graph and initial resolution
  ref_name <- scseq@metadata$ref_name
  is_ref <- !is.null(ref_name)

  snn_possible <-
    (!multisample & 'PCA' %in% red.names) |
    ( multisample & 'corrected' %in% red.names)

  if (is_ref) {
    # will grab other resolutions from scseq@colData
    resoln <- scseq@metadata$resoln
    snn_graph <- NULL

  } else if (snn_possible) {
    # allow resolution changes if supplied corrected
    resoln <- scseq@metadata$resoln
    if (is.null(resoln)) resoln <- 1
    snn_graph <- get_snn_graph(scseq)

  } else {
    # prevent resolution changes when provided UMAP/TSNE and clusters
    # but not corrected for multi-sample or PCA for uni-sample
    snn_graph <- NULL
    resoln <- 'provided.clusters'
  }

  # provided clusters were removed if not provided:
  # corrected for multi-sample or PCA for uni-sample
  use_provided_clusters <- !is.null(scseq$cluster)
  if (!use_provided_clusters) {
    resoln <- 1
    scseq$cluster <- get_clusters(snn_graph, resolution = resoln)
  }

  annot <- levels(scseq$cluster)
  levels(scseq$cluster) <- seq_along(levels(scseq$cluster))

  progress$set(message = "saving", value = value + 3)

  # store what is stable with resolution change
  scseq_data <- list(scseq = scseq,
                     snn_graph = snn_graph,
                     species = species,
                     founder = dataset_name,
                     resoln = resoln,

                     # metadata to reproduce exports
                     ref_name = ref_name,
                     meta = scseq@metadata$meta,
                     prev_groups = scseq@metadata$prev_groups)

  # remove NULLs
  scseq_data <- scseq_data[!sapply(scseq_data, is.null)]
  save_scseq_data(scseq_data, dataset_name, sc_dir, add_integrated = multisample)

  # save multi-sample specific data
  if (multisample) {
    args <- list()

    if (is_ref) {
      args$ref_name <- ref_name
      args$integration_type <- get_ref_type(ref_name)
      save_ref_clusters(scseq@colData, dataset_name, sc_dir)

    } else {
      args$integration_type <- 'harmony'
    }

    save_scseq_args(args = args, dataset_name, sc_dir)
  }

  # run what depends on resolution
  run_post_cluster(scseq, dataset_name, sc_dir, resoln = resoln)


  # save things in resolution sub-directory
  dataset_subname <- file.path(dataset_name, get_resoln_dir(resoln))
  scseq_subdata <- list(annot = annot)

  # flag to indicate resolution of provided clusters
  if (use_provided_clusters) scseq_subdata$provided_clusters <- TRUE

  save_scseq_data(scseq_subdata, dataset_subname, sc_dir, overwrite = FALSE)
}

get_need_multisample <- function(scseq) {
  red.names <- SingleCellExperiment::reducedDimNames(scseq)

  # get reduction if missing
  # will remove supplied clusters
  need_reduction <- !any(c('TSNE', 'UMAP') %in% red.names)

  # need corrected to get reduction and clusters
  # if have clusters/reduction but not corrected will disable resolution changes
  no_clusters <- is.null(scseq$cluster)
  need_corrected <- !'corrected' %in% red.names & (need_reduction | no_clusters)

  # normalize to get logcounts if missing
  need_normalize <- !'logcounts' %in% names(scseq@assays)

  need_add_qc <- is.null(scseq$mito_percent)

  need_integrate <- need_corrected | need_normalize

  return(list(
    reduction = need_reduction,
    integrate = need_integrate,
    add_qc = need_add_qc,
    process_samples = need_integrate | need_add_qc
  ))
}


process_robject_multisample <- function(scseq, scseqs, need) {

  # use harmony if had to get logcounts per sample or need corrected to get UMAP/TSNE
  if (need$integrate) {
    project <- scseq$project[1]
    message('getting corrected ...')
    message('supplied clusters will be replaced')

    scseq <- integrate_scseqs(scseqs, type='harmony')
    scseq$project <- project
  }

  # transfer QC metrics from individual datasets
  if (need$integrate | need$add_qc)
    scseq <- add_combined_metrics(scseq, scseqs)

  if (need$reduction) {
    message('getting UMAP/TSNE ...')
    scseq <- run_reduction(scseq, dimred = 'corrected')
  }

  return(scseq)
}


process_robject_samples <- function(scseq, tx2gene_dir, species, metrics) {
  samples <- unique(scseq$batch)
  unisample <- length(samples) == 1

  rdata <- SummarizedExperiment::rowData(scseq)
  red.names <- SingleCellExperiment::reducedDimNames(scseq)

  ## checks required for unisample ---

  # get doublet scores if unisample and missing
  need_doublets <- unisample & is.null(scseq$doublet_score)

  # run QC if unisample and have QC metrics
  need_run_qc <- unisample & !is.null(metrics)

  # get reduction if unisample and missing
  need_reduction <- unisample & !any(c('TSNE', 'UMAP') %in% red.names)

  # get PCA if unisample, missing, and need UMAP or no clusters
  no_clusters <- is.null(scseq$cluster)
  need_run_pca <- unisample & !'PCA' %in% red.names & (need_reduction | no_clusters)
  if (need_run_pca) scseq$cluster <- NULL

  # get HVGs and BIO if need run PCA or unisample and BIO missing
  need_hvgs <- need_run_pca | (unisample & !'bio' %in% names(rdata))

  ## checks required for uni/multi-sample ---

  # normalize to get logcounts if missing
  need_normalize <- !'logcounts' %in% names(scseq@assays)

  # add QC metrics if they are missing
  need_add_qc <- is.null(scseq$mito_percent)

  # add mito/ribo genes if adding QC metrics and missing
  if (need_add_qc && is.null(scseq@metadata$mrna)) {
    message('adding mito/ribo genes ...')
    tx2gene <- load_tx2gene(species, tx2gene_dir)
    qcgenes <- load_scseq_qcgenes(species, tx2gene)

    scseq@metadata$mrna <- qcgenes$mrna
    scseq@metadata$rrna <- qcgenes$rrna
  }

  # get list of scseqs
  samples <- unique(scseq$batch)
  scseqs <- list()

  for (sample in samples) {
    scseqi <- scseq[, scseq$batch %in% sample]

    if (need_doublets) scseqi <- add_doublet_score(scseqi)
    if (need_add_qc) scseqi <- add_scseq_qc_metrics(scseqi, for_qcplots = TRUE)
    if (need_run_qc) scseqi <- run_scseq_qc(scseqi, metrics)
    if (need_normalize) scseqi <- normalize_scseq(scseqi)

    # reduction related
    if (need_hvgs) scseqi <- add_hvgs(scseqi)
    if (need_run_pca) scseqi <- run_pca(scseqi)
    if (need_reduction) scseqi <- run_reduction(scseqi)

    scseqs[[sample]] <- scseqi
  }

  return(scseqs)
}

#' Convenience utility to run import_scseq in background by callr::r_bg
#'
#' Calls import_scseq
#'
#' @keywords internal
#' @noRd
run_import_scseq <- function(opts, uploaded_data_dir, sc_dir, tx2gene_dir, indices_dir, species, ref_name = NULL) {

  for (opt in opts) {
    import_scseq(opt$dataset_name,
                 uploaded_data_dir,
                 sc_dir,
                 tx2gene_dir,
                 indices_dir,
                 metrics = opt$metrics,
                 founder = opt$founder,
                 species = species,
                 ref_name = ref_name)

  }

  return(TRUE)
}


#' Process Count Data for App
#'
#' Performs the following:
#' - normalization
#' - adds HVGs
#' - dimensionality reduction
#' - clustering
#'
#' Then calls run_post_cluster
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
                              tx2gene_dir,
                              cluster_alg = 'leiden',
                              npcs = 30,
                              resoln = 1,
                              hvgs = NULL,
                              ref_name = NULL,
                              founder = NULL,
                              progress = NULL,
                              value = 0) {

  if (is.null(progress)) {
    progress <- list(set = function(value, message = '', detail = '') {
      cat(value, message, detail, '...\n')
    })
  }

  progress$set(message = "normalizing", value = value)
  scseq@metadata$npcs <- npcs
  scseq <- normalize_scseq(scseq)

  # add HVGs
  scseq <- add_hvgs(scseq, hvgs)
  species <- scseq@metadata$species

  is_ref <- !is.null(ref_name)
  if (is_ref) ref_type <- get_ref_type(ref_name)

  if (!is_ref) {
    progress$set(message = "reducing", value = value + 1)
    scseq <- run_pca(scseq)
    scseq <- run_reduction(scseq)
    gc()

    progress$set(message = "clustering", detail = '', value = value + 2)
    snn_graph <- get_snn_graph(scseq, npcs)
    scseq$cluster <- get_clusters(snn_graph, cluster_alg, resoln)
    gc()

    # save items that are independent of resolution
    anal <- list(scseq = scseq, snn_graph = snn_graph, founder = founder, resoln = resoln, species = species)
    save_scseq_data(anal, dataset_name, sc_dir)

  } else if (ref_type == 'Azimuth') {
    progress$set(message = "running Azimuth", detail = '', value = value + 2)

    resoln <- get_ref_resoln(ref_name)
    azres <- run_azimuth(list(one = scseq), ref_name, species, tx2gene_dir)
    scseq <- transfer_azimuth(azres, scseq, resoln)
    rm(azres); gc()

    anal <- list(scseq = scseq, founder = founder, resoln = resoln, ref_name = ref_name, species = species)
    save_scseq_data(anal, dataset_name, sc_dir)
    save_ref_clusters(scseq@colData, dataset_name, sc_dir)

  } else if (ref_type == 'symphony') {
    progress$set(message = "running symphony", detail = '', value = value + 2)

    resoln <- get_ref_resoln(ref_name)
    counts <- SingleCellExperiment::counts(scseq)
    logcounts <- SingleCellExperiment::logcounts(scseq)

    symres <- run_symphony(counts, logcounts, ref_name, scseq$project, species, tx2gene_dir)
    scseq <- transfer_symphony(symres, scseq)
    rm(symres); gc()

    anal <- list(scseq = scseq, founder = founder, resoln = resoln, ref_name = ref_name, species = species)
    save_scseq_data(anal, dataset_name, sc_dir)
    save_ref_clusters(scseq@colData, dataset_name, sc_dir)

  }

  # run what depends on resolution
  run_post_cluster(scseq,
                   dataset_name,
                   sc_dir,
                   resoln,
                   progress,
                   value + 3,
                   reset_annot = !is_ref)

  progress$set(value = value + 7)
}

transfer_symphony <- function(symres, scseq) {

  scseq$cluster <- symres$cluster
  scseq$predicted.celltype <- symres$predicted.celltype
  scseq$predicted.celltype.score <- symres$predicted.celltype.score


  SingleCellExperiment::reducedDim(scseq, 'UMAP') <-
    SingleCellExperiment::reducedDim(symres, 'UMAP')

  return(scseq)
}

transfer_azimuth <- function(azres, scseq, resoln) {

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

  scseq$mapping.score <- get_meta(azres, 'mapping.score')

  cols <- colnames(azres[[1]]@meta.data)
  azi_cols <- get_ref_cols(cols)
  for (col in azi_cols) scseq[[col]] <- get_meta(azres, col)

  clus <- factor(scseq[[resoln]])
  scseq$cluster <- factor(as.numeric(clus))
  return(scseq)
}

get_resoln_dir <- function(resoln) {
  is.num <- is.numstring(resoln)
  ifelse(is.num,
         paste0('snn', resoln),
         resoln)
}

is.numstring <- function(x) !is.na(suppressWarnings(as.numeric(x)))

get_ref_cols <- function(cols, type = c('both', 'score', 'cluster')) {

  score_cols <- grep('^predicted[.].+?[.]score$', cols, value = TRUE)
  if (type[1] == 'score') return(score_cols)

  if (length(score_cols)) clust_cols <- gsub('^(predicted[.].+?)[.]score$', '\\1', score_cols)
  else clust_cols <- grep('^predicted[.].+?$', cols, value = TRUE)

  if (type[1] == 'cluster') return(clust_cols)
  return(c(score_cols, clust_cols))
}


run_azimuth <- function(scseqs, azimuth_ref, species, tx2gene_dir) {

  # prevents error about it being 500MB
  old <- options(future.globals.maxSize = 30*1024*1024^2)
  on.exit(options(old), add = TRUE)

  ref_species <- refs$species[refs$name == azimuth_ref]
  reference <- dseqr.data::load_data(paste0(azimuth_ref, '.qs'))

  pat <- '^celltype|^annotation|^class$|^cluster$|^subclass$|^cross_species_cluster$'
  refnames <- grep(pat, colnames(reference$map@meta.data), value = TRUE)

  refdata <- lapply(refnames, function(x) {
    reference$map[[x, drop = TRUE]]
  })

  names(refdata) <- refnames

  if ('ADT' %in% names(reference$map@assays)) {
    refdata[["impADT"]] <- Seurat::GetAssayData(
      object = reference$map[['ADT']],
      slot = 'data'
    )
  }

  queries <- list()
  for (ds in names(scseqs)) {

    scseq <- scseqs[[ds]]
    counts <- SingleCellExperiment::counts(scseq)

    # convert to reference species gene names
    if (species != ref_species)
      counts <- convert_species(counts, tx2gene_dir, species, ref_species)

    # error if query and ref have identical cell names
    orig.cells <- colnames(counts)
    colnames(counts) <- paste0('query-', orig.cells)

    query <- Seurat::CreateSeuratObject(counts = counts,
                                        min.cells = 1, min.features = 1)

    is.sct <- 'SCTModel.list' %in% slotNames(reference$map[['refAssay']])

    if (is.sct) {
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
    } else {
      # log normalize
      query <- Seurat::NormalizeData(object = query)
      query <- Seurat::FindVariableFeatures(query)
      query <- SeuratObject::RenameAssays(query, 'RNA' = 'refAssay')
    }


    # error if more than ncells
    ncells <- ncol(scseq)
    k.map <- min(100, round(ncells/4))
    k.score <- min(30, ncells-1)

    # Find anchors between query and reference
    anchors <- Seurat::FindTransferAnchors(
      reference = reference$map,
      query = query,
      k.filter = NA,
      reference.neighbors = "refdr.annoy.neighbors",
      reference.assay = "refAssay",
      query.assay = "refAssay",
      reference.reduction = "refDR",
      normalization.method = ifelse(is.sct, 'SCT', 'LogNormalize'),
      features = intersect(rownames(reference$map), Seurat::VariableFeatures(query)),
      dims = 1:50,
      n.trees = 20,
      k.score = k.score,
      mapping.score.k = k.map,
      verbose = TRUE
    )

    # Transfer cell type labels and impute protein expression
    #
    # Transferred labels are in metadata columns named "predicted.*"
    # The maximum prediction score is in a metadata column named "predicted.*.score"
    # The prediction scores for each class are in an assay named "prediction.score.*"
    # The imputed assay is named "impADT" if computed
    query <- tryCatch(Seurat::TransferData(
      reference = reference$map,
      query = query,
      dims = 1:50,
      anchorset = anchors,
      refdata = refdata,
      n.trees = 20,
      k.weight = round(k.map/2),
      store.weights = TRUE,
      verbose = FALSE
    ),
    error = function(e) return(NULL))

    if (is.null(query)) next()

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
    meta.data <- reference$map[[]]
    if ("ori.index" %in% colnames(meta.data)) {
      query <- NNTransform(
        object = query,
        meta.data = reference$map[[]]
      )
    }

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
      metadata = Seurat::MappingScore(anchors = anchors, ksmooth = k.map, kanchors = round(k.map/2)),
      col.name = "mapping.score"
    )

    # restore cell names
    query <- Seurat::RenameCells(query, new.names = gsub('^query-', '', colnames(query)))
    queries[[ds]] <- query
  }

  return(queries)
}

save_ref_clusters <- function(meta, dataset_name, sc_dir) {

  # cluster columns
  cols <- get_ref_cols(colnames(meta), 'cluster')

  # save annotation and clusters for each
  for (col in cols) {
    dataset_subname <- file.path(dataset_name, col)
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
create_scseq <- function(data_dir, tx2gene_dir, project, type = c('kallisto', 'cellranger')) {

  # load counts
  #TODO suport mouse for kallisto
  if (type[1] == 'kallisto') {
    data_dir <- file.path(data_dir, 'bus_output')
    counts <- load_kallisto_counts(data_dir)
    species <- 'Homo sapiens'
    tx2gene <- load_tx2gene(species, tx2gene_dir)

  } else if (type[1] == 'cellranger') {
    counts <- load_cellranger_counts(data_dir)
    species <- get_species(counts)

    # use supplied names as alternative
    alt_genes <- load_cellranger_genes(data_dir)

    tx2gene <- load_tx2gene(species, tx2gene_dir)
    counts <- process_cellranger_counts(counts, tx2gene, alt_genes)
  }

  # get ambience expression profile/determine outlier genes
  # if pre-filtered cellranger, can't determine outliers/empty droplets
  ncount <- Matrix::colSums(counts)
  qcgenes <- load_scseq_qcgenes(species, tx2gene)

  # will be ~98% of drops if filtered
  pct.lt10 <- sum(ncount < 10)/length(ncount)

  if (pct.lt10 < 0.20) {
    keep_cells <- seq_len(ncol(counts))
    rowData <- NULL

  } else {
    ambience <- DropletUtils::estimateAmbience(counts, good.turing = FALSE, round = FALSE)
    keep_cells <- detect_cells(counts, qcgenes)
    rowData <- S4Vectors::DataFrame(ambience)
  }

  project <- rep(project, length(keep_cells))

  # add ambience metadata for genes
  colData <- S4Vectors::DataFrame(project = project, batch = project)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts[, keep_cells]),
    colData = colData,
    rowData = rowData
  )

  sce@metadata$species <- species
  sce@metadata$mrna <- qcgenes$mrna
  sce@metadata$rrna <- qcgenes$rrna
  return(sce)
}

load_tx2gene <- function(species, tx2gene_dir) {
  # simplify species name
  ensdb_species <- strsplit(species, " ")[[1]]
  ensdb_species[1] <- tolower(substr(ensdb_species[1], 1, 1))
  ensdb_species <- paste(ensdb_species, collapse = '')

  # load if previously saved otherwise save
  fname <- paste0(ensdb_species, '_tx2gene.qs')
  fpath <- file.path(tx2gene_dir,  fname)

  if (file.exists(fpath)) {
    tx2gene <- qs::qread(fpath)
  } else {
    tx2gene <- dseqr.data::load_tx2gene(species, release = NULL, with_hgnc = TRUE)
    qs::qsave(tx2gene, fpath)
  }

  return(tx2gene)
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

  sample_group <- meta$group
  is.test <- sample_group == groups[1]
  is.ctrl <- sample_group == groups[2]
  sample_group[is.test] <- 'test'
  sample_group[is.ctrl] <- 'ctrl'
  other <- setdiff(sample_group, c('test', 'ctrl'))
  names(sample_group) <- row.names(meta)
  cell_group <- unname(sample_group[scseq$batch])

  scseq$orig.ident <- factor(cell_group, levels = c('test', 'ctrl', other))
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
#' @export
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
  tmp <- file.path(dirname(fpath), 'tmp.qs')

  if (file.copy(fpath, tmp)) {
    res <- try(qs::qread(tmp), silent = TRUE)
    if (class(res) != "try-error") {
      unlink(fpath)
      file.move(tmp, fpath)
    }
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
  h5file <- list.files(data_dir, '[.]h5$|[.]hdf5$', full.names = TRUE)
  if (length(h5file)) {
    counts <- Seurat::Read10X_h5(h5file, use.names = FALSE)
  } else {
    counts <- Seurat::Read10X(data_dir, gene.column = 1)
  }

  # remove species prefix from ensembl ids
  ensids <- gsub('^[^_]+[_]+(ENS.+?)$', '\\1', row.names(counts))

  enids <- gsub('[.]\\d+', '', ensids)
  row.names(counts) <- make.unique(enids)

  if (methods::is(counts, 'list')) counts <- counts$`Gene Expression`
  return(counts)
}


load_cellranger_genes <- function(data_dir) {
  # read the data in using ENSG features
  h5file <- list.files(data_dir, '[.]h5$|[.]hdf5$', full.names = TRUE)
  gene.file <- list.files(data_dir, 'features.tsv|genes.tsv', full.names = TRUE)[1]

  if (length(h5file)) {
    infile <- hdf5r::H5File$new(h5file, 'r')
    genomes <- names(infile)
    slot <- ifelse('matrix' %in% genomes, 'features/name', 'gene_names')
    genes <- infile[[file.path(genomes[1], slot)]][]

  } else {
    genes <- read.table(gene.file)[[2]]
  }

  return(genes)
}

process_cellranger_counts <- function(counts, tx2gene, alt_genes) {
  gene_id <- gene_name <- NULL

  # name genes by tx2gene so that best match with cmap/l1000 data
  map <- tx2gene %>%
    dplyr::filter(gene_id %in% row.names(counts)) %>%
    dplyr::select(gene_id, gene_name) %>%
    dplyr::distinct()

  idx <- match(map$gene_id, row.names(counts))
  row.names(counts)[idx] <- map$gene_name

  # use supplied gene names where no match
  no.match <- setdiff(seq_len(nrow(counts)), idx)
  row.names(counts)[no.match] <- alt_genes[no.match]

  # remove non-expressed genes
  # counts <- counts[Matrix::rowSums(counts) > 0, ]

  # sum counts in rows with same gene
  counts <- Matrix.utils::aggregate.Matrix(counts, row.names(counts), fun = 'sum')
  return(counts)
}

#' Get species for human or mouse depending on gene ids
#'
#' @param x dgCMatrix or data.frame where rownames are gene ids
#'
#' @return either 'Homo sapiens' or 'Mus musculus'
#' @export
#' @keywords internal
#'
get_species <- function(x) {

  # remove species prefix
  ensid <- gsub('^[^_]+[_]+(ENS.+?)$', '\\1', row.names(x))

  ensid <- grep('^ENS[A-Z]*G[0-9.]+$', ensid, value = TRUE)[1]
  ensid <- gsub('^(ENS[A-Z]*)G[0-9.]+$', '\\1', ensid)
  if (is.na(ensid)) stop('Need Ensembl IDs')

  return(ensmap[ensid, 'species'])

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
  h5.file <- grepl('[.]h5$|[.]hdf5$', files)

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
load_scseq_qcgenes <- function(species = 'Homo sapiens', tx2gene = NULL) {

  # load mito and ribo genes
  if (species == 'Homo sapiens') {
    rrna <- readLines(system.file('extdata', 'rrna.csv', package = 'dseqr', mustWork = TRUE))
    mrna <- readLines(system.file('extdata', 'mrna.csv', package = 'dseqr', mustWork = TRUE))

  } else if (species == 'Mus musculus') {
    rrna <- readLines(system.file('extdata', 'rrna_mouse.csv', package = 'dseqr', mustWork = TRUE))
    mrna <- readLines(system.file('extdata', 'mrna_mouse.csv', package = 'dseqr', mustWork = TRUE))

  } else {
    rrna <- NULL
    mrna <- tx2gene$gene_name[tx2gene$seq_name == 'MT']
  }

  return(list(rrna=rrna, mrna=mrna))
}

#' Utility wrapper to run normalization and log transformation
#'
#'
#' @param scseq \code{SingleCellExperiment} object
#'
#' @return Normalized and log transformed \code{scseq}.
#' @export
normalize_scseq <- function(scseq) {

  set.seed(100)
  preclusters <- scran::quickCluster(scseq)
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

  # always add bio even if supplied hvgs
  dec <- scran::modelGeneVar(sce)
  SummarizedExperiment::rowData(sce)$bio <- dec$bio

  if (is.null(hvgs)) {
    hvgs <- scran::getTopHVGs(dec, prop=0.1)
  }

  SummarizedExperiment::rowData(sce)$hvg <- row.names(sce) %in% hvgs
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

  # make sure > 200 counts to dodge errors
  counts <- SingleCellExperiment::counts(scseq)
  pass.check <- Matrix::colSums(counts) > 200
  counts <- counts[, pass.check]

  # is error prone
  res <- tryCatch(
    scDblFinder::scDblFinder(counts)@colData,

    error = function(e) {
      message('scDblFinder failed. Proceeding anyway.')

      res <- data.frame(row.names = colnames(counts))
      res$scDblFinder.class <- factor('singlet', levels = c('singlet', 'doublet'))
      res$scDblFinder.score <- 0
      return(res)
    })

  # assume droplets with <= 200 counts are singlets
  scseq$doublet_score <- 0
  scseq$doublet_class <- factor('singlet', levels = c('singlet', 'doublet'))

  scseq$doublet_score[pass.check] <- res$scDblFinder.score
  scseq$doublet_class[pass.check] <- res$scDblFinder.class

  return(scseq)
}



#' Calculate QC metrics for SingleCellExperiment
#'
#' @param sce \code{SingleCellExperiment}
#' @param for_qcplots Are the QC metrics being added for QC plots? Used by
#' \link{load_raw_scseq}.
#'
#' @return \code{sce} with qc metrics added by \code{\link[scater]{addPerCellQC}}
add_scseq_qc_metrics <- function(sce, for_qcplots = FALSE) {

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


#' Get cluster markers using presto
#'
#' @param scseq SingleCellExperiment object
#'
#' @return list of data.frames, one for each cluster
#' @export
#' @importFrom presto wilcoxauc
#' @importFrom magrittr "%>%"
#'
#' @keywords internal
get_presto_markers <- function(scseq) {
  markers <- presto::wilcoxauc(scseq, group_by = 'cluster', assay = 'logcounts', verbose = TRUE)
  markers <- markers %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(-auc) %>%
    dplyr::group_split() %>%
    as.list()

  res <- list()
  for (df in markers) {
    df <- as.data.frame(df)
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
integrate_scseqs <- function(scseqs, species, tx2gene_dir, type = c('harmony', 'fastMNN', 'Azimuth', 'symphony'), pairs = NULL, hvgs = NULL, ref_name = NULL) {

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

  no_correct <- function(assay.type) function(...) batchelor::noCorrect(..., assay.type = assay.type)
  combined <- do.call(no_correct('logcounts'), scseqs)
  logcounts <- SummarizedExperiment::assay(combined, 'merged')
  gc()

  # integration
  if (type[1] == 'harmony') {
    cor.out <- run_harmony(logcounts, hvgs, combined$batch, pairs = pairs)

  } else if (type[1] == 'fastMNN') {
    cor.out <- run_fastmnn(logcounts, hvgs, scseqs)

  } else if (type[1] == 'Azimuth') {

    azres <- run_azimuth(scseqs, ref_name, species, tx2gene_dir)

    # exclude unsuccessful samples
    combined <- combined[, combined$batch %in% names(azres)]
    scseqs <- scseqs[names(azres)]

    resoln <- get_ref_resoln(ref_name)
    cor.out <- transfer_azimuth(azres, combined, resoln)
    rm(azres); gc()

    # re-name merged to logcounts
    SummarizedExperiment::assayNames(cor.out) <- 'logcounts'

    # corrected is used to check if integrated
    SingleCellExperiment::reducedDim(cor.out, 'corrected') <- matrix(nrow = ncol(cor.out))

  } else if (type[1] == 'symphony') {
    all.counts <- do.call(no_correct('counts'), scseqs)
    all.counts <- SummarizedExperiment::assay(all.counts, 'merged')
    cor.out <- run_symphony(all.counts, logcounts, ref_name, combined$batch, species, tx2gene_dir)
  }

  # bio is used for sorting markers when no cluster selected
  SummarizedExperiment::rowData(cor.out)$bio <- decs$bio

  rm(combined); gc()

  # get counts for pseudobulk
  if (!exists('all.counts')) {
    all.counts <- do.call(no_correct('counts'), scseqs)
    dimnames(all.counts) <- dimnames(cor.out)
    SummarizedExperiment::assay(cor.out, 'counts') <- SummarizedExperiment::assay(all.counts, 'merged')
  }

  rm(all.counts, scseqs); gc()
  return(cor.out)
}

run_symphony <- function(counts, logcounts, ref_name, batch, species, tx2gene_dir) {

  meta_data <- data.frame(batch = batch)
  vars <- 'batch'

  # load reference
  reference <- dseqr.data::load_data(paste0(ref_name, '_reference.qs'))

  # convert to reference species if need to
  counts_use <- counts
  ref_species <- refs[refs$name == ref_name, 'species']
  if (species != ref_species) {
    counts_use <- convert_species(
      counts,
      tx2gene_dir,
      species,
      other_species = ref_species)
  }

  # get uwot model
  uwot_name <- paste0(ref_name, '_uwot_model')
  dl_dir <- system.file(package = "dseqr.data", mustWork = TRUE)
  uwot_path <- file.path(dl_dir, 'extdata', uwot_name)
  if (!file.exists(uwot_path)) dseqr.data::dl_data(uwot_name)

  # indicate model path for mapQuery
  reference$save_uwot_path <- uwot_path

  # map query
  query <- mapQuery(counts_use,
                              meta_data,
                              reference,
                              vars = vars,
                              do_normalize = TRUE)

  query <- knn_predict(query, reference, reference$meta_data$cell_type)
  meta <- query$meta_data

  cor.out <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts, logcounts = logcounts),
    reducedDims = list(corrected = t(query$Z), UMAP = query$umap),
    colData =  S4Vectors::DataFrame(
      batch = batch,
      cluster = meta$predicted.celltype,
      predicted.celltype = meta$predicted.celltype,
      predicted.celltype.score = meta$predicted.celltype.score))



  return(cor.out)
}


knn_predict <- function(query_obj, ref_obj, train_labels, k = 5, save_as = "predicted.celltype",
                        confidence = TRUE, seed = 0) {

  na.label <- is.na(train_labels)

  ref_obj$Z_corr <- ref_obj$Z_corr[, !na.label]
  train_labels <- train_labels[!na.label]

  set.seed(seed)
  if (confidence) {
    knn_pred = FNN::knn(t(ref_obj$Z_corr), t(query_obj$Z),
                        train_labels, k = k, prob = TRUE)
    knn_prob = attributes(knn_pred)$prob
    query_obj$meta_data[save_as] = knn_pred
    query_obj$meta_data[paste0(save_as, ".score")] = knn_prob
  }
  else {
    knn_pred = FNN::knn(t(ref_obj$Z_corr), t(query_obj$Z),
                        train_labels, k = k, prob = FALSE)
    query_obj$meta_data[save_as] = knn_pred
  }
  return(query_obj)
}


get_ref_resoln <- function(ref_name) {
  switch(ref_name,
         'human_pbmc' = 'predicted.celltype.l2',
         'human_lung' = 'predicted.annotation.l1',
         'human_lung_v2' = 'predicted.ann_level_4',
         'human_bonemarrow' = 'predicted.celltype.l2',
         'human_differentiated_tcell' = 'predicted.celltype.cytokines',
         'human_motorcortex' = 'predicted.subclass',
         'mouse_motorcortex' = 'predicted.subclass',
         'predicted.celltype'
  )
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
  if (ncol(counts) != ncol(ambience)) return(NULL)
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

