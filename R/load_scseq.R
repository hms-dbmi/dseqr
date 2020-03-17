#' Process raw single cell fastq or cellranger files for app
#'
#' @param dataset_name Name of dataset
#' @param fastq_dir Directory with fastq or cellranger files
#' @param sc_dir Single cell directory for app. Will store results in \code{dataset_name} subdirectory
#' @param progress Optional shiny \code{Progress} object. Default will print progress.
#'
#' @return NULL
#' @export
#'
load_raw_scseq <- function(dataset_name, fastq_dir, sc_dir, indices_dir, progress = NULL, recount = FALSE, value = 0, founder = dataset_name, metrics = c('low_lib_size',
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
  if (is.cellranger) standardize_cellranger(fastq_dir)

  progress$set(message = "Quantifying files", value = value + 1)
  if (!is.cellranger) run_kallisto_scseq(indices_dir, fastq_dir, recount = recount)

  progress$set(message = "Loading", value + 2)
  type <- ifelse(is.cellranger, 'cellranger', 'kallisto')
  scseq <- create_scseq(fastq_dir, project = dataset_name, type = type)
  gc()

  progress$set(message = "Running QC", value = value + 3)
  scseq <- normalize_scseq(scseq)
  scseq <- add_hvgs(scseq)
  gc()

  scseq <- add_doublet_score(scseq)
  scseq <- add_scseq_qc_metrics(scseq, scseq@metadata$species, for_qcplots = TRUE)

  scseq <- run_scseq_qc(scseq, metrics)

  process_raw_scseq(scseq, dataset_name, sc_dir, founder, progress = progress, value = value + 3)
}

#' Process Count Data for App
#'
#' @param scseq \code{SingleCellExperiment}
#' @param dataset_name Name of dataset to save
#' @param sc_dir Directory to save dataset to
#' @param score_doublets Should doublet scores be computed? Intended for non-subsetted data within \code{\link{load_raw_scseq}}.
#'  Default is \code{FALSE}.
#' @param progress Shiny progress object. Default (\code{NULL}) prints to stdout.
#' @param value Initial value of progress.
#'
#' @return NULL
#' @export
process_raw_scseq <- function(scseq, dataset_name, sc_dir, founder = NULL, progress = NULL, value = 0) {

  if (is.null(progress)) {
    progress <- list(set = function(value, message = '', detail = '') {
      cat(value, message, detail, '...\n')
    })
  }

  progress$set(message = "Clustering", detail = '', value = value + 1)
  scseq <- normalize_scseq(scseq)
  scseq <- add_hvgs(scseq)
  scseq <- add_scseq_clusters(scseq)
  gc()


  progress$set(message = "Reducing dimensions", value = value + 2)
  scseq <- run_tsne(scseq)
  gc()

  progress$set(message = "Getting markers", value = value + 3)
  tests <- pairwise_wilcox(scseq)
  markers <- get_scseq_markers(tests)

  # top markers for SingleR
  top_markers <- scran::getTopMarkers(tests$statistics, tests$pairs)

  progress$set(message = "Saving", value = value + 4)
  anal <- list(scseq = scseq, markers = markers, tests = tests, annot = names(markers), top_markers = top_markers, founder = founder)

  save_scseq_data(anal, dataset_name, sc_dir)

  # don't save raw counts for loom (saved as non-sparse)
  SummarizedExperiment::assay(scseq, 'counts') <- NULL; gc()

  progress$set(message = "Saving loom", value = value + 5)
  save_scle(scseq, file.path(sc_dir, dataset_name))
  progress$set(value = value + 7)
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
#' @examples
#'
#' data_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'
#' load_scseq(data_dir)
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

  # get ambient expression profile/determine outlier genes
  # if pre-filtered cellranger, can't determine outliers/empty droplets
  ncount <- Matrix::colSums(counts)
  if (min(ncount) > 10) {
    pct_ambient <- rep(0, nrow(counts))
    out_ambient <- rep(FALSE, nrow(counts))
    keep_cells <- seq_len(ncol(counts))

  } else {
    pct_ambient <- get_pct_ambient(counts)
    out_ambient <- get_outliers(pct_ambient)
    keep_cells <- detect_cells(counts, species = species)
  }

  # add ambient metadata for genes
  rowData <- S4Vectors::DataFrame(pct_ambient, out_ambient)
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
#' Loads scle.loom if exists otherwise scseq.rds
#'
#' @param dataset_dir Directory with scle.loom or scseq.rds
#'
#' @return \code{SingleCellExperiment} or \code{SingleCellLoomExperiment}
#' @export
#' @keywords internal
load_scseq <- function(dataset_dir) {
  scle_path <- file.path(dataset_dir, 'scle.loom')
  scseq_path <- file.path(dataset_dir, 'scseq.rds')

  # load loom if available (faster and less memory)
  scseq <- tryCatch(LoomExperiment::import(scle_path, type = 'SingleCellLoomExperiment'),
                    error = function(e) {
                      unlink(scle_path)
                      return(NULL)
                    })

  if (!is.null(scseq)) {
    # fixes for SCLE
    colnames(SingleCellExperiment::reducedDim(scseq, 'TSNE')) <- c('TSNE1', 'TSNE2')
    scseq$cluster <- factor(as.numeric(scseq$cluster))

    is.integrated <- !is.null(scseq$orig.ident) && all(scseq$orig.ident %in% c('test', 'ctrl'))
    if (is.integrated) {
      scseq$orig.ident <- factor(scseq$orig.ident, levels = c('test', 'ctrl'))
      scseq$orig.cluster <- factor(as.numeric(scseq$orig.cluster))
    }

  } else {
    scseq <- readRDS(scseq_path)
  }

  return(scseq)
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

#' Determine ambient percent for each gene
#'
#' Looks at droplets with counts less than or equal to 10 calculates the total percent for each gene
#'
#' @param counts \code{dgTMatrix} of counts. Rows are genes, columns are droplets.
#'
#' @return Named numeric vector of percentages for each gene.
#' @export
#' @keywords internal
get_pct_ambient <- function(counts) {

  # get drops with less than 10 counts
  ncount <- Matrix::colSums(counts)
  ambient <- counts[, ncount <= 10]

  # percentage of counts per gene
  nambient <- Matrix::rowSums(ambient)
  pct_ambient <- nambient / sum(nambient) * 100
  return(pct_ambient)
}

#' Flag outliers
#'
#'
#' @param x Named numeric vector
#'
#' @return Boolean vector with \code{length(x)} indicating if values of \code{x} are outliers (TRUE) or not (FALSE).
#' @export
#' @keywords internal
get_outliers <- function(x) {
  outliers <- graphics::boxplot(x, plot = FALSE)$out
  is.outlier <- names(x) %in% names(outliers)
  return(is.outlier)
}



#' Read kallisto/bustools market matrix and annotations
#'
#'
#' @inheritParams load_scseq
#'
#' @return sparse dgTMatrix with barcodes in columns and genes in rows.
#' @export
load_kallisto_counts <- function(data_dir) {

  # read sparse matrix
  counts <- Matrix::readMM(file.path(data_dir, 'genecount', 'genes.mtx'))
  counts <- Matrix::t(counts)
  counts <- as(counts, 'dgCMatrix')

  # read annotations
  row.names(counts) <- readLines(file.path(data_dir, 'genecount', 'genes.genes.txt'))
  colnames(counts) <- readLines(file.path(data_dir, 'genecount', 'genes.barcodes.txt'))

  # remove non-expressed genes
  counts <- counts[Matrix::rowSums(counts) > 0, ]

  return(counts)
}

#' Load cell ranger counts
#'
#' Mainly to avoid having to download massive datasets that have already been quantified.
#'))))
#' @inheritParams load_scseq
#' @importFrom magrittr "%>%"
#'
#' @return dgCMatrix
#' @export
load_cellranger_counts <- function(data_dir) {
  # read the data in using ENSG features
  h5file <- list.files(data_dir, '.h5$', full.names = TRUE)
  if (length(h5file)) {
    counts <- Read10X_h5(h5file, use.names = FALSE)
  } else {
    counts <- Read10X(data_dir, gene.column = 1)
  }

  if (class(counts) == 'list') counts <- counts$`Gene Expression`
  return(counts)
}

process_cellranger_counts <- function(counts, species) {

  if (species == 'Homo sapiens') tx2gene <- tx2gene
  else if (species == 'Mus musculus') tx2gene <- tx2gene_mouse
  else stop('Species not supported')

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
#' @export
#'
#' @examples
#' \dontrun{
#' # For output from CellRanger < 3.0
#' data_dir <- 'path/to/data/directory'
#' list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
#' expression_matrix <- Read10X(data.dir = data_dir)
#' seurat_object = CreateSeuratObject(counts = expression_matrix)
#'
#' # For output from CellRanger >= 3.0 with multiple data types
#' data_dir <- 'path/to/data/directory'
#' list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
#' data <- Read10X(data.dir = data_dir)
#' seurat_object = CreateSeuratObject(counts = data$`Gene Expression`)
#' seurat_object[['Protein']] = CreateAssayObject(counts = data$`Antibody Capture`)
#' }
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
    feature.names <- read.delim(
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
    list_of_data[[j]] <- as(object = list_of_data[[j]], Class = "dgCMatrix")
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
#' @export
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
      giveCsparse = FALSE
    )
    if (unique.features) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = 'dgCMatrix')
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
#' @export
check_is_cellranger <- function(data_dir) {

  # cellranger file names
  files <- list.files(data_dir)
  mtx.file <- grep('.mtx', files, fixed = TRUE, value = TRUE)
  genes.file <- grep('features.tsv|genes.tsv', files, value = TRUE)
  barcodes.file <-  grep('barcodes.tsv', files, fixed = TRUE, value = TRUE)

  if (length(mtx.file) & length(genes.file) & length(barcodes.file)) return(TRUE)
  return(FALSE)
}

#' Rename CellRanger files for loading by Read10X
#'
#' @param data_dir Path to folder with CellRanger files.
#'
#' @return NULL
#' @export
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
#' This are the genes used by alevin for whitelisting.
#'
#' @return Named list with \code{rrna} and \code{mrna} character vectors.
#' @export
load_scseq_qcgenes <- function(species = 'Homo sapiens') {

  # load mito and ribo genes
  if (species == 'Homo sapiens') {
    rrna <- readLines(system.file('extdata', 'rrna.csv', package = 'drugseqr', mustWork = TRUE))
    mrna <- readLines(system.file('extdata', 'mrna.csv', package = 'drugseqr', mustWork = TRUE))

  } else if (species == 'Mus musculus') {
    rrna <- readLines(system.file('extdata', 'rrna_mouse.csv', package = 'drugseqr', mustWork = TRUE))
    mrna <- readLines(system.file('extdata', 'mrna_mouse.csv', package = 'drugseqr', mustWork = TRUE))

  } else {
    stop("Only 'Homo sapiens' and 'Mus musculus' supported")
  }

  return(list(rrna=rrna, mrna=mrna))
}

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
#' @return
#' @export
add_hvgs <- function(sce) {

  dec <- scran::modelGeneVar(sce)
  hvg <- row.names(sce) %in% scran::getTopHVGs(dec, prop=0.1)
  SummarizedExperiment::rowData(sce)$hvg <- hvg

  return(sce)
}


#' Run TSNE
#'
#' @param sce \code{SingleCellExperiment}
#' @param dimred reducedDim to run TSNE on
#'
#' @return \code{sce} with \code{'TSNE'} \code{reducedDim}
#' @export
run_tsne <- function(sce, dimred = 'PCA') {
  set.seed(1100101001)
  sce <- scater::runTSNE(sce, dimred = dimred, n_dimred = sce@metadata$npcs)
  colnames(SingleCellExperiment::reducedDim(sce, 'TSNE')) <- c('TSNE1', 'TSNE2')
  return(sce)
}


#' Get number of clusters different number of PCs
#'
#' Used to pick number of PCs to retain
#'
#' @param sce \code{SingleCellExperiement}
#'
#' @return result of \code{scran::getClusteredPCs}
#' @export
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
  pcs <- pcs[, head(seq_len(ncol(pcs)), 50)]

  choices <- scran::getClusteredPCs(pcs, FUN = FUN)
  names(choices$clusters) <- choices$n.pcs

  npcs <- S4Vectors::metadata(choices)$chosen
  cluster <- factor(choices$clusters[[as.character(npcs)]])


  return(list(npcs = npcs, cluster = cluster))
}

#' Cluster SingleCellExperiment
#'
#' @param sce \code{SingleCellExperiment}
#'
#' @return \code{sce} with column \code{cluster} in colData and \code{'npcs'} in metadata
#' @export
add_scseq_clusters <- function(sce) {

  # run PCA on HVGs
  rdata <- SummarizedExperiment::rowData(sce)
  subset_row <- row.names(rdata[rdata$hvg, ])

  set.seed(100)
  sce <- scater::runPCA(sce, subset_row = subset_row)

  # pick number of PCs
  choices <- get_npc_choices(sce)
  sce@metadata$npcs <- choices$npcs

  # add clusters
  sce$cluster <- choices$cluster
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

  dbl.dens <- scran::doubletCells(scseq, subset.row=hvgs)
  scseq$doublet_score <- log10(dbl.dens+1)

  return(scseq)
}



#' Calculate QC metrics for SingleCellExperiment
#'
#' @param sce \code{SingleCellExperiment}
#'
#' @return \code{sce} with qc metrics added by \code{\link[scater]{addPerCellQC}}
#' @export
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
#' @export
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
#'
#' @return List of \code{data.frame}s, one for each cluster.
#' @export
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
#' @param pval.type
#'
#' @return List of data.frames
#' @export
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


#' Run fastMNN integration
#'
#' @param logcounts dgCMatrix of logcounts.
#' @inheritParams  batchelor::fastMNN
#' @param scseqs list of \code{SingleCellExperiment} objects to integrate.
#'
#' @return Integrated \code{SingleCellExperiment} with logcounts assay, corrected reducedDim, and batch annotation.
#' @export
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

#' Run liger integration
#'
#' @inheritParams run_fastmnn
#' @param batch Character vector indicating the original datasets for each cell in scseqs.
#'  Must be in same order as \code{scseqs}.
#'
#' @inherit run_fastmnn return
#' @export
#' @keywords internal.
run_liger <- function(logcounts, scseqs, batch) {
  countl <- list()
  for (i in seq_along(scseqs)) {
    dataset_name <- names(scseqs)[i]
    counts <- SingleCellExperiment::counts(scseqs[[i]])
    colnames(counts) <- paste(colnames(counts), dataset_name, sep = '::')
    countl[[dataset_name]] <- counts
  }

  countl <- liger::createLiger(countl); gc()
  countl <- liger::normalize(countl); gc()
  countl <- liger::selectGenes(countl); gc()
  countl <- liger::scaleNotCenter(countl); gc()
  countl <- liger::optimizeALS(countl, k=20); gc()
  countl <- liger::quantile_norm(countl); gc()

  corrected <- countl@H.norm
  row.names(corrected) <- gsub('::.+?$', '', row.names(corrected))

  cor.out <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = logcounts),
    reducedDims = list(corrected = corrected),
    colData =  S4Vectors::DataFrame(batch = batch)
  )

  return(cor.out)


}


#' Run harmony integration
#'
#' @inheritParams scater::calculatePCA
#' @inheritParams run_liger
#' @param pairs
#'
#' @inherit run_fastmnn return
#' @export
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
#' @param type One of \code{'harmony'} (default), \code{'liger'} or \code{'fastMNN'} specifying integration to use.
#'
#' @return Integrated \code{SingleCellExperiment} object.
#' @export
integrate_scseqs <- function(scseqs, type = c('harmony', 'liger', 'fastMNN'), pairs = NULL) {

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
  hvgs <- decs$bio > 0

  # integration
  no_correct <- function(assay.type) function(...) batchelor::noCorrect(..., assay.type = assay.type)
  combined <- do.call(no_correct('logcounts'), scseqs)
  logcounts <- SummarizedExperiment::assay(combined, 'merged')
  gc()

  if (type[1] == 'harmony') {
    cor.out <- run_harmony(logcounts, hvgs, combined$batch, pairs = pairs)

  } else if (type[1] == 'liger') {
    cor.out <- run_liger(logcounts, scseqs, combined$batch)

  } else if (type[1] == 'fastMNN') {
    cor.out <- run_fastmnn(logcounts, hvgs, scseqs)
  }

  cor.out$orig.ident <- unlist(lapply(scseqs, `[[`, 'orig.ident'), use.names = FALSE)
  cor.out$orig.cluster <- unlist(lapply(scseqs, `[[`, 'cluster'), use.names = FALSE)


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
#' @export
#' @keywords internal
get_integrated_ambient <- function(scseqs) {

  # datasets that are test samples
  is.test <- sapply(scseqs, function(x) levels(x$orig.ident) == 'test')

  # genes that are ambient in at least one test sample
  ambient.test <- lapply(scseqs[is.test], function(x) SingleCellExperiment::rowData(x))
  ambient.test <- lapply(ambient.test, function(x) row.names(x)[x$out_ambient])
  ambient.test <- unique(unlist(ambient.test))

  # genes that are ambient in at least one ctrl sample
  ambient.ctrl <- lapply(scseqs[!is.test], function(x) SingleCellExperiment::rowData(x))
  ambient.ctrl <- lapply(ambient.ctrl, function(x) row.names(x)[x$out_ambient])
  ambient.ctrl <- unique(unlist(ambient.ctrl))

  return(list(test = ambient.test, ctrl = ambient.ctrl))
}

#' Mark ambient outliers in combined dataset
#'
#' A gene is marked as an ambient outlier if it is an ambient outlier in at least one of the datasets.
#'
#' @param scseqs the original scseqs
#' @param combined the combined scseqs
#'
#' @return \code{combined} with \code{out_ambient} column added to \code{meta.features} slot of \code{SCT} assay.
#' @export
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
#' @export
exist_clusters <- function(scseq) {
  length(unique(scseq$cluster)) > 1
}

