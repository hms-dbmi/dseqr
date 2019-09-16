#' Load bulk RNA-Seq data into an ExpressionSet.
#'
#' @param data_dir Directory with raw and quantified RNA-Seq files.
#' @param type Either \code{'salmon'} or \code{'kallisto'}. The package used for quantification. Must be on the PATH.
#' @param species Character vector indicating species. Genus and species should be space seperated, not underscore. Default is \code{Homo sapiens}.
#' @param release EnsemblDB release. Should be same as used in \code{\link{build_ensdb_index}}.
#' @param load_saved If TRUE (default) and a saved \code{ExpressionSet} exists, will load from disk.
#' @param save_eset If TRUE (default) and either \code{load_saved} is \code{FALSE} or a saved \code{ExpressionSet} does not exist,
#'   then an ExpressionSet will be saved to disk. Will overwrite if already exists.
#' @param save_dgel If TRUE, will save the \code{DGEList} object. Used for testing purposes. Default is \code{FALSE}.
#' @param filter if \code{TRUE} (default), low count genes are filtered using \code{\link[edgeR]{filterByExpr}}.
#'
#' @return \code{\link[Biobase]{ExpressionSet}} with attributes/accessors:
#' \itemize{
#'   \item \code{sampleNames} from names of raw RNA-Seq files (excluding .fastq.gz suffix).
#'   \item \code{annotation} Character vector of annotation package used.
#'   \item \code{exprs} Length scaled counts generated from abundances for use in
#'     \code{\link[limma]{voom}} (see \code{vignette("tximport", package = "tximport")}).
#'   \item \code{phenoData} from \code{pdata_path} annotation file with added columns:
#'     \itemize{
#'       \item \code{lib.size} library size from \code{\link[edgeR]{calcNormFactors}}.
#'       \item \code{norm.factors} library normalization factors from \code{\link[edgeR]{calcNormFactors}}.
#'     }
#' }
#'
#' @export
#'
#' @examples
#'
#' data_dir <- 'data-raw/example-data'
#' eset <- load_seq(data_dir, load_saved = FALSE, save_eset = FALSE)
#'
load_seq <- function(data_dir, type = 'kallisto', species = 'Homo sapiens', release = '94', load_saved = TRUE, save_eset = TRUE, save_dgel = FALSE, filter = TRUE) {

  # check if already have
  eset_path  <- file.path(data_dir, 'eset.rds')
  if (load_saved & file.exists(eset_path))
    return(readRDS(eset_path))

  # import quants and filter low counts
  quants <- import_quants(data_dir, filter, type = type)

  # construct eset
  annot <- get_ensdb_package(species, release)
  fdata <- setup_fdata()
  eset <- construct_eset(quants, fdata, annot)

  if (save_dgel) {
    # used for testing purposes
    dgel_path <- gsub('eset', 'dgel', eset_path)
    saveRDS(quants, dgel_path)
  }

  # save eset and return
  if (save_eset) saveRDS(eset, eset_path)
  return(eset)
}

#' Construct expression set
#'
#' @param quants \code{DGEList} with RNA-seq counts.
#' @param fdata \code{data.table} returned from \code{\link{setup_fdata}}.
#' @param annot Character vector with ensembldb package name. e.g. \code{'EnsDb.Hsapiens.v94'}. Returned from \code{\link{get_ensdb_package}}.
#'
#' @return \code{ExpressionSet} with input arguments are in corresponding atributes.
#' @keywords internal
#' @export
#'
construct_eset <- function(quants, fdata, annot) {
  # remove duplicate rows of counts
  rn <- row.names(quants$counts)
  mat <- unique(data.table::data.table(quants$counts, rn, key = 'rn'))

  # merge exprs and fdata
  dt <- merge(fdata, mat, by.y = 'rn', by.x = 'SYMBOL', all.y = TRUE, sort = FALSE)
  dt <- as.data.frame(dt)
  row.names(dt) <- make.unique(dt[[1]])

  # remove edgeR assigned group column
  pdata <- quants$samples
  pdata$group <- NULL

  # seperate fdata and exprs and transfer to eset
  eset <- Biobase::ExpressionSet(as.matrix(dt[, row.names(pdata), drop=FALSE]),
                                 phenoData=Biobase::AnnotatedDataFrame(pdata),
                                 featureData=Biobase::AnnotatedDataFrame(dt[, colnames(fdata), drop=FALSE]),
                                 annotation=annot)

  return(eset)
}

#' Setup feature annotation data
#'
#' @return \code{data.table} with columns \code{SYMBOL} and \code{ENTREZID_HS} corresponding to
#'   HGNC symbols and human entrez ids respectively.
#' @keywords internal
#' @export
#'
setup_fdata <- function() {

  # unlist entrezids
  fdata <- data.table::data.table(tx2gene)
  fdata <- fdata[, list(ENTREZID_HS = as.character(unlist(entrezid)))
                 , by = c('tx_id', 'gene_name')]

  fdata[, tx_id := NULL]
  fdata <- unique(fdata)

  # setup so that will work with crossmeta
  fdata <- fdata[, .(SYMBOL = gene_name, ENTREZID_HS)]
  return(fdata)
}

#' Import kallisto quants
#'
#' @inheritParams setup_fdata
#' @inheritParams load_seq
#'
#' @return \code{DGEList} with length scaled counts. Lowly expressed genes are filtered.
#' @keywords internal
#' @export
#'
import_quants <- function(data_dir, filter, type) {

  # don't ignoreTxVersion if dots in tx2gene
  ignore <- TRUE
  if (any(grepl('[.]', tx2gene$tx_id))) ignore <- FALSE

  # import quants using tximport
  # using limma::voom for differential expression (see tximport vignette)
  pkg_version <- get_pkg_version(type)
  qdirs   <- list.files(file.path(data_dir, paste0(type, pkg_version, 'quants', sep = '_')))


  quants_paths <- ifelse(type == 'kallisto',
                         file.path(data_dir, 'quants', qdirs, 'abundance.h5'),
                         file.path(data_dir, 'quants', qdirs, 'quant.sf'))

  # use folders as names (used as sample names)
  names(quants_paths) <- qdirs

  txi <- tximport::tximport(quants_paths, tx2gene = tx2gene, type = type,
                            ignoreTxVersion = ignore, countsFromAbundance = "lengthScaledTPM")

  quants <- edgeR::DGEList(txi$counts)

  # filtering low counts (as in tximport vignette)
  if (filter) {
    keep <- edgeR::filterByExpr(quants)
    quants <- quants[keep, ]
    if (!nrow(quants)) stop("No genes with reads after filtering")
  }

  quants <- edgeR::calcNormFactors(quants)
  return(quants)
}

#' Get ensembldb package name
#'
#' @inheritParams load_seq
#'
#' @keywords internal
#' @return Character vector with ensembldb package name. e.g. \code{'EnsDb.Hsapiens.v94'}.
#' @export
#'
get_ensdb_package <- function(species, release) {
  ensdb_species    <- strsplit(species, ' ')[[1]]
  ensdb_species[1] <- toupper(substr(ensdb_species[1], 1, 1))

  ensdb_package <- paste('EnsDb', paste0(ensdb_species, collapse = ''), paste0('v', release), sep='.')
  return(ensdb_package)
}

#' Get transcript to gene map.
#'
#' @inheritParams load_seq
#' @param columns Character vector of columns from ensdb package to return or \code{'list'} to print available options.
#'
#' @return \code{data.frame} with columns \code{tx_id}, \code{gene_name}, and \code{entrezid}
#' @export
#'
#' @keywords internal
#'
get_tx2gene <- function(species = 'Homo sapiens', release = '94', columns = c("tx_id", "gene_name", "entrezid", "gene_id", "seq_name", "description")) {
  is.installed('ensembldb', level = 'error')

  # load EnsDb package
  ensdb_package <- get_ensdb_package(species, release)
  if (!require(ensdb_package, character.only = TRUE)) {
    build_ensdb(species, release)
    require(ensdb_package, character.only = TRUE)
  }

  if (columns[1] == 'list') {
    print(ensembldb::listColumns(get(ensdb_package)))
    return(NULL)
  }

  # map from transcripts to genes
  tx2gene <- ensembldb::transcripts(get(ensdb_package), columns=columns, return.type='data.frame')
  tx2gene[tx2gene == ""] <- NA
  tx2gene$description <- gsub(' \\[Source.+?\\]', '', tx2gene$description)
  tx2gene <- tx2gene[!is.na(tx2gene$gene_name), ]
  return(tx2gene)
}


#' Build ensembldb annotation package.
#'
#' @inheritParams load_seq
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' # build ensembldb annotation package for human
#' build_ensdb()
#'
build_ensdb <- function(species = 'Homo sapiens', release = '94') {
  # packages in suggests for faster install time
  suggests <- c('ensembldb', 'AnnotationHub', 'AnnotationDbi')
  is.installed(suggests, level = 'error')

  # store ensembl databases in built package
  ensdb_dir <- 'EnsDb'
  unlink('EnsDb', recursive = TRUE)
  dir.create(ensdb_dir)

  # format is genus_species in multiple other functions but not here
  species <- gsub('_', ' ', species)

  # generate new ensembl database from specified release
  ah <- AnnotationHub::AnnotationHub()
  ahDb <- AnnotationHub::query(ah, pattern = c(species, "EnsDb", release))

  if (!length(ahDb)) stop('Specified ensemble species/release not found in AnnotationHub.')

  ahEdb <- ahDb[[1]]

  ensembldb::makeEnsembldbPackage(AnnotationDbi::dbfile(ensembldb::dbconn(ahEdb)),
                                  '0.0.1', 'Alex Pickering <alexvpickering@gmail.com>',
                                  'Alex Pickering',
                                  ensdb_dir)

  # install new ensemble database
  ensdb_name <- list.files(ensdb_dir)
  ensdb_path <- file.path(ensdb_dir, ensdb_name)
  install.packages(ensdb_path, repos = NULL)

  # remove source files
  unlink(ensdb_dir, recursive = TRUE)
}

is.installed <- function(packages, level = c('quiet', 'error', 'message')) {
  level <- level[1]

  for (package in packages) {

    if (!requireNamespace(package, quietly = TRUE)) {
      msg <- paste0("Package \"", package, "\" needed for this function to work. Please install it.")
      if (level == 'error') stop(msg, call. = FALSE)
      if (level == 'message') message(msg)
      return(FALSE)
    }
  }

  return(TRUE)
}
