#' Load raw RNA-Seq data into an ExpressionSet.
#'
#' @param data_dir Directory with raw and quantified RNA-Seq files.
#' @param species Character vector indicating species. Genus and species should be space seperated, not underscore. Default is \code{Homo sapiens}.
#' @param release EnsemblDB release. Should be same as used in \code{\link{build_index}}.
#' @param overwrite If FALSE (default) and a saved \code{ExpressionSet} exists, will load from disk.
#' @param dgel If TRUE, will also save the \code{DGEList} object. Used for testing purposes.
#'
#' @return \code{\link[Biobase]{ExpressionSet}} with attributes/accessors:
#' \itemize{
#'   \item \code{sampleNames} from names of raw RNA-Seq files (excluding .fastq.gz suffix).
#'   \item \code{annotation} Character vector of annotation package used.
#'   \item \code{exprs} Length scaled counts generated from abundances for use in
#'     \code{\link[limma]{voom}} (see \code{vignette("tximport", package = "tximport")}).
#'   \item \code{phenoData} from \code{pdata_path} annotation file with added columns:
#'     \itemize{
#'       \item \code{quants_dir} directory within '\code{data_dir}/quants' containing salmon quantification results.
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
#' eset <- load_seq(data_dir, overwrite = TRUE)
#'
load_seq <- function(data_dir, species = 'Homo sapiens', release = '94', overwrite = FALSE, dgel = FALSE) {

  # check if already have
  eset_path  <- file.path(data_dir, 'eset.rds')
  if (!overwrite & file.exists(eset_path))
    return(readRDS(eset_path))

  pdata_path <- file.path(data_dir, 'pdata.rds')
  if (!file.exists(pdata_path)) stop("No 'pdata.rds' file found in data_dir. Did you run_salmon first?")
  pdata <- readRDS(pdata_path)

  if (species != 'Homo sapiens') stop('only implemented for Homo sapiens')


  # transcript to gene map
  tx2gene <- get_tx2gene(species, release)

  # import quant.sf files and filter low counts
  quants <- import_quants(data_dir, tx2gene)

  # add library normalization
  pdata <- add_norms(quants, pdata)

  # construct eset
  annot <- get_ensdb_package(species, release)
  fdata <- setup_fdata(tx2gene)
  eset <- construct_eset(quants, fdata, pdata, annot)

  if (dgel) {
    dgel_path <- gsub('eset', 'dgel', eset_path)
    saveRDS(quants, dgel_path)
  }

  # save eset and return
  saveRDS(eset, eset_path)
  return(eset)
}

#' Construct expression set
#'
#' @param quants \code{DGEList} with RNA-seq counts.
#' @param fdata \code{data.table} returned from \code{\link{setup_fdata}}.
#' @param pdata \code{data.frame} of sample annotation data processed by \code{\link{match_pdata}} followed by \code{\link{add_norms}}.
#' @param annot Character vector with ensembldb package name. e.g. \code{'EnsDb.Hsapiens.v94'}. Returned from \code{\link{get_ensdb_package}}.
#'
#' @return \code{ExpressionSet} with input arguments are in corresponding atributes.
#' @keywords internal
#' @export
#'
construct_eset <- function(quants, fdata, pdata, annot) {
  # remove duplicate rows of counts
  rn <- row.names(quants$counts)
  mat <- unique(data.table::data.table(quants$counts, rn, key = 'rn'))

  # merge exprs and fdata
  dt <- merge(fdata, mat, by.y = 'rn', by.x = 'SYMBOL', all.y = TRUE, sort = FALSE)
  dt <- as.data.frame(dt)
  row.names(dt) <- make.unique(dt[[1]])

  # seperate fdata and exprs and transfer to eset
  row.names(pdata) <- pdata$quants_dir

  # Replicate column no longer needed (one row kept per replicate)
  pdata$Replicate <- NULL

  eset <- Biobase::ExpressionSet(as.matrix(dt[, row.names(pdata), drop=FALSE]),
                                 phenoData=Biobase::AnnotatedDataFrame(pdata),
                                 featureData=Biobase::AnnotatedDataFrame(dt[, colnames(fdata), drop=FALSE]),
                                 annotation=annot)

  return(eset)
}

#' Setup feature annotation data
#'
#' @param tx2gene \code{data.frame} mapping transcripts to gene names. Returned from from \code{\link{get_tx2gene}}.
#'
#' @return \code{data.table} with columns \code{SYMBOL} and \code{ENTREZID_HS} corresponding to
#'   HGNC symbols and human entrez ids respectively.
#' @keywords internal
#' @export
#'
setup_fdata <- function(tx2gene) {

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

#' Import salmon quant.sf files
#'
#' @param data_dir Directory with a folder named 'quants' that contains salmon quantification folders for each sample.
#' @inheritParams setup_fdata
#'
#' @return \code{DGEList} with length scaled counts. Lowly expressed genes are filtered.
#' @keywords internal
#' @export
#'
import_quants <- function(data_dir, tx2gene) {

  # don't ignoreTxVersion if dots in tx2gene
  ignore <- TRUE
  if (any(grepl('[.]', tx2gene$tx_id))) ignore <- FALSE

  # import quants using tximport
  # using limma::voom for differential expression (see tximport vignette)
  qdirs   <- list.files(file.path(data_dir, 'quants'))
  quants_paths <- file.path(data_dir, 'quants', qdirs, 'quant.sf')

  # use folders as names (used as sample names)
  names(quants_paths) <- qdirs

  txi <- tximport::tximport(quants_paths, tx2gene = tx2gene, type = "salmon",
                            ignoreTxVersion = ignore, countsFromAbundance = "lengthScaledTPM", importer=utils::read.delim)

  quants <- edgeR::DGEList(txi$counts)

  # filtering low counts (as in tximport vignette)

  keep <- edgeR::filterByExpr(quants)
  quants <- quants[keep, ]
  if (!nrow(quants)) stop("No genes with reads after filtering")

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
#'
#' @return \code{data.frame} with columns \code{tx_id}, \code{gene_name}, and \code{entrezid}
#' @export
#'
#' @keywords internal
#'
get_tx2gene <- function(species, release) {
  # load EnsDb package
  ensdb_package <- get_ensdb_package(species, release)
  if (!require(ensdb_package, character.only = TRUE)) {
    build_ensdb(species, release)
    require(ensdb_package, character.only = TRUE)
  }

  # map from transcripts to genes
  tx2gene <- ensembldb::transcripts(get(ensdb_package), columns=c("tx_id", "gene_name", "entrezid"),
                                    return.type='data.frame')
  tx2gene[tx2gene == ""] <- NA
  tx2gene <- tx2gene[!is.na(tx2gene$gene_name), ]
  return(tx2gene)
}

#' Add edgeR normalization factors to pdata.
#'
#' @inheritParams construct_eset
#' @param pdata \code{data.frame} with sample annotations.
#'
#' @return \code{pdata} with \code{lib.size} and \code{norm.factors} columns added.
#' @export
#' @keywords internal
#'
add_norms <- function(quants, pdata) {
  idxs <- match(colnames(quants), pdata$quants_dir)
  if (any(is.na(idxs))) stop("add_norms failed to match quant directory names and file names")
  pdata <- pdata[idxs,, drop=FALSE]
  pdata <- cbind(pdata, quants$samples[, c('lib.size', 'norm.factors')])
  return(pdata)
}

#' Match and add file names for pdata rows.
#'
#' Attempts to find exact match between file names (excluding fastq.gz suffix) and first pdata column.
#' If an exact match is not possible, \code{grep} is used to search for a unique substring.
#'
#' @param pdata \code{data.table} with sample annotations.
#' @param fastqs Character vector of fastq.gz files to find matches for.
#'
#' @return \code{pdata} with \code{file} column identifying file for each sample (row).
#' @export
#' @keywords internal
#'
match_pdata <- function(pdata, fastqs) {
  ids <- pdata[[1]]

  if (length(ids) != length(fastqs)) stop('Must be one row in sample annotation per fastq.gz file.')

  # check if every id has a perfect match in fastq.gz prefixes
  idxs <- match(ids, gsub('.fastq.gz$', '', fastqs))

  # otherwise try unique grep
  if (sum(is.na(idxs)) != 0) {

    # remove _ and - as word seperation (for grep below)
    spaced_ids <- gsub('_|-', ' ', ids)
    spaced_qdirs <- gsub('_|-', ' ', fastqs)

    tomatch <- spaced_ids[is.na(idxs)]
    gmatch  <- sapply(tomatch, function(id) {

      res <- grep(paste0('\\b', id, '\\b'), spaced_qdirs, ignore.case = TRUE)

      if (length(res) > 1)
        stop(shQuote(id), ' is a substring of more than one fastq.gz file. Fix in first column of sample annotation file.')

      if (!length(res))
        stop(shQuote(id), ' is not a substring of a fastq.gz file. Fix in first column of sample annotation file.')

      return(res)
    })

    # fill in grep result where NA in idxs
    idxs[is.na(idxs)] <- gmatch[is.na(idxs)]
  }

  # append file names to pdata
  pdata <- tibble::add_column(pdata, 'File Name' = fastqs[idxs], .after = 1)
  return(pdata)
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
