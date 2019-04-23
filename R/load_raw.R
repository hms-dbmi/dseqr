#' Build ensembldb annotation package.
#'
#' @param species Character vector indicating species. Genus and species should be space seperated, not underscore. Default is \code{Homo sapiens}.
#' @param release EnsemblDB release. Should be same as used in \code{\link{build_index}}.
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


#' Load raw RNA-Seq data into an ExpressionSet.
#'
#' @param data_dir Directory with raw RNA-Seq files.
#' @param pdata_path Path to text file with sample annotations.
#' @inheritParams build_ensdb
#'
#' @return
#' @export
#'
#' @examples
load_seq <- function(data_dir, pdata_path, species = 'Homo sapiens', release = '94') {

  # load EnsDb package
  ensdb_species    <- strsplit(species, ' ')[[1]]
  ensdb_species[1] <- substr(ensdb_species[1], 1, 1)

  ensdb_package <- paste('EnsDb', paste0(ensdb_species, collapse = ''), paste0('v', release), sep='.')

  if (!require(ensdb_package, character.only = TRUE)) {
    build_ensdb(species, release)
    require(ensdb_package, character.only = TRUE)
  }

}

#' Match and add file names for pdata rows.
#'
#' Attempts to find exact match between file names (excluding fastq.gz suffix) and first pdata column.
#' If an exact match is not possible, \code{grep} is used to search for a unique substring.
#'
#' @param pdata \code{data.table} with sample annotations.
#' @param files Character vector of file names to find matches for.
#'
#' @return \code{pdata} with \code{file} column identifying file for each sample (row).
#' @export
#'
#' @examples
match_pdata <- function(pdata, files) {
  # TODO: write some unit tests

  ids <- pdata[[1]]

  if (length(ids) != length(files)) stop('Must be one row in sample annotation per raw data file.')

  # check if every id has a perfect match in files without fastq.gz suffix
  idxs <- match(ids, gsub('.fastq.gz$', '', files))

  # otherwise try unique grep
  if (sum(is.na(idxs)) != 0) {

    # remove _ and - as word seperation (for grep below)
    spaced_ids <- gsub('_|-', ' ', ids)
    spaced_files <- gsub('_|-', ' ', files)

    tomatch <- spaced_ids[is.na(idxs)]
    gmatch  <- sapply(tomatch, function(id) {

      res <- grep(paste0('\\b', id, '\\b'), spaced_files, ignore.case = TRUE)

      if (length(res) > 1)
        stop(shQuote(id), ' is a substring of more than one raw rna-seq data file. Fix in first column of sample annotation file.')

      if (!length(res))
        stop(shQuote(id), ' is not a substring of a raw rna-seq data file. Fix in first column of sample annotation file.')

      return(res)
    })

    # fill in grep result where NA in idxs
    idxs[is.na(idxs)] <- gmatch[is.na(idxs)]
  }

  # append file names to pdata
  pdata$file <- files[idxs]
  return(pdata)
}
