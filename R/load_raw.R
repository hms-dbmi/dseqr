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
