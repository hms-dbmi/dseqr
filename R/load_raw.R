#' Build ensembldb annotation package.
#'
#' @param species Character vector indicating species. Should be space seperated, not underscore. Default is \code{Homo sapiens}.
#' @param version EnsemblDB version. Should be same as used in \code{\link{build_index}}.
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' # build ensembldb annotation package for human
#' build_ensdb()
#'
build_ensdb <- function(species = 'Homo sapiens', version = '94') {
  # store ensembl databases in built package
  ensdb_dir <- system.file('EnsDb', package = 'drugseqr')

  # format is genus_species in multiple other functions but not here
  species <- gsub('_', ' ', species)

  # generate new ensembl database from specified version
  ah <- AnnotationHub::AnnotationHub()
  ahDb <- AnnotationHub::query(ah, pattern = c(species, "EnsDb", version))

  if (!length(ahDb)) stop('Specified ensemble species/version not found in AnnotationHub.')

  ahEdb <- ahDb[[1]]

  ensembldb::makeEnsembldbPackage(AnnotationDbi::dbfile(ensembldb::dbconn(ahEdb)),
                                  '0.0.1', 'Alex Pickering <alexvpickering@gmail.com>',
                                  'Alex Pickering',
                                  ensdb_dir)

  # build new ensemble database
  ensdb_name <- list.files(ensdb_dir)
  ensdb_path <- file.path(ensdb_dir, ensdb_name)

  work_dir <- getwd()
  setwd(ensdb_dir)
  system(paste('R CMD build', ensdb_name))

  # install new ensemble database
  ensdb_build <- paste0(ensdb_name, '_0.0.1.tar.gz')
  system(paste('R CMD INSTALL', ensdb_build))

  # remove build file
  setwd(work_dir)
  unlink(ensdb_dir, recursive = TRUE)
}
