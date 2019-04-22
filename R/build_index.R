#' Download ensembl transcriptome and build index for salmon quantification
#'
#' @param species The species. Default is \code{homo_sapiens.}
#' @param version ensembl version. Default is \code{96} (latest as of 04-22-2019).
#'
#' @return
#' @export
#'
#' @examples
build_index <- function(species = 'homo_sapiens', version = '96') {

  indices_dir <- system.file('indices', package = 'drugseqr')

  # construct ensembl url for transcriptome
  ensembl_species <- gsub(' ', '_', tolower(species))
  ensembl_version <- paste0('release-', version)
  ensembl_url <- paste0('ftp://ftp.ensembl.org/pub/', ensembl_version, '/fasta/', ensembl_species, '/cdna/')

  # get list of all files
  handle <- curl::new_handle(dirlistonly=TRUE)
  con <- curl::curl(ensembl_url, "r", handle)
  tbl <- utils::read.table(con, stringsAsFactors=FALSE)
  close(con)

  # get transcripts cdna.all file
  ensembl_all <- grep('cdna.all.fa.gz$', tbl[[1]], value = TRUE)
  ensembl_url <- paste0(ensembl_url, ensembl_all)

  work_dir <- getwd()
  setwd(indices_dir)
  curl::curl_download(ensembl_url, ensembl_all)

  # build index
  tryCatch(system2(paste('salmon index -t', ensembl_all, '-i', ensembl_species), stdout=TRUE, stderr=TRUE),
           error = function(e) stop('Is salmon installed and on the PATH?'))

  unlink(ensembl_all)
  setwd(work_dir)
}



