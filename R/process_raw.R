#' Download ensembl transcriptome and build index for salmon quantification
#'
#' @param species The species. Default is \code{homo_sapiens.}
#' @param version ensembl version. Default is \code{94} (latest in release for AnnotationHub - needs to match with \code{\link{build_ensdb}}).
#'
#' @return
#' @export
#'
#' @examples
build_index <- function(species = 'homo_sapiens', version = '94') {

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
  tryCatch(system2('salmon', args=c('index',
                                    '-t', ensembl_all,
                                    '-i', ensembl_species)),
           error = function(e) stop('Is salmon installed and on the PATH?'))

  unlink(ensembl_all)
  setwd(work_dir)
}

#' Runs salmon quantification.
#'
#' @param data_dir Directory with .fastq.gz RNA-Seq files.
#' @param species Species name. Default is \code{homo_sapiens}.
#' Used to determine transcriptome index to use.
#'
#' @return
#' @export
#'
#' @examples
run_salmon <- function(data_dir, species = 'homo_sapiens') {
  # TODO: make it handle single and paired-end data
  # now assumes single end

  # location of index
  salmon_idx <- system.file('indices', species, package = 'drugseqr')
  if (!dir.exists(salmon_idx)) stop('No index found. See ?build_index')

  # save quants here
  quants_dir <- file.path(data_dir, 'quants')
  dir.create(quants_dir)

  # loop through fastq files and quantify
  fastq_files <- list.files(data_dir, '.fastq.gz$')
  for (fastq_file in fastq_files) {
    fastq_path <- shQuote(file.path(data_dir, fastq_file))

    # save each sample in it's own folder
    sample_name <- gsub('.fastq.gz$', '', fastq_file)
    out_dir <- file.path(quants_dir, sample_name)
    dir.create(out_dir)

    # run salmon
    system2('salmon',
            args=c('quant',
                   '-i', salmon_idx,
                   '-l', 'A',
                   '--validateMappings',
                   '-r', fastq_path, # single-end flag
                   '--gcBias',
                   '-o', out_dir))
  }
}



