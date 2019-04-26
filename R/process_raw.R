#' Download ensembl transcriptome and build index for salmon quantification
#'
#' @param species The species. Default is \code{homo_sapiens.}
#' @param release ensembl release. Default is \code{94} (latest in release for AnnotationHub - needs to match with \code{\link{build_ensdb}}).
#'
#' @return NULL
#' @export
#'
#' @examples
#' # build salmon index for humans
#' build_index()
#'
build_index <- function(species = 'homo_sapiens', release = '94') {

  indices_dir <- system.file('indices', package = 'drugseqr')

  # construct ensembl url for transcriptome
  ensembl_species <- gsub(' ', '_', tolower(species))
  ensembl_release <- paste0('release-', release)
  ensembl_url <- paste0('ftp://ftp.ensembl.org/pub/', ensembl_release, '/fasta/', ensembl_species, '/cdna/')

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
           error = function(err) {err$message <- 'Is salmon installed and on the PATH?'; stop(err)})

  unlink(ensembl_all)
  setwd(work_dir)
}

#' Runs salmon quantification.
#'
#' For pair-ended experiments, reads for each pair should be in a seperate file.
#'
#' @param data_dir Directory with raw RNA-Seq fastq.gz files.
#' @param species Species name. Default is \code{homo_sapiens}.
#' Used to determine transcriptome index to use.
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' # first place IBD data in data-raw/example-data
#' data_dir <- file.path('data-raw', 'example-data')
#' run_salmon(data_dir)
#'
run_salmon <- function(data_dir, species = 'homo_sapiens') {
  # TODO: make it handle single and paired-end data
  # now assumes single end

  # location of index
  salmon_idx <- system.file('indices', species, package = 'drugseqr')
  if (!dir.exists(salmon_idx)) stop('No index found. See ?build_index')

  # save quants here
  quants_dir <- file.path(data_dir, 'quants')
  dir.create(quants_dir)

  fastq_files <- list.files(data_dir, '.fastq.gz$')

  # check if pair-ended
  fastq_paths <- file.path(data_dir, fastq_files)
  fastq_id1s <- get_fastq_id1s(fastq_paths)
  paired <- detect_paired(fastq_id1s)
  if (paired) stop('Have not yet implemented pair-ended experiments. Please contact author.')


  # loop through fastq files and quantify
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
                   '-o', shQuote(out_dir)))
  }
}

#' Get first sequence identifiers for fastq.gz files
#'
#' Function is used to determine if an experiment is single or paired.
#'
#' @param fastq_paths Character vector of paths to fastq.gz files
#'
#' @return Named character vector of first sequence id lines (start with @) in \code{fastq_paths}.
#'   Names are the \code{fastq_paths}.
#' @keywords internal
#' @export
#'
#' @example
#' # required fastq.gz files in data-raw/example-data (e.g. IBD example)
#' fastq_paths <- list.files(file.path('data-raw', 'example-data'), '.fastq.gz$', full.names = TRUE)
#' fastq_id1s <- get_fastq_id1s(fastq_paths)
#'
get_fastq_id1s <- function(fastq_paths) {

  # get first line with @ symbol (sequence identifier)
  fastq_id1s <- sapply(fastq_paths, function(f) {
    incon <- gzfile(f)
    while (TRUE) {
      line = readLines(incon, n = 1)
      if (grepl('^@', line)) {
        break
      }
    }
    close(incon)
    return(line)
  })
  return(fastq_id1s)
}

#' Detect if experiment is pair-ended.
#'
#' @param fastq_id1s Character vector of first sequence identifiers from fastq.gz files. Returned from \code{\link{get_fastq_id1s}}.
#'
#' @return boolean indicating if experiement is pair-ended (\code{TRUE}) or single-ended (\code{FALSE}).
#' @export
#'
#' @examples
detect_paired <- function(fastq_id1s) {

  # TODO: handle SRA fastq files
  # for now not implemented
  if (any(grepl('^@SRR\\d+', fastq_id1s)))
    stop('Detecting if SRA fastq files are paired is not yet implemented.')

  # older illumina sequence identifiers have 1 part
  # newer illumina sequence identifiers have 2 space-seperated parts
  id_parts <- strsplit(fastq_id1s, ' ')
  older <- all(sapply(id_parts, length) == 1)
  newer <- all(sapply(id_parts, length) == 2)

  if (older) {
    # pair is 1 or 2 at end of sequence id after /
    pairs <- gsub('^.+?/([12])$', '\\1', fastq_id1s)

  } else if (newer) {
    # pair is 1 or 2 followed by : followed by N or Y at beginning of second part
    id_parts2 <- sapply(id_parts, `[`, 2)
    pairs <- gsub('^([12]):[YN]:.+$', '\\1', id_parts2)

  } else {
    stop("fastq.gz files don't appear to be from older/newer Illumina software. Please contact package author.")
  }

  # paired experiments will have '1' and '2'
  paired <- length(unique(pairs)) == 2
  return(paired)
}



