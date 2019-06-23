#' Download ensembl transcriptome and build index for kalliso quantification
#'
#' This index is used for RNA seq quantification.
#'
#' @param species The species. Default is \code{homo_sapiens}.
#' @param release ensembl release. Default is \code{94} (latest in release for AnnotationHub -
#'   needs to match with \code{\link{build_ensdb}}).
#' @param indices_dir Directory to put kallisto indices.
#'
#' @return NULL
#' @export
#'
#' @examples
#' # build kallist index for humans
#' indices_dir <- 'data-raw/indices/kallisto'
#' build_ensdb_index(indices_dir)
#'
build_kallisto_index <- function(indices_dir, species = 'homo_sapiens', release = '94') {

  if (!dir.exists(indices_dir))
    dir.create(indices_dir, recursive = TRUE)

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
  ensembl_fasta <- grep('cdna.all.fa.gz$', tbl[[1]], value = TRUE)
  ensembl_url <- paste0(ensembl_url, ensembl_fasta)

  work_dir <- getwd()
  setwd(indices_dir)
  curl::curl_download(ensembl_url, ensembl_fasta)

  # build index
  index_fname <- gsub('fa.gz$', paste0(ensembl_release, '_k31.idx'), ensembl_fasta)
  index_fname <- tolower(index_fname)
  tryCatch(system2('kallisto', args=c('index',
                                   '-i', index_fname,
                                   ensembl_fasta)),
           error = function(err) {err$message <- 'Is kallisto installed and on the PATH?'; stop(err)})

  unlink(ensembl_fasta)
  setwd(work_dir)
}



#' Runs salmon quantification.
#'
#' For pair-ended experiments, reads for each pair should be in a seperate file.
#'
#' @param indices_dir Directory with salmon indices. See \code{\link{build_ensdb_index}}.
#' @param data_dir Directory with raw fastq.gz RNA-Seq files.
#' @param pdata_path Path to text file with sample annotations. Must be readable by \code{\link[data.table]{fread}}.
#' The first column should contain sample ids that match a single raw rna-seq data file name.
#' @param pdata Previous result of call to \code{run_salmon} or \code{\link{select_pairs}}. Used to bypass another call to \code{select_pairs}.
#' @param species Species name. Default is \code{homo_sapiens}.
#' Used to determine transcriptome index to use.
#' @param command System command to invoke salmon with.
#' @param flags Character vector of flags to pass to salmon.
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' # first place IBD data in data-raw/example-data
#' indices_dir <- 'data-raw/indices'
#' data_dir <- file.path('data-raw', 'example-data')
#' pdata_path <- file.path(data_dir, 'Phenotypes.csv')
#' run_salmon(indices_dir, data_dir, pdata_path)
#'
run_salmon <- function(indices_dir, data_dir, pdata_path = NULL, pdata = NULL, species = 'homo_sapiens', command = 'salmon',
                       flags = c('--validateMappings', '--posBias', '--seqBias', '--gcBias')) {
  # TODO: make it handle single and paired-end data
  # now assumes single end

  # location of index
  salmon_version <- get_salmon_version(command)
  salmon_idx <- file.path(indices_dir, salmon_version, 'ensdb', species)
  if (!dir.exists(salmon_idx)) stop('No index found. See ?build_ensdb_index')

  if (is.null(pdata_path) & is.null(pdata)) stop('One of pdata_path or pdata must be supplied.')

  # prompt user to select pairs/validate file names etc
  if (is.null(pdata)) pdata <- select_pairs(data_dir, pdata_path)
  pdata$quants_dir <- gsub('.fastq.gz$', '', pdata$`File Name`)

  # save selections
  saveRDS(pdata, file.path(data_dir, 'pdata.rds'))

  # save quants here
  quants_dir <- file.path(data_dir, paste0('quants_', salmon_version))
  unlink(quants_dir, recursive = TRUE)
  dir.create(quants_dir)

  # gcBias is experimental for single-end experiments
  paired <- 'Pair' %in% colnames(pdata)
  if (!paired) flags <- setdiff(flags, '--gcBias')

  # loop through fastq files and quantify
  fastq_files <- pdata$`File Name`
  if (is.null(fastq_files)) stop("File names should be in 'File Name' column")

  while (length(fastq_files)) {

    # grab next
    fastq_file <-fastq_files[1]
    row_num <- which(pdata$`File Name` == fastq_file)

    # include any replicates
    rep_num <- pdata$Replicate[row_num]
    if (!is.null(rep_num) && !is.na(rep_num)) {
      fastq_file <- pdata[pdata$Replicate %in% rep_num, 'File Name']
    }

    # remove fastq_files for next quant loop
    fastq_files <- setdiff(fastq_files, fastq_file)

    # possibly spaces in file names
    fastq_path <- paste(shQuote(file.path(data_dir, fastq_file)), collapse = ' ')

    # save each sample in it's own folder
    # use the first file name in any replicates
    sample_name <- gsub('.fastq.gz$', '', fastq_file[1])
    out_dir <- file.path(quants_dir, sample_name)
    dir.create(out_dir)

    # run salmon
    system2(command,
            args=c('quant',
                   '-i', salmon_idx,
                   '-l', 'A',
                   '-r', fastq_path, # single-end flag
                   flags,
                   '-o', shQuote(out_dir)))
  }
  return(NULL)
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
#' @examples
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



