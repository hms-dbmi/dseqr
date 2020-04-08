#' Runs kallisto quantification for bulk samples.
#'
#' For pair-ended experiments, reads for each pair should be in a seperate file.
#'
#' @param indices_dir Directory with kallisto indices. See \code{\link{build_kallisto_index}}.
#' @param data_dir Directory with raw fastq.gz RNA-Seq files.
#' @param pdata Previous result of call to \code{run_kallisto_bulk} or \code{\link{select_pairs}}. Used to bypass another call to \code{select_pairs}.
#' @param species Species name. Default is \code{homo_sapiens}.
#' Used to determine transcriptome index to use.
#' @param fl.mean Estimated average fragment length (only relevant for single-end reads). Default (\code{NULL}) uses 200.
#' @param fl.sd Estimated standard deviation of fragment length (only relevant for single-end reads). Default (\code{NULL}) uses 20.
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' # first place IBD data in data-raw/example-data
#' indices_dir <- '/srv/drugseqr/indices'
#' data_dir <- file.path('data-raw', 'example-data')
#' run_kallisto_bulk(indices_dir, data_dir)
#'
run_kallisto_bulk <- function(indices_dir, data_dir, paired = NULL, pdata = NULL, species = 'homo_sapiens', release = '94', fl.mean = NULL, fl.sd = NULL, updateProgress = NULL) {

  # TODO: implement interaction between pairs and duplicates
  data_dir <- path.expand(data_dir)

  # default updateProgress and number of steps
  if (is.null(updateProgress)) updateProgress <- function(...) {NULL}

  # get index_path
  kallisto_version <- get_pkg_version('kallisto')
  index_path <- get_kallisto_index(indices_dir, species, release)

  if (is.null(pdata)) pdata <- select_pairs(data_dir)

  # save quants here
  quants_dir <- file.path(data_dir, paste('kallisto', kallisto_version, 'quants', sep = '_'))
  unlink(quants_dir, recursive = TRUE)
  dir.create(quants_dir)

  # if not supplied from app, then from select_pairs pdata
  if (is.null(paired)) paired <- sum(!is.na(pdata$Pair)) > 0

  # specific flags for single end experiments
  flags <- NULL
  if (!paired) {
    if (is.null(fl.mean)) {
      message('Single-end experiment but estimated average fragment length not provided. Setting to 200.')
      fl.mean <- 200
    }

    if (is.null(fl.sd)) {
      message('Single-end experiment but estimated standard deviation of fragment length not provided. Setting to 20.')
      fl.sd <- 20
    }
    flags <- c('--single', '-l', fl.mean, '-s', fl.sd)
  }

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
      fastq_file <- pdata[pdata$Replicate %in% rep_num, 'File Name', drop = TRUE]
    }

    # include any pairs
    pair_num <- pdata$Pair[row_num]
    if (!is.null(pair_num) && !is.na(pair_num)) {
      fastq_file <- pdata[pdata$Pair %in% pair_num, 'File Name', drop = TRUE]
    }

    # update progress
    updateProgress(amount = length(fastq_file))

    # remove fastq_files for next quant loop
    fastq_files <- setdiff(fastq_files, fastq_file)

    # possibly spaces in file names
    fastq_path <- paste(shQuote(file.path(data_dir, fastq_file)), collapse = ' ')

    # save each sample in it's own folder
    # use the first file name in any replicates/pairs
    sample_name <- gsub('.fastq.gz$', '', fastq_file[1])
    out_dir <- file.path(quants_dir, sample_name)

    # run kallisto
    system2('kallisto',
            args=c('quant',
                   '-i', index_path,
                   '-o', shQuote(out_dir),
                   '-t', 6,
                   flags,
                   fastq_path))
  }
  return(NULL)
}


#' Get path to kallisto index
#'
#' @param indices_dir Path to folder with indices dir
#' @param species Character vector giving species
#' @param release Character vector of EnsDB release number
#'
#' @return Path to kallisto index
#' @export
#' @keywords internal
#'
get_kallisto_index <- function(indices_dir, species = 'homo_sapiens', release = '94') {
  kallisto_version <- get_pkg_version('kallisto')
  species <- gsub(' ', '_', tolower(species))

  index_path <- file.path(indices_dir, paste0('kallisto_', kallisto_version))
  index_path <- list.files(index_path, paste0(species, '.+?.cdna.all.release-', release, '_k31.idx'), full.names = TRUE)

  return(index_path)
}
