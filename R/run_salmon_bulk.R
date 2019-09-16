#' Runs salmon quantification.
#'
#' For pair-ended experiments, reads for each pair should be in a seperate file.
#'
#' @param indices_dir Directory with salmon indices. See \code{\link{build_ensdb_index}}.
#' @param data_dir Directory with raw fastq.gz RNA-Seq files.
#' @param pdata Previous result of call to \code{run_salmon} or \code{\link{select_pairs}}. Used to bypass another call to \code{select_pairs}.
#' @param species Species name. Default is \code{homo_sapiens}.
#' Used to determine transcriptome index to use.
#' @param flags Character vector of flags to pass to salmon.
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' # first place IBD data in data-raw/example-data
#' indices_dir <- '~/Documents/Batcave/zaklab/drugseqr.data/inst/indices'
#' data_dir <- file.path('data-raw', 'example-data')
#' pdata_path <- file.path(data_dir, 'Phenotypes.csv')
#' run_salmon_bulk(indices_dir, data_dir, pdata_path)
#'
run_salmon_bulk <- function(indices_dir, data_dir, pdata = NULL, species = 'homo_sapiens',
                       flags = c('--validateMappings', '--posBias', '--seqBias', '--gcBias')) {
  # TODO: make it handle single and paired-end data
  # now assumes single end

  species <- gsub(' ', '_', tolower(species))

  # location of index
  salmon_version <- get_pkg_version('salmon')
  salmon_idx <- file.path(indices_dir, paste0('salmon_', salmon_version), 'ensdb', species)
  if (!dir.exists(salmon_idx)) stop('No index found. See ?drugseqr.data::build_salmon_index')

  # prompt user to select pairs/validate file names etc
  if (is.null(pdata)) pdata <- select_pairs(data_dir)
  pdata$quants_dir <- gsub('.fastq.gz$', '', pdata$`File Name`)

  # save selections
  saveRDS(pdata, file.path(data_dir, 'pdata.rds'))

  # save quants here
  quants_dir <- file.path(data_dir, paste0('salmon_quants_', salmon_version))
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

#' Get version of salmon/kallisto from system command.
#'
#' @param type Either \code{'salmon'} or \code{'kallisto'}.
#'
#' @return Version of salmon/kallisto.
#' @export
#' @keywords internal
#'
#' @examples
get_pkg_version <- function(type) {
  # possibly use older salmon with version appended to executable name
  if (type == 'salmon') {
    version <- system(paste(type, '--version'), intern = TRUE)
    version <- gsub('^salmon ', '', version)

  } else if (type == 'kallisto') {
    version <- system(paste(type, 'version'), intern = TRUE)
    version <- gsub('^kallisto, version ', '', version)

  } else {
    stop('type must be either salmon or kallisto')
  }

  return(version)
}
