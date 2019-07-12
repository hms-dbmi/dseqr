#' Runs kallisto quantification for bulk samples.
#'
#' For pair-ended experiments, reads for each pair should be in a seperate file.
#'
#' @param indices_dir Directory with kallisto indices. See \code{\link{build_kallisto_index}}.
#' @param data_dir Directory with raw fastq.gz RNA-Seq files.
#' @param pdata Previous result of call to \code{run_kallisto_bulk} or \code{\link{select_pairs}}. Used to bypass another call to \code{select_pairs}.
#' @param species Species name. Default is \code{homo_sapiens}.
#' Used to determine transcriptome index to use.
#' @param fl.mean Estimated average fragment length (only relevant for single-end reads).
#' @param fl.sd Estimated standard deviation of fragment length (only relevant for single-end reads).
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' # first place IBD data in data-raw/example-data
#' indices_dir <- 'data-raw/indices/kallisto'
#' data_dir <- file.path('data-raw', 'example-data')
#' pdata_path <- file.path(data_dir, 'Phenotypes.csv')
#' run_kallisto_bulk(indices_dir, data_dir, pdata_path)
#'
run_kallisto_bulk <- function(indices_dir, data_dir, pdata = NULL, species = 'homo_sapiens', release = '94', fl.mean = NULL, fl.sd = NULL) {
  # TODO: make it handle single and paired-end data
  # now assumes single end

  # get index_path
  index_path <- file.path(indices_dir, paste0(species, '.grch38.cdna.all.release-', release, '_k31.idx'))

  # prompt user to select pairs/validate file names etc
  if (is.null(pdata)) pdata <- select_pairs(data_dir)
  pdata$quants_dir <- gsub('.fastq.gz$', '', pdata$`File Name`)

  # save selections
  saveRDS(pdata, file.path(data_dir, 'pdata.rds'))

  # save quants here
  quants_dir <- file.path(data_dir, 'quants')
  unlink(quants_dir, recursive = TRUE)
  dir.create(quants_dir)

  # specific flags for single end experiments
  paired <- 'Pair' %in% colnames(pdata)
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

    # remove fastq_files for next quant loop
    fastq_files <- setdiff(fastq_files, fastq_file)

    # possibly spaces in file names
    fastq_path <- paste(shQuote(file.path(data_dir, fastq_file)), collapse = ' ')

    # save each sample in it's own folder
    # use the first file name in any replicates
    sample_name <- gsub('.fastq.gz$', '', fastq_file[1])
    out_dir <- file.path(quants_dir, sample_name)

    # run kallisto
    system2('kallisto',
            args=c('quant',
                   '-i', index_path,
                   '-o', shQuote(out_dir),
                   flags,
                   fastq_path))
  }
  return(NULL)
}
