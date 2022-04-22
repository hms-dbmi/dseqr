#' Run kallisto/bustools for quantifying 10X scRNA-seq data
#'
#' @param indices_dir Directory with kallisto indices built with \code{build_kallisto_index}.
#' @param data_dir Path to folder with 10X fastq.gz scRNA-seq files
#' @param out_dir Path to save output to.
#' @param species Species to use index from. Currently supports 'human' or 'mouse'.
#' @param threads Number of threads to use. Default is 1.
#'
#' @return called for side effects.
#' @seealso \link{download_kb_index}
#' @export
#'
run_kb_scseq <- function(indices_dir, data_dir, out_dir = file.path(data_dir, 'bus_output'), species = 'human', threads = 4) {

  # kallisto needs expanded paths (no tilde's)
  indices_dir <- path.expand(indices_dir)
  out_dir <- path.expand(out_dir)

  # get index_path
  kb_version <- get_kb_version()
  index_dir <- file.path(path.expand(indices_dir),
                         paste0('kb_', kb_version),
                         species)

  # get CB and read fastq files
  cb_fastqs <- identify_sc_files(data_dir, 'R1')
  read_fastqs <- identify_sc_files(data_dir, 'R2')

  if (length(cb_fastqs) != length(read_fastqs))
    stop("Different number of cell-barcode and read fastqs.")

  # alternate CB and read fastqs
  fqs <- c(rbind(cb_fastqs, read_fastqs))

  # detect v1/v2/v3 chemistry
  chemistry <- detect_10x_chemistry(index_dir, data_dir, fqs)
  message("\nChemistry detected: ", chemistry)

  # run quantification/bustools
  run_kb_count(index_dir, out_dir, file.path(data_dir, fqs), chemistry, threads = threads)
  return(NULL)
}

# get kb-python version
get_kb_version <- function() {
  version <- suppressWarnings(system2("kb", stderr = TRUE))
  version <- version[grepl('kb_python', version)]
  version <- gsub("^kb_python ", "", version)
  return(version)
}

#' Download pre-built index using kb-python
#'
#' Downloads index.idx and t2g.txt.
#'
#' @param indices_dir directory with indices. Output files will be saved in
#' indices_dir/kb_version/species/
#' @param species either 'human' or 'mouse'
#'
#' @return Called for side effects
#' @export
#'
download_kb_index <- function(indices_dir, species = c('human', 'mouse')) {
  kb_version <- get_kb_version()
  index_dir <- file.path(path.expand(indices_dir),
                         paste0('kb_', kb_version),
                         species[1])

  dir.create(index_dir, recursive = TRUE)

  system2("kb",
          args = c(
            "ref -d", species,
            "-i", shQuote(file.path(index_dir, 'index.idx')),
            "-g", shQuote(file.path(index_dir, 't2g.txt'))
          ))
}


#' Auto Detect 10x Chemistry
#'
#' Creates samples of 10,000 reads per lane and then determins the number of reads
#' with barcodes in agreement with the chemistry whitelist.
#'
#' @param index_dir Path to kallisto index.
#' @param data_dir Path to folder with 10x fastq.gz files.
#' @param bus_args Additional arguments to kallisto bus.
#' @param fqs Vector of fastq file names.
#' @param techs 10x chemistries to check. Passed to kallisto -x argument.
#' @param threads Number of threads to use. Default is 1.
#'
#' @return one of \code{techs} corresponding to chemistry with most reads that agree with whitelist.
#' @keywords internal
detect_10x_chemistry <- function(index_dir, data_dir, fqs, techs = c('10XV2', '10XV3'), threads = 1) {

  # run quanitification on first 10000 reads
  out_dir <- file.path(data_dir, 'tmp_output')
  fqn <- file.path(data_dir, paste0('N=10000_', fqs))
  fqs <- file.path(data_dir, fqs)

  # generate sample fastqs for detection
  for (i in seq_along(fqs))
    system(paste('zcat', shQuote(fqs[i]), '| head -n 40000 | gzip >', shQuote(fqn[i])))

  # run kallisto on samples and get number of reads with barcode in agreement with whitelist
  nreads <- c()
  for (tech in techs) {
    message('\nTesting if chemistry is: ', tech)

    ntech <- 0
    try({
      run_kb_count(index_dir, out_dir, fqn, tech, threads, stderr = NULL)
      info <- jsonlite::fromJSON(file.path(out_dir, 'inspect.json'))
      ntech <- info$numReadsOnWhitelist
    }, silent = TRUE)


    message('Reads on whitelist: ', ntech)
    nreads <- c(nreads, ntech)
    unlink(out_dir, recursive = TRUE)
  }

  # return tech with largest percent
  unlink(fqn)
  if (all(nreads == 0)) stop("Could not detect 10X chemistry.")

  techs[which.max(nreads)]
}

# runs kb count from kb-python
run_kb_count <- function(index_dir, out_dir, fqs, chemistry, threads = 1, stderr = '') {

  system2("kb",
          args = c(
            "count",
            "-i", shQuote(file.path(index_dir, 'index.idx')),
            "-g", shQuote(file.path(index_dir, 't2g.txt')),
            "-o", shQuote(out_dir),
            "-x", chemistry,
            "-t", threads,
            "--tmp", file.path(tempdir(), 'kb_tmp'),
            "--cellranger",
            shQuote(fqs)
          ), stderr = stderr)
}

#' Get 10X FastQ paths for a read type and sort by increasing lane number
#'
#' @inheritParams run_kb_scseq
#' @param read_type One of \code{'R1'} or \code{'R2'} for CB+UMI and read sequences respectively.
#'
#' @return Character vector of lane ordered paths to FastQ read type files.
#' @keywords internal
identify_sc_files <- function(data_dir, read_type = 'R1') {
  # see https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input for description of file name format

  # get single cell files for read type
  sc_files <- list.files(data_dir, paste0('_', read_type, '_.+?.fastq.gz'))

  if (!length(sc_files))
    stop('No fastq.gz files identified in data_dir in read type of ', read_type)

  if (length(sc_files) == 1) return(sc_files)

  # order by increasing lane if multiple
  # first try format from mkfastq
  lanes <- gsub('^.+?_L([0-9]+)_.+?.fastq.gz$', '\\1', sc_files)

  # try format from (deprecated)
  if (!length(lanes))
    lanes <- gsub('^.+?_lane-([0-9]+)-.+?.fastq.gz$', '\\1', sc_files)

  if (!length(lanes))
    stop('10X fastq.gz file name formats not recognized. Should be output by mkfastq or possibly demux (deprecated).')

  sc_files <- sc_files[order(as.integer(lanes))]
  return(sc_files)
}

