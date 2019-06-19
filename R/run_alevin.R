
#' Runs salmon alevin single-cell quantification.
#'
#' @param indices_dir Directory with salmon indices. See \code{\link{build_gencode_index}}.
#' @param data_dir Directory with raw 10X single-cell fastq.gz RNA-Seq files.
#' @param species Species name. Default is \code{human}.
#' Used to determine transcriptome index to use.
#' @param overwrite Do you want to overwrite results of previous run? Default is \code{FALSE}.
#' @param command System command to invoke salmon. Can be used to invoked different versions of salmon.
#' @param soup Should alevin be run to obtain ambient counts for \code{SoupX}? In development.
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' indices_dir <- 'data-raw/indices'
#'
#' # first place data in data-raw/single-cell/example-data
#' data_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'
#' run_alevin(data_dir, indices_dir)
#'
#' # another sample
#' data_dir <- 'data-raw/single-cell/example-data/Run2643-10X-Lung/10X_FID12518_Diseased_3hg'
#' run_alevin(data_dir, indices_dir)
#'
run_alevin <- function(data_dir, indices_dir, species = 'human', overwrite = FALSE, command = 'salmon', soup = FALSE) {

  # possibly use older salmon with version appended to executable name
  salmon_version <- get_salmon_version(command)

  # need .mtx if using salmon 0.14.0 until tximport supports
  flags <- c()
  if(salmon_version == '0.14.0') flags <- '--dumpMtx'

  # make sure its 10X
  sc_method <- detect_sc_method(data_dir)

  # location of tgMap
  tgmap_path <- system.file('extdata', 'txp2hgnc.tsv', package = 'drugseqr')

  # location of index
  alevin_idx <- file.path(indices_dir, salmon_version, 'gencode', species)
  if (!dir.exists(alevin_idx)) stop('No index found. See ?build_gencode_index')

  # save alevin output here
  out_dir <- file.path(data_dir, paste0('alevin_output_', salmon_version))

  if (soup) {
    soup_whitelist_path <- make_soup_whitelist(out_dir)
    out_dir <- paste0(out_dir, '_soup')
    flags <- c(flags, '--whitelist', soup_whitelist_path, '--maxNumBarcodes', 4294967295, '--freqThreshold', 0)

  } else {
    # location of ribosomal/mitochondrial gene files (used for whitelist model)
    rrna_path <- system.file('extdata', 'rrna.csv', package = 'drugseqr')
    mrna_path <- system.file('extdata', 'mrna.csv', package = 'drugseqr')

    flags <- c(flags, '--mrna', mrna_path, '--rrna', rrna_path)
  }

  # make sure not overwriting previous result by mistake
  if (file.exists(file.path(out_dir, 'alevin', 'quants_mat.gz')) & !overwrite)
    stop('Output from previous alevin run already exists. Use overwrite = TRUE to re-run analysis.')

  unlink(out_dir, recursive = TRUE)

  # get R1 (CB+UMI sequences) and R2 (read sequences) order by lane number
  cb_fastqs <- identify_sc_files(data_dir, 'R1')
  read_fastqs <- identify_sc_files(data_dir, 'R2')

  if (length(cb_fastqs) != length(read_fastqs))
    stop("Detected different number of cell barcode and read sequence FastQs.")


  # run salmon alevin
  system2(command,
          args=c('alevin',
                 '-l', 'ISR',
                 '-1', paste(cb_fastqs, collapse = ' '),
                 '-2', paste(read_fastqs, collapse = ' '),
                 '--chromium',
                 '-i', alevin_idx,
                 '-o', shQuote(out_dir),
                 '-p', 6,
                 '--dumpFeatures',
                 flags,
                 '--tgMap', tgmap_path))

  return(NULL)
}

make_soup_whitelist <- function(out_dir) {
  freq_path <- file.path(out_dir, 'alevin', 'raw_cb_frequency.txt')

  if (!file.exists(freq_path))
    stop('Run alevin once with soup=FALSE first is required to quantify soup.')

  freqs <- data.table::fread(freq_path, col.names = c('barcode', 'umis'))
  soup_whitelist <- freqs[umis > 0 & umis <= 10, .(barcode)]
  soup_whitelist <- soup_whitelist[sample(.N, min(.N, 50000))]

  soup_whitelist_path <- file.path(out_dir, 'alevin', 'soup_whitelist.txt')

  data.table::fwrite(soup_whitelist, soup_whitelist_path, col.names = FALSE)
  return(soup_whitelist_path)
}

#' Determine Method used for Single Cell RNA Seq
#'
#' Currently just looks for 10X data.
#'
#' @inheritParams run_alevin
#'
#' @return Character vector of detected method.
#' @export
#'
#' @examples
detect_sc_method <- function(data_dir) {
  fastq_paths <- list.files(data_dir, '.fastq.gz', full.names = TRUE)
  is_10x <- any(grepl('10X', fastq_paths))

  if (!is_10x) stop('only 10X implemented - 10X should be somewhere in the fastq.gz file name or path')
  return('10X')
}

#' Get 10X FastQ paths for a read type and sort by increasing lane number
#'
#' @inheritParams run_alevin
#' @param read_type One of \code{'R1'} or \code{'R2'} for CB+UMI and read sequences respectively.
#'
#' @return Character vector of lane ordered paths to FastQ read type files.
#' @export
#'
#' @examples
identify_sc_files <- function(data_dir, read_type = 'R1') {
  # see https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input for description of file name format

  # get single cell files for read type
  sc_files <- list.files(data_dir, paste0('_', read_type, '_.+?.fastq.gz'))

  if (!length(sc_files))
    stop('No fastq.gz files identified in data_dir in read type of ', read_type)

  if (length(sc_files) == 1) return(file.path(data_dir, sc_files))

  # order by increasing lane if multiple
  # first try format from mkfastq
  lanes <- gsub('^.+?_L([0-9]+)_.+?.fastq.gz$', '\\1', sc_files)

  # try format from (deprecated)
  if (!length(lanes))
    lanes <- gsub('^.+?_lane-([0-9]+)-.+?.fastq.gz$', '\\1', sc_files)

  if (!length(lanes))
    stop('10X fastq.gz file name formats not recognized. Should be output by mkfastq or possibly demux (deprecated).')

  sc_files <- sc_files[order(as.integer(lanes))]
  return(file.path(data_dir, sc_files))
}
