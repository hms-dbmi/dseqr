#' Run kallisto/bustools for quantifying 10X scRNA-seq data
#'
#' @param indices_dir Directory with kallisto indices. See \code{\link{build_kallisto_index}}.
#' @param data_dir Path to folder with 10X fastq.gz scRNA-seq files
#' @param bus_args Character vector of arguments to bustools.
#' @inheritParams build_kallisto_index
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' indices_dir <- 'data-raw/indices/kallisto'
#' data_dir <- 'data-raw/single-cell/example-data/Run2643-10X-Lung/10X_FID12518_Diseased_3hg'
#'
#' run_kallisto_scseq(indices_dir, data_dir)
#'
run_kallisto_scseq <- function(indices_dir, data_dir, bus_args = c('-x 10xv2', '-t 4'), species = 'homo_sapiens', release = '94') {

  # make sure that have whitelist and get path
  dl_10x_whitelists(indices_dir)
  whitepath <- get_10x_whitepath(indices_dir, bus_args)

  # get index_path
  index_path <- file.path(indices_dir, paste0(species, '.grch38.cdna.all.release-', release, '_k31.idx'))

  # get CB and read fastq files
  cb_fastqs <- identify_sc_files(data_dir, 'R1')
  read_fastqs <- identify_sc_files(data_dir, 'R2')

  if (length(cb_fastqs) != length(read_fastqs))
    stop("Different number of cell-barcode and read fastqs.")

  # alternate CB and read fastqs
  fastqs <- c(rbind(cb_fastqs, read_fastqs))

  # run quantification
  out_dir <- file.path(data_dir, 'bus_output')

  system2('kallisto',
          args = c('bus',
                   '-i', index_path,
                   '-o', out_dir,
                   bus_args,
                   fastqs))

  # location of map from transcript names to gene names
  tgmap_path <- system.file('extdata', 'txp2hgnc.tsv', package = 'drugseqr')

  # run bustools
  tmp_dir <- file.path(out_dir, 'tmp')
  gct_dir <- file.path(out_dir, 'genecount', 'genes')
  dir.create(tmp_dir)
  dir.create(gct_dir, recursive = TRUE)

  system2('bustools',
          args = c('correct',
                   '-w', whitepath,
                   '-p', file.path(out_dir, 'output.bus'),
                   '| bustools sort -T', tmp_dir,
                   '-t 4 -p - | bustools count -o', gct_dir,
                   '-g', tgmap_path,
                   '-e', file.path(out_dir, 'matrix.ec'),
                   '-t', file.path(out_dir, 'transcripts.txt'),
                   '--genecounts -'))

  return(NULL)
}

#' Get path to 10x whitelist
#'
#'
#' @inheritParams run_kallisto_scseq
#'
#' @return Path to 10x whitelist.
#' @export
#' @keywords internal
#'
#' @examples
get_10x_whitepath <- function(indices_dir, bus_args) {
  if (any(grepl('-x 10xv2', bus_args))) {
    whitepath <- file.path(indices_dir, '10xv2_whitelist.txt')

  } else if (any(grepl('-x 10xv3', bus_args))) {
    whitepath <- file.path(indices_dir, '10xv3_whitelist.txt')

  } else {
    stop('only 10X v2 and v3 currently implemented')
  }

  return(whitepath)
}

#' Download 10X whitelists for kallisto bustools correction
#'
#' @inheritParams run_kallisto_scseq
#'
#' @return
#' @export
#' @keywords internal
#'
#' @examples
dl_10x_whitelists <- function(indices_dir) {

  whilelist_urls <- c(
    'https://github.com/bustools/getting_started/releases/download/getting_started/10xv2_whitelist.txt',
    'https://github.com/BUStools/getting_started/releases/download/species_mixing/10xv3_whitelist.txt')

  whitelist_paths <- file.path(indices_dir, c('10xv2_whitelist.txt', '10xv3_whitelist.txt'))

  # download what don't have
  for (i in seq_along(whitelist_paths)) {
    if (!file.exists(whitelist_paths[i]))
      system2('wget', args=c('--output-document', whitelist_paths[i],
                             whilelist_urls[i]))
  }
}

#' Get 10X FastQ paths for a read type and sort by increasing lane number
#'
#' @inheritParams run_kallisto_scseq
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
