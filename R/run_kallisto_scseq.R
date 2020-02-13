#' Run kallisto/bustools for quantifying 10X scRNA-seq data
#'
#' @param indices_dir Directory with kallisto indices.
#' @param data_dir Path to folder with 10X fastq.gz scRNA-seq files
#' @param bus_args Character vector of arguments to bustools.
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' indices_dir <- 'data-raw/indices'
#' data_dir <- 'data-raw/single-cell/example-data/Run2643-10X-Lung/10X_FID12518_Diseased_3hg'
#'
#' run_kallisto_scseq(indices_dir, data_dir)
#'
run_kallisto_scseq <- function(indices_dir, data_dir, bus_args = '-t 4', species = 'homo_sapiens', release = '94', recount = TRUE) {

  out_dir <- file.path(data_dir, 'bus_output')
  if (dir.exists(out_dir) & !recount) return(NULL)

  # make sure that have whitelist and get path
  dl_10x_whitelists(indices_dir)

  # get index_path
  kal_version <- get_pkg_version('kallisto')
  index_path <- file.path(indices_dir, paste0('kallisto_', kal_version), paste0(species, '.grch38.cdna.all.release-', release, '_k31.idx'))

  # get CB and read fastq files
  cb_fastqs <- identify_sc_files(data_dir, 'R1')
  read_fastqs <- identify_sc_files(data_dir, 'R2')

  if (length(cb_fastqs) != length(read_fastqs))
    stop("Different number of cell-barcode and read fastqs.")

  # alternate CB and read fastqs
  fqs <- c(rbind(cb_fastqs, read_fastqs))

  # detect v1/v2/v3 chemistry
  chemistry <- detect_10x_chemistry(indices_dir, index_path, data_dir, bus_args, fqs)
  whitepath <- file.path(indices_dir, paste0(chemistry, '_whitelist.txt'))

  # run quantification/bustools
  bus_args <- c('bus',
                '-i', index_path,
                '-o', out_dir,
                '-x', chemistry,
                bus_args,
                file.path(data_dir, fqs))

  run_kallisto_scseq_commands(bus_args, whitepath, out_dir)
  return(NULL)
}

#' Runs kallisto scseq commands
#'
#' @param bus_args arguments to \code{system2} for kallisto bus
#' @param whitepath Path to whitelist
#' @param out_dir Directory to write to
#' @param inspection Boolean indicating if kallisto inspect should be run. Used for \code{detect_10x_chemistry}.
#'
#' @return stdout from kallisto inspect if \code{inspection = TRUE}, otherwise \code{NULL}.
#' @export
run_kallisto_scseq_commands <- function(bus_args, whitepath, out_dir, inspection = FALSE) {


  system2('kallisto', args = bus_args)

  # run bustools
  tmp_dir <- file.path(out_dir, 'tmp')
  gct_dir <- file.path(out_dir, 'genecount', 'genes')
  dir.create(tmp_dir)
  dir.create(gct_dir, recursive = TRUE)

  # inspection used to auto-detect chemistry
  if (inspection) {
    out <- system2(
      'bustools',
      args = c('correct',
               '-w', whitepath,
               '-p', file.path(out_dir, 'output.bus'),
               '| bustools sort -T', tmp_dir,
               '-t 4 -p - | bustools inspect -w', whitepath,
               '-p -'),
      stdout = TRUE)
    return(out)

  } else {

    # location of map from transcript names to gene names
    tgmap_path <- system.file('extdata', 'txp2hgnc.tsv', package = 'drugseqr', mustWork = TRUE)

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
  }

}

#' Auto Detect 10x Chemistry
#'
#' Creates samples of 10,000 reads per lane and then determins the number of reads
#' with barcodes in agreement with the chemistry whitelist.
#'
#' @param indices_dir Folder that containes whitelist for each \code{tech}.
#' @param index_path Path to kallisto index.
#' @param data_dir Path to folder with 10x fastq.gz files.
#' @param bus_args Additional arguments to kallisto bus.
#' @param fqs Vector of fastq file names.
#' @param techs 10x chemistries to check. Passed to kallisto -x argument.
#'
#' @return one of \code{techs} corresponding to chemistry with most reads that agree with whitelist.
#' @export
detect_10x_chemistry <- function(indices_dir, index_path, data_dir, bus_args, fqs, techs = c('10xv2', '10xv3')) {

  # run quanitification on first 10000 reads
  out_dir <- file.path(data_dir, 'tmp_output')
  fqn <- file.path(data_dir, paste0('N=10000_', fqs))
  fqs <- file.path(data_dir, fqs)

  # generate sample fastqs for detection
  for (i in seq_along(fqs))
    system(paste('zcat', fqs[i], '| head -n 40000 | gzip >', fqn[i]))

  # run kallisto on samples and get number of reads with barcode in agreement with whitelist
  nreads <- c()
  for (tech in techs) {
    whitepath <- file.path(indices_dir, paste0(tech, '_whitelist.txt'))

    bus_argsi <- c('bus',
                   '-i', index_path,
                   '-o', out_dir,
                   '-x', tech,
                   bus_args,
                   fqn)

    out <- run_kallisto_scseq_commands(bus_argsi, whitepath, out_dir, inspection = TRUE)
    nread <- gsub('^.+? whitelist: (\\d+) .+?$', '\\1', out[16])
    nreads <- c(nreads, as.numeric(nread))
    unlink(out_dir, recursive = TRUE)
  }

  # return tech with largest percent
  unlink(fqn)
  techs[which.max(nreads)]
}

#' Get path to 10x whitelist
#'
#'
#' @inheritParams run_kallisto_scseq
#'
#' @return Path to 10x whitelist.
#' @export
#' @keywords internal
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
#' @return NULL
#' @export
#' @keywords internal
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
  return(NULL)
}

#' Get 10X FastQ paths for a read type and sort by increasing lane number
#'
#' @inheritParams run_kallisto_scseq
#' @param read_type One of \code{'R1'} or \code{'R2'} for CB+UMI and read sequences respectively.
#'
#' @return Character vector of lane ordered paths to FastQ read type files.
#' @export
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
