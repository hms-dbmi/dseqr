#' Run kallisto/bustools for quantifying scRNA-seq
#'
#' Primarily for testing purposes. Recommendation is to use \code{run_alevin}.
#'
#' @param index_path Path to kallisto index file
#' @param whitelist_path path to 10X Chromium version 2 whitelist
#' @param data_dir Path to folder with fastq.gz scRNA-seq files
#'
#' @return NULL
#' @export
#'
#' @examples
run_kallisto <- function(index_path, whitelist_path, data_dir) {

  out_dir <- file.path(data_dir, 'bus_output')

  cb_fastqs <- identify_sc_files(data_dir, 'R1')
  read_fastqs <- identify_sc_files(data_dir, 'R2')

  # run quantification
  system2('/home/alex/miniconda3/envs/kallisto/bin/kallisto',
          args = c('bus',
                   '-i', index_path,
                   '-o', out_dir,
                   '-x', '10xv2',
                   '-t', 4,
                   cb_fastqs,
                   read_fastqs))

  # location of tgMap
  tgmap_path <- system.file('extdata', 'txp2hgnc.tsv', package = 'drugseqr')

  # run bustools
  tmp_dir <- file.path(out_dir, 'tmp')
  gct_dir <- file.path(out_dir, 'genecount', 'genes')
  dir.create(tmp_dir)
  dir.create(gct_dir, recursive = TRUE)

  system2('bustools',
          args = c('correct',
                   '-w', whitelist_path,
                   '-p', file.path(out_dir, 'output.bus'),
                   '| bustools sort -T', tmp_dir,
                   '-t 4 -p - | bustools count -o', gct_dir,
                   '-g', tgmap_path,
                   '-e', file.path(out_dir, 'matrix.ec'),
                   '-t', file.path(out_dir, 'transcripts.txt'),
                   '--genecounts -'))

  return(NULL)
}

#' Read kallisto/bustools market matrix and annotations
#'
#' Primarily for testing purposes. Recommendation is to use \code{run_alevin} and \code{load_scseq}.
#'
#' @inheritParams run_kallisto
#'
#' @return sparse dgTMatrix with barcodes in columns and genes in rows.
#' @export
#'
#' @examples
#'
read_kallisto <- function(data_dir) {

  count_dir <- file.path(data_dir, 'genecount')

  # read sparse matrix
  counts <- Matrix::readMM(file.path(count_dir, 'genes.mtx'))
  counts <- Matrix::t(counts)

  # read annotations
  row.names(counts) <- read.delim1(file.path(count_dir, 'genes.genes.txt'))
  colnames(counts) <- read.delim1(file.path(count_dir, 'genes.barcodes.txt'))

  return(counts)
}
