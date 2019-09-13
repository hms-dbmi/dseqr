#' Download ensembl transcriptome and build index for salmon quantification
#'
#' This index is used for bulk RNA seq quantification. See \code{\link{build_gencode_index}} for single cell RNA-seq equivalent.
#'
#' @param species The species. Default is \code{homo_sapiens.}
#' @param release ensembl release. Default is \code{94} (latest in release for AnnotationHub -
#'   needs to match with \code{\link{build_ensdb}}) and corresponds to Gencode release 29 for \code{\link{build_gencode_index}}.
#' @param indices_dir Directory with salmon indices.
#' @param command System command for salmon.
#'
#' @return NULL
#' @export
#'
#' @examples
#' # build salmon index for humans
#' build_ensdb_index()
#'
build_ensdb_index <- function(indices_dir, species = 'homo_sapiens', release = '94', command = 'salmon') {

  salmon_version <- get_salmon_version(command)
  if (salmon_version == '0.14.0')
    stop('EnsemblDb indices not yet implemented for salmon 0.14.0')

  indices_dir <- file.path(indices_dir, salmon_version, 'ensdb')
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
  ensembl_all <- grep('cdna.all.fa.gz$', tbl[[1]], value = TRUE)
  ensembl_url <- paste0(ensembl_url, ensembl_all)

  work_dir <- getwd()
  setwd(indices_dir)
  curl::curl_download(ensembl_url, ensembl_all)

  # build index
  tryCatch(system2(command, args=c('index',
                                   '-t', ensembl_all,
                                   '-i', ensembl_species)),
           error = function(err) {err$message <- 'Is salmon installed and on the PATH?'; stop(err)})

  unlink(ensembl_all)
  setwd(work_dir)
}

#' Download ensembl transcriptome and build index for salmon quantification
#'
#' This index is used for single cell RNA-seq quantification. See \code{\link{build_ensdb_index}} for bulk RNA-seq equivalent.
#'
#' @param indices_dir directory to place indices in.
#' @param species The species. Default is \code{homo_sapiens.}
#' @param release gencode release. Default is \code{29} (matches ensembl release 94 for \code{\link{build_ensdb_index}})).
#'
#' @return NULL
#' @export
#'
#' @examples
#' # build salmon alevin index for humans
#' indices_dir <- 'data-raw/indices'
#' build_gencode_index(indices_dir)
#'
build_gencode_index <- function(indices_dir, species = 'human', release = '29', command = 'salmon') {

  salmon_version <- get_salmon_version(command)
  salmon_old <- salmon_lt_0.14.0(salmon_version)

  # newer salmon uses decoys in index
  if (salmon_old)
    build_gencode_index_no_decoys(indices_dir, species, release)
  else
    build_gencode_index_decoys(indices_dir, species, release)
}

build_gencode_index_decoys <- function(indices_dir, species, release) {

  indices_dir <- file.path(indices_dir, '0.14.0', 'gencode')
  if (!dir.exists(indices_dir)) dir.create(indices_dir, recursive = TRUE)

  if (species != 'human' | release != '29')
    stop('Only implemented for human gencode v29.')

  gencode_file <- 'human_GENCODEv29.tar.gz'
  gencode_url <- paste0('https://s3.us-east-2.amazonaws.com/drugseqr/', gencode_file)

  work_dir <- getwd()
  setwd(indices_dir)

  download.file(gencode_url, gencode_file)
  utils::untar(gencode_file)

  # build index
  tryCatch(system2('salmon', args=c('index',
                                    '-t', 'human_GENCODEv29/gentrome.fa',
                                    '--gencode',
                                    '-d', 'human_GENCODEv29/decoys.txt',
                                    '-i', species)),
           error = function(err) {err$message <- 'Is salmon installed and on the PATH?'; stop(err)})


  unlink(c(gencode_file, 'human_GENCODEv29'), recursive = TRUE)
  setwd(work_dir)
}


build_gencode_index_no_decoys <- function(indices_dir, species, release) {

  indices_dir <- file.path(indices_dir, '0.13.1', 'gencode')
  if (!dir.exists(indices_dir)) dir.create(indices_dir, recursive = TRUE)

  # construct ensembl url for protein coding transcriptome
  gencode_file <- paste0('gencode.v', release, '.pc_transcripts.fa.gz')
  gencode_url <- paste0('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_', species, '/release_', release, '/', gencode_file)

  work_dir <- getwd()
  setwd(indices_dir)
  curl::curl_download(gencode_url, gencode_file)

  # build index
  tryCatch(system2('salmon', args=c('index',
                                    '-t', gencode_file,
                                    '--gencode',
                                    '-i', species)),
           error = function(err) {err$message <- 'Is salmon installed and on the PATH?'; stop(err)})

  unlink(gencode_file)
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
#' run_salmon_bulk(indices_dir, data_dir, pdata_path)
#'
run_salmon_bulk <- function(indices_dir, data_dir, pdata_path = NULL, pdata = NULL, species = 'homo_sapiens', command = 'salmon',
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

#' Get version of salmon from system command.
#'
#' @param command System command for salmon. Used to determine salmon version.
#'
#' @return Version of salmon.
#' @export
#' @keywords internal
#'
#' @examples
get_salmon_version <- function(command) {
  # possibly use older salmon with version appended to executable name
  salmon_version <- system(paste(command, '--version'), intern = TRUE)
  salmon_version <- gsub('^salmon ', '', salmon_version)
  return(salmon_version)
}
