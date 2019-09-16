
#' Get metadata needed to download RNA-seq data for GSE
#'
#' Goes to GSE page to get GSMs then goes to each GSM page to get SRX then to each SRX page to get some more metadata
#' Adapted from Alex's personal SRAdb package.
#'
#' @param gse_name GEO study name to get metadata for
#'
#' @return \code{data.frame} with sample annotations for each GSE
#' @export
#'
#' @examples
#'
#' srp_meta <- get_srp_meta('GSE117570')
get_srp_meta <- function(gse_name, data_dir = getwd()) {

  gse_dir <- file.path(data_dir, gse_name)
  if (!dir.exists(gse_dir)) dir.create(gse_dir)

  # load srp meta from file if exists
  srp_meta_path <- file.path(gse_dir, 'srp_meta.rds')
  if (file.exists(srp_meta_path)) return(readRDS(srp_meta_path))

  # get GSM names ----

  # get html text for GSE page
  gse_url  <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gse_name)

  gse_html <- NULL
  attempt <- 1
  while(is.null(gse_html) && attempt <= 3) {
    try(gse_html <- xml2::read_html(gse_url))
    if(is.null(gse_html)) Sys.sleep(15)
  }
  gse_text <- rvest::html_text(gse_html)

  # GSM names
  samples <- stringr::str_extract(gse_text, stringr::regex('\nSamples.+?\nRelations', dotall = TRUE))
  samples <- stringr::str_extract_all(samples, 'GSM\\d+')[[1]]

  # get SRX for each GSM ----

  # save in srp_meta
  srp_meta <- data.frame(stringsAsFactors = FALSE)

  cat(length(samples), 'GSMs to process\n')
  for (j in 1:length(samples)) {
    cat('Working on GSM number', j, '\n')

    gsm_name <- samples[j]
    # get html text
    gsm_url  <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gsm_name)

    gsm_html <- NULL
    attempt <- 1
    while(is.null(gsm_html) && attempt <= 3) {
      try(gsm_html <- xml2::read_html(gsm_url))
      if(is.null(gsm_html)) Sys.sleep(5)
    }
    gsm_text <- rvest::html_text(gsm_html)

    # get SRA number for this GSM
    experiment <- stringr::str_extract(gsm_text, 'SRA\nSRX\\d+\n')
    experiment <- gsub('SRA\n(.+?)\n', '\\1', experiment)

    # extract other GSM data while we are here
    title <- stringr::str_extract(gsm_text, '\nTitle\n.+?\n')
    sample_type <- stringr::str_extract(gsm_text, '\nSample type\n.+?\n')
    source_name <- stringr::str_extract(gsm_text, '\nSource name\n.+?\n')
    organism <- stringr::str_extract(gsm_text, '\nOrganism\n.+?\n')
    characteristics <- stringr::str_extract(gsm_text, '\nCharacteristics\n.+?\n')
    treatment_protocol <- stringr::str_extract(gsm_text, '\nTreatment protocol\n.+?\n')
    growth_protocol <- stringr::str_extract(gsm_text, '\nGrowth protocol\n.+?\n')
    extracted_molecule <- stringr::str_extract(gsm_text, '\nExtracted molecule\n.+?\n')
    extraction_protocol <- stringr::str_extract(gsm_text, stringr::regex("\nExtraction protocol\n.+?\nLibrary strategy", dotall = TRUE))
    extraction_protocol <- gsub('\nExtraction protocol\n(.+?)\n\nLibrary strategy', '\\1', extraction_protocol)
    library_strategy <- stringr::str_extract(gsm_text, '\nLibrary strategy\n.+?\n')
    library_selection <- stringr::str_extract(gsm_text, '\nLibrary selection\n.+?\n')
    instrument_model <- stringr::str_extract(gsm_text, '\nInstrument model\n.+?\n')
    description <- stringr::str_extract(gsm_text, '\nDescription\n.+?\n')
    data_processing <- stringr::str_extract(gsm_text, stringr::regex('\nData processing\n.+?\nSubmission date\n', dotall = TRUE))
    data_processing <- gsub('\nData processing\n(.+?)\n\nSubmission date\n', '\\1', data_processing)
    platform_id <- stringr::str_extract(gsm_text, '\nPlatform ID\n.+?\n')

    if (!is.na(experiment)) {

      # extract SRR runs info -----
      srx_url  <- paste0("https://www.ncbi.nlm.nih.gov/sra/", experiment, '[accn]?report=FullXml')

      srx_html <- NULL
      attempt <- 1
      while(is.null(srx_html) && attempt <= 3) {
        try(srx_html <- xml2::read_html(srx_url))
        if(is.null(srx_html)) Sys.sleep(5)
      }
      srx_text <- rvest::html_text(srx_html)

      runs <- stringr::str_extract_all(srx_text, '<PRIMARY_ID>SRR\\d+</PRIMARY_ID>')[[1]]
      runs <- gsub('<PRIMARY_ID>(SRR\\d+)</PRIMARY_ID>', '\\1', unique(runs))
      taxon_id <- stringr::str_extract(srx_text, '<TAXON_ID>\\d+</TAXON_ID>')
      library_source <- stringr::str_extract(srx_text, '<LIBRARY_SOURCE>.+?</LIBRARY_SOURCE>')
      library_layout <- stringr::str_extract(srx_text, '<LIBRARY_LAYOUT>.+?</LIBRARY_LAYOUT>')
      library_layout <- stringr::str_extract(library_layout, 'SINGLE|PAIRED')

      if (length(runs)) {
        # add info to srp_meta
        srp_meta[runs, 'run'] <- runs
        srp_meta[runs, 'experiment'] <- experiment
        srp_meta[runs, 'gsm_name'] <- gsm_name
        srp_meta[runs, 'title'] <- gsub('\nTitle\n(.+?)\n', '\\1', title)
        srp_meta[runs, 'sample_type'] <- gsub('\nSample type\n(.+?)\n', '\\1', sample_type)
        srp_meta[runs, 'source_name'] <- gsub('\nSource name\n(.+?)\n', '\\1', source_name)
        srp_meta[runs, 'organism'] <- gsub('\nOrganism\n(.+?)\n', '\\1', organism)
        srp_meta[runs, 'characteristics'] <- gsub('\nCharacteristics\n(.+?)\n', '\\1', characteristics)
        srp_meta[runs, 'treatment_protocol'] <- gsub('\nTreatment protocol\n(.+?)\n', '\\1', treatment_protocol)
        srp_meta[runs, 'growth_protocol'] <- gsub('\nGrowth protocol\n(.+?)\n', '\\1', growth_protocol)
        srp_meta[runs, 'extracted_molecule'] <- gsub('\nExtracted molecule\n(.+?)\n', '\\1', extracted_molecule)
        srp_meta[runs, 'extraction_protocol'] <- substr(extraction_protocol, 1, nchar(extraction_protocol)-2)
        srp_meta[runs, 'library_strategy'] <- gsub('\nLibrary strategy\n(.+?)\n', '\\1', library_strategy)
        srp_meta[runs, 'library_selection'] <- gsub('\nLibrary selection\n(.+?)\n', '\\1', library_selection)
        srp_meta[runs, 'instrument_model'] <- gsub('\nInstrument model\n(.+?)\n', '\\1', instrument_model)
        srp_meta[runs, 'description'] <- gsub('\nDescription\n(.+?)\n', '\\1', description)
        srp_meta[runs, 'data_processing'] <- substr(data_processing, 1, nchar(data_processing)-2)
        srp_meta[runs, 'platform_id'] <- gsub('\nPlatform ID\n(.+?)\n', '\\1', platform_id)

        srp_meta[runs, 'library_source'] <- gsub('<LIBRARY_SOURCE>(.+?)</LIBRARY_SOURCE>', '\\1', library_source)
        srp_meta[runs, 'library_layout'] <- library_layout
        srp_meta[runs, 'taxon_id'] <- gsub('<TAXON_ID>(\\d+)</TAXON_ID>', '\\1', taxon_id)
        srp_meta[runs, 'ebi_dir']  <- sapply(runs, get_dldir, 'ebi')
        srp_meta[runs, 'ncbi_dir'] <- sapply(runs, get_dldir, 'ncbi')
      }
    }
  }

  saveRDS(srp_meta, srp_meta_path)
  return(srp_meta)
}


#' Gets dldir for add_srp_contacts
#'
#' @param srr
#' @param type
#' @param prefix
#'
#' @return
#' @export
#'
#' @examples
get_dldir <- function(srr, type = 'ebi', prefix = '^SRR') {


  dir1 <- substr(srr, 1, 6)

  if (type == 'ebi') {
    digits  <- gsub(prefix, '', srr)
    ndigits <- nchar(digits)

    if (ndigits == 7) {
      dir2 <- paste0('00', substr(digits, 7, 7))

    } else if (ndigits == 8) {
      dir2 <- paste0('0', substr(digits, 7, 8))

    } else if (ndigits == 9) {
      dir2 <- substr(digits, 7, 9)

    } else if (ndigits == 6) {
      return(file.path(dir1, srr))
    }

    return(file.path(dir1, dir2, srr))

  } else if (type == 'ncbi') {
    return(file.path(dir1, srr))
  }

}


#' Download and fasterq-dump RNA-seq data from GEO
#'
#' First tries to get RNA-Seq fastq files from EBI. If not successful, gets SRA from GEO and converts to fastq.gz.
#'
#' @param gse_name GSE name. Will create folder with this name in \code{data_dir} and download data there.
#' @param srp_meta \code{data.frame} with SRP meta info. Returned from \code{\link{get_srp_meta}}.
#' @param data_dir Path to folder that \code{gse_name} folder will be created in.
#'
#' @export
get_fastqs <- function(gse_name, srp_meta, data_dir) {

  # setup gse directory
  gse_dir  <- file.path(data_dir, gse_name)
  dir.create(gse_dir)


  # seperate runs based on GSM (can be multiple per GSM)
  srr_names <- srp_meta$run
  gsm_names <- unique(srp_meta$gsm_name)
  srr_names_list <- lapply(gsm_names, function(gsm_name) srr_names[srp_meta$gsm_name %in% gsm_name])
  names(srr_names_list) <- gsm_names

  # download everything
  for (i in seq_along(srr_names_list)) {
    gsm_name <- names(srr_names_list)[i]
    gsm_dir <- file.path(gse_dir, gsm_name)
    dir.create(gsm_dir)

    srr_names <- srr_names_list[[1]]

    for (srr_name in srr_names) {
      # try to get fastq from ebi
      get_ebi_fastqs(srp_meta, srr_name, gsm_dir)

      # download/fasterq-dump SRA files
      # run_fasterq_dump(srr_name, out_dir = file.path(gse_dir, srr_name))
    }
  }
  return(NULL)
}

#' Download fastqs from EBI
#'
#' Much faster to use aspera than ftp
#'
#' @param srp_meta Result from \code{get_srp_meta}.
#' @param srr_name Run accession as string.
#' @param gsm_dir Folder to save fastq.gz files in
#' @param method One of either \code{'aspera'} (Default) or \code{'ftp'}.
#'
#' @export
#'
get_ebi_fastqs <- function(srp_meta, srr_name, gsm_dir, method = c('aspera', 'ftp')) {
  url <- paste0('ftp://ftp.sra.ebi.ac.uk/vol1/fastq/', srp_meta[srr_name, 'ebi_dir'], '/.')
  fnames <- unlist(strsplit(getURL(url, dirlistonly = TRUE), '\n'))

  if (method[1] == 'aspera') {
    ascp_path <- system('which ascp', intern = TRUE)
    ascp_pubkey <- gsub('bin/ascp$', 'etc/asperaweb_id_dsa.openssh', ascp_path)

    ascpCMD <- paste('ascp -QT -l 1g -P33001 -i', ascp_pubkey)
    files <- paste0('era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/', srp_meta[srr_name, 'ebi_dir'], '/', fnames)
    ascpR(ascpCMD, files, gsm_dir)

  } else if (method[1] == 'ftp') {
    files <- paste0('ftp://ftp.sra.ebi.ac.uk/vol1/fastq/', srp_meta[srr_name, 'ebi_dir'], '/', fnames)
    download.file(files, file.path(gsm_dir, fnames))
  }
}

#' Utility function to run aspera
#'
#' @param ascpCMD aspera command as string.
#' @param files Urls to aspera files to download.
#' @param destDir Path to directory to download \code{files} into.
#'
#' @export
ascpR <- function (ascpCMD, files, destDir = getwd()) {
  ret <- c()
  for (src in files) {
    ascp_cmd <- paste(ascpCMD, src, destDir, sep = " ")
    ret <- c(ret, system(ascp_cmd))
  }
  return(ret)
}


#' Run fasterq-dump from ncbi/sra-tools
#'
#' Downloads SRR.sra file and converts to fastq files.
#'
#' @param srr_name Run name starting with \code{'SRR'}
#' @param out_dir Directory to dump fastq files into.
#' @param args Character vector of additional arguments to fasterq-dump
#'
#' @return
#' @export
#' @keywords internal
#'
#' @examples
run_fasterq_dump <- function(srr_name, out_dir, args = NULL) {

  # update fastq result if fastq-dump works
  system2('fasterq-dump',
          args = c(
            srr_name,
            '--outdir', out_dir,
            args))
}
