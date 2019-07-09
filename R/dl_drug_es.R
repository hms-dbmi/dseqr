#' Download drug effect size data.
#'
#' @param files Character vector of file names to download.
#' @param check Check that existing drug effect size data is loadable? Default is FALSE.
#'
#' @return NULL
#' @export
dl_drug_es <- function(files = c('cmap_es_ind.rds', 'l1000_es.rds'), check = FALSE) {

  # make sure doesn't already exist
  dest_dir <- system.file('extdata', package = 'drugseqr')
  if (!dir.exists(dest_dir)) stop('Make sure you build drugseqr first.')


  can_load <- c()
  exist_files <- file.exists(file.path(dest_dir, files))
  exist_files <- files[exist_files]
  if (length(exist_files) & check) {
    message(paste(exist_files, collapse = ' and '), ' already exists.')

    # check that can load
    for (exist_file in exist_files) {
      message('Checking that ', exist_file, ' can be loaded.')

      # store file name if can load
      fname <- tryCatch({
        drug_es <- readRDS(file.path(dest_dir, exist_file))
        exist_file
      },
      error = function(err) {
        message("Couldn't load ", exist_file, '. Will download.')
        unlink(file.path(dest_dir, exist_file))
        return(NULL)
      })
      can_load <- c(can_load, fname)
    }
    exist_files <- can_load
  }

  need_files <- setdiff(files, can_load)
  if (!length(need_files)) return(NULL)

  for (need_file in need_files) {
    message('downloading: ', need_file)
    dl_url <- paste0('https://s3.us-east-2.amazonaws.com/drugseqr/', need_file)
    download.file(dl_url, file.path(dest_dir, need_file))
  }
}
