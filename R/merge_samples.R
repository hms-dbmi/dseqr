merge_samples <- function(data_dir, merge_list) {
  # open shell and merge batch and ambience
  shell <- qs::qread(file.path(data_dir, 'shell.qs'))
  rdata <- SummarizedExperiment::rowData(shell)

  for (merge_name in names(merge_list)) {
    merge_samples <- merge_list[[merge_name]]
    shell$batch[shell$batch %in% merge_samples] <- merge_name

    # TODO: allow for samples without ambience (was pre-filtered)
    amb.merge <- paste0(merge_samples, '_ambience')
    rdata[[paste0(merge_name, '_ambience')]] <- rowSums(as.matrix(rdata[, amb.merge]))
    rdata[amb.merge] <- NULL
  }
  SummarizedExperiment::rowData(shell) <- rdata
  qs::qsave(shell, file.path(data_dir, 'shell.qs'))

  # remove base files related to differential analyses
  unlink(list.files(data_dir, '^lm_fit_grid_|^meta.qs$|^pairs.qs$|^prev_groups.qs$|^summed_grid.qs$', full.names = TRUE))

  # remove saved plots
  unlink(file.path(data_dir, 'plots'), recursive = TRUE)

  # remove analysis sub-directories
  dirs <- list.dirs(data_dir, full.names = FALSE, recursive = TRUE)
  subdirs <- grep('/', dirs, value = TRUE)
  unlink(file.path(data_dir, subdirs), recursive = TRUE)

  # remove resolution files related to differential analyses
  resoln_dirs <- setdiff(dirs, c(subdirs, ''))
  resoln_files <- list.files(file.path(data_dir, resoln_dirs), full.names = TRUE)
  resoln_basenames <- basename(resoln_files)
  keep <- resoln_basenames %in% c('annot.qs', 'clusters.qs')
  unlink(resoln_files[!keep])



}

merge_list <- list(
  'SJIALD_20200716' = c('SJIALD_1-20200716_QC1.5', 'SJIALD_2-20200716_QC2'),
  'SJIALD_FID12518' = c('SJIALD_FID12518_Diseased_QC0.5', 'SJIALD_FID12518_Normal_QC0.5')
)
