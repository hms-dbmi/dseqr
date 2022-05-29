# samples_merge_list <- list(
#   'SJIALD_20200716' = c('SJIALD_1-20200716_QC2.0', 'SJIALD_2-20200716_QC3'),
#   'SJIALD_FID12518' = c('SJIALD_FID12518_Diseased_QC1.5', 'SJIALD_FID12518_Normal_QC3.5')
# )

# merge samples that are from the sample subject
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

# merge clusters that want to analyze together
merge_clusters <- function(data_dir, merge_list, resoln_dir = NULL) {

  resoln_fpath <- file.path(data_dir, 'resoln.qs')
  if (is.null(resoln_dir)) {
    resoln <- qs::qread(resoln_fpath)
    resoln_dir <- get_resoln_dir(resoln)
  }

  # merge annot/clusters
  resoln_path <- file.path(data_dir, resoln_dir)
  clusters <- qs::qread(file.path(resoln_path, 'clusters.qs'))
  annot <- qs::qread(file.path(resoln_path, 'annot.qs'))

  for (i in seq_along(merge_list)) {
    merge_nums <- merge_list[[i]]
    merge_nums <- as.numeric(merge_nums)
    annot[merge_nums] <- annot[merge_nums[1]]
  }

  levels(clusters) <- annot
  annot <- levels(clusters)
  levels(clusters) <- seq_along(annot)
  annot <- remove.unique(annot)
  annot <- pretty.unique(annot)

  orig_path <- file.path(data_dir, paste0(resoln_dir, '_orig'))
  if (!dir.exists(orig_path)) {
    # store original resolution
    file.rename(resoln_path, orig_path)
    dir.create(resoln_path)

  } else {
    # remove files from previous merge
    prev_files <- list.files(resoln_path, full.names = TRUE, include.dirs = FALSE, recursive = TRUE)
    unlink(prev_files)
  }

  qs::qsave(clusters, file.path(resoln_path, 'clusters.qs'))
  qs::qsave(annot, file.path(resoln_path, 'annot.qs'))
}

