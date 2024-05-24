#' Merge samples that are from the sample subject
#'
#' @param data_dir of single-cell
#' @param merge_list Names are new name, values are current names to merge
#' @param new_founder Name of new founder for dataset for grouping in select menu.
#' Default (\code{NULL}) keeps previous founder.
#'
#' @return called for side effects
#'
#' @examples
#'
#' if (interactive()) {
#'
#'   # FILL THIS IN
#'   sc_dir <- ''
#'
#'   dataset_names <- c('SJIALD_VS_HEALTHY_merged_QC1.0_fastMNN')
#'
#'   merge_list <- list(
#'     'SJIALD_20200716' = c('SJIALD_1-20200716_QC3.0_NO_NEUTRO', 'SJIALD_2-20200716_QC3.0_NO_NEUTRO'),
#'     'SJIALD_FID12518' = c('SJIALD_FID12518_Diseased_QC2.0'))
#'
#'   new_founder <- 'SJIALD_VS_HEALTHY_merged'
#'
#'   for (dataset_name in dataset_names) {
#'     data_dir <- file.path(sc_dir, dataset_name)
#'     merge_samples(data_dir, merge_list, new_founder)
#'   }
#'
#' }
#'
merge_samples <- function(data_dir, merge_list, new_founder = NULL) {
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

  # change founder
  if (!is.null(new_founder)) {
    qs::qsave(new_founder, file.path(data_dir, 'founder.qs'))
  }
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
  if (!dir_exists(orig_path)) {
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

find_merged_clusters <- function(orig_dir, resoln_dir) {

  # check if any of selected clusters differ from original
  orig_clusters <- qs::qread(file.path(orig_dir, 'clusters.qs'))
  curr_clusters <- qs::qread(file.path(resoln_dir, 'clusters.qs'))

  merged_clusters <- c()
  for (cluster in unique(curr_clusters)) {

    # merged if selected cells form more than one cluster in orig
    is.cluster <- curr_clusters == cluster
    nclus.orig <- length(unique(orig_clusters[is.cluster]))

    if (nclus.orig > 1)
      merged_clusters <- c(merged_clusters, cluster)
  }
  return(merged_clusters)
}



# unmerge previously merged clusters
unmerge_clusters <- function(resoln_dir, selected) {

  orig_path <- file.path(paste0(resoln_dir, '_orig'), 'clusters.qs')
  orig_clusters <- qs::qread(orig_path)

  clusters_path <- file.path(resoln_dir, 'clusters.qs')
  clusters <- qs::qread(clusters_path)

  annot_path <- file.path(resoln_dir, 'annot.qs')
  annot <- qs::qread(annot_path)

  next_clus <- length(annot) + 1
  for (cluster in selected) {

    # original values of cluster to undo
    is.cluster <- clusters == cluster
    orig <- orig_clusters[is.cluster]

    # convert to next available integers
    orig.levels <- unique(orig)
    new_next_clus <- next_clus + length(orig.levels)
    new.levels <- as.character(seq(next_clus, new_next_clus-1))
    new <- plyr::mapvalues(orig, from = orig.levels, new.levels)

    # use to overwrite current merged cluster
    levels(clusters) <- c(levels(clusters), new.levels)
    clusters[is.cluster] <- new

    # use annotation of previously merged cluster for restored clusters
    cluster_name <- annot[as.numeric(cluster)]
    annot <- c(annot, rep(cluster_name, length(new.levels)))

    next_clus <- new_next_clus
  }

  # drop unused levels (previously merged)
  clusters <- droplevels(clusters)

  # remove annotation for previously merged clusters
  keep_clusters <- levels(clusters)
  annot <- annot[as.numeric(keep_clusters)]

  # relevel
  annot <- remove.unique(annot)
  annot <- pretty.unique(annot)
  levels(clusters) <- seq_along(annot)

  # save
  qs::qsave(clusters, clusters_path)
  qs::qsave(annot, annot_path)
}

