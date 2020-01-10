# example script showing how to process multiple single-cell datasets without interaction

sc_dir <- '/home/alex/Documents/Batcave/zaklab/drugseqr/data-raw/patient_data/sjia/single-cell'

new_dir <- file.path(sc_dir, '10X-data')
status_dirs <- list.files(new_dir, full.names = TRUE)
names(status_dirs) <- c('active', 'healthy', 'inactive', 'lung_disease')

sample_dirs_list <- lapply(status_dirs, list.files)
indices_dir <- '/srv/drugseqr/indices'

# through each status
for (i in seq_along(status_dirs)) {
  cat('working on status', i, 'of', length(status_dirs), '...\n')
  status_dir <- status_dirs[i]
  status <- names(status_dir)

  # through each sample
  sample_dirs <- sample_dirs_list[[i]]

  for (j in seq_along(sample_dirs)) {
    cat('working on sample', j, 'of', length(sample_dirs), '\n')

    fastq_dir <- file.path(status_dir, sample_dirs[j], 'fastqs')
    dataset_name <- paste0('sjia_pbmcs_', status, j)

    if (dir.exists(file.path(sc_dir, dataset_name))) next()

    run_kallisto_scseq(indices_dir, fastq_dir)

    cat('loading...\n')
    scseq <- load_scseq(fastq_dir, project = dataset_name, type = 'kallisto')
    scseq <- scseq[, scseq$whitelist]
    gc()

    cat('normalizing...\n')
    scseq <- normalize_scseq(scseq)
    gc()

    cat('clustering...\n')
    scseq <- add_hvgs(scseq)
    scseq <- add_scseq_clusters(scseq)
    gc()

    cat('reducing dims...\n')
    scseq <- run_tsne(scseq)
    gc()

    cat('getting markers...\n')
    wilcox_tests <- pairwise_wilcox(scseq)
    markers <- get_scseq_markers(wilcox_tests)

    # top markers for SingleR
    top_markers <- scran::getTopMarkers(wilcox_tests$statistics, wilcox_tests$pairs)

    cat('saving...\n')
    anal <- list(scseq = scseq, markers = markers, annot = names(markers), top_markers = top_markers)
    save_scseq_data(anal, dataset_name, sc_dir)
  }
}
