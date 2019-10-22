library(Seurat)
# create small test-datasets for testing scseq shiny app

# where the original analyses are
source_dir <- 'data-raw/patient_data/example/'

# where we will put a few of them
static_dir <- 'inst/app/tests/data/static/'
unlink(static_dir, recursive = TRUE)
dir.create(static_dir)

# copy over the original analyses

file.copy(source_dir, static_dir, recursive = TRUE)

scseq_paths <- list.files(static_dir, 'scseq.rds', recursive = TRUE, full.names = TRUE)

# make smaller versions
for (scseq_path in scseq_paths) {

  # keep at most 20 cells in each cluster
  scseq <- readRDS(scseq_path)
  clusters <- scseq$seurat_clusters
  keep <- c()
  for (cluster in levels(clusters)) {
    in.cluster <- colnames(scseq)[clusters == cluster]
    keep_cluster <- sample(in.cluster, min(length(in.cluster), 20))
    keep <- c(keep, keep_cluster)
  }

  scseq <- scseq[, keep]

  saveRDS(scseq, scseq_path)
}
