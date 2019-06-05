library(drugseqr)

data_dir <- 'data-raw/single-cell/example-data/Run2643-10X-Lung/10X_FID12518_Diseased_3hg'
scseq.full <- load_scseq(data_dir)

# qc plots
hist_scseq_whitelist(scseq.full)

# run sctransform
scseq.full <- preprocess_scseq(scseq.full)
tsne_scseq_whitelist(scseq.full)

# subset by alevin whitelist
scseq <- scseq.full[, scseq.full$whitelist]
scseq <- preprocess_scseq(scseq)

# add clusters
scseq <- add_scseq_clusters(scseq)

# run tnse and explore clusters
scseq <- run_tsne(scseq)
markers <- get_scseq_markers(scseq)
explore_scseq_clusters(scseq, markers)

