

# load counts and soup
data_dir <- 'data-raw/single-cell/example-data/Run2643-10X-Lung/10X_FID12518_Diseased_3hg'
counts <- load_scseq_counts(data_dir)
soup <- load_scseq_counts(data_dir, soup = TRUE)


# only keep soup with 1-10 counts
soup <- soup[Matrix::rowSums(soup) <= 10, ]

# concatenate
empty <- c(rep(FALSE, nrow(counts)), rep(TRUE, nrow(soup)))
counts <- Matrix::t(rbind(counts, soup))

# strain
scseq <- strain_scseq(counts, empty, 'test')
whitelist <- readLines('data-raw/single-cell/example-data/Run2643-10X-Lung/10X_FID12518_Diseased_3hg/alevin_output_0.14.0/alevin/whitelist.txt')
scseq <- scseq[, whitelist]

# sctransform norm/variance stabilization
scseq <- preprocess_scseq(scseq)

# add clusters
scseq <- add_scseq_clusters(scseq)

# run tnse
scseq <- run_tsne(scseq)

# get markers and explore
markers <- get_scseq_markers(scseq)

explore_scseq_clusters(scseq, markers)
