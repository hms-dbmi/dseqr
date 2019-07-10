library(drugseqr)

ctrl_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'
test_dir <- 'data-raw/single-cell/example-data/Run2643-10X-Lung/10X_FID12518_Diseased_3hg'

# load raw counts and remove non expressed genes
ctrl_counts <- load_scseq_counts(ctrl_dir)
test_counts <- load_scseq_counts(test_dir)

ctrl_counts <- ctrl_counts[Matrix::rowSums(ctrl_counts) > 0, ]
test_counts <- test_counts[Matrix::rowSums(test_counts) > 0, ]

# calculate empty droplets
ctrl_empty <- get_empty(ctrl_counts)
test_empty <- get_empty(test_counts)

# strain soup
ctrl_strained <- strain_scseq(ctrl_counts, ctrl_empty, 'ctrl')
test_strained <- strain_scseq(test_counts, test_empty, 'test')

# quality control cells
# generate/load whitelist
whitelist <- get_scseq_whitelist(counts, data_dir)
whitelist <- data.frame(whitelist = colnames(counts) %in% whitelist, row.names = colnames(counts))
kneelist  <- readLines(file.path(data_dir, 'bus_output', 'kneelist.txt'))

# get ambient expression profile/determine outlier genes
pct_ambient <- get_pct_ambient(counts)
out_ambient <- get_outliers(pct_ambient)

# covert to Seurat object
srt <- Seurat::CreateSeuratObject(counts[, kneelist], meta.data = whitelist, project = project)

# sctransform norm/variance stabilization
ctrl_scseq <- preprocess_scseq(ctrl_scseq)
test_scseq <- preprocess_scseq(test_scseq)

# add clusters
ctrl_scseq <- add_scseq_clusters(ctrl_scseq, resolution = 1.6)
test_scseq <- add_scseq_clusters(test_scseq)

# run tnse
ctrl_scseq <- run_tsne(ctrl_scseq)
test_scseq <- run_tsne(test_scseq)

# get markers and explore
ctrl_markers <- get_scseq_markers(ctrl_scseq)
test_markers <- get_scseq_markers(test_scseq)

explore_scseq_clusters(ctrl_scseq, ctrl_markers)
explore_scseq_clusters(test_scseq, test_markers)

# integrated analysis ------

# perform integration
anchors <- FindIntegrationAnchors(object.list = list(ctrl_scseq, test_scseq))
combined <- IntegrateData(anchors)
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)

# get clusters/tsne/markers
combined <- add_scseq_clusters(combined)
combined <- run_tsne(combined)
markers <- get_scseq_markers(combined)

# look at how integration changes cell groupings
explore_scseq_clusters(combined, markers)
