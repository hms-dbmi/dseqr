library(drugseqr)

ctrl_dir <- 'data-raw/single-cell/example-data/Run2644-10X-Lung/10X_FID12518_Normal_3hg'
test_dir <- 'data-raw/single-cell/example-data/Run2643-10X-Lung/10X_FID12518_Diseased_3hg'
sjia_dir <- '~/Documents/Batcave/zaklab/drugseqr/data-raw/single-cell/example-anals/sjia'

# load raw counts and remove non expressed genes
ctrl_scseq <- load_scseq(ctrl_dir, project = 'ctrl', soupx = TRUE)
test_scseq <- load_scseq(test_dir, project = 'test', soupx = TRUE)

ctrl_scseq <- ctrl_scseq[, ctrl_scseq$whitelist]
test_scseq <- test_scseq[, test_scseq$whitelist]

ctrl_scseq <- preprocess_scseq(ctrl_scseq)
test_scseq <- preprocess_scseq(test_scseq)

# add clusters and run tsne
ctrl_scseq <- add_scseq_clusters(ctrl_scseq, resolution = 1.4)
test_scseq <- add_scseq_clusters(test_scseq)

ctrl_scseq <- run_umap(ctrl_scseq)
test_scseq <- run_umap(test_scseq)

# get markers
ctrl_markers <- get_scseq_markers(ctrl_scseq)
test_markers <- get_scseq_markers(test_scseq)

ctrl_anal <- list(scseq = ctrl_scseq, markers = ctrl_markers, annot = names(ctrl_markers))
test_anal <- list(scseq = test_scseq, markers = test_markers, annot = names(test_markers))

save_scseq_data(ctrl_anal, 'sjia_lung_healthy_soupx', sjia_dir)
save_scseq_data(test_anal, 'sjia_lung_diseased_soupx', sjia_dir)

explore_scseq_clusters(sjia_dir, test_data = FALSE)
