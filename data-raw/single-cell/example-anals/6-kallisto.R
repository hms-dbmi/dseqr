library(drugseqr)

test_dir <- 'data-raw/single-cell/example-data/Run2643-10X-Lung/10X_FID12518_Diseased_3hg'
whitelist_path <- file.path(test_dir, 'alevin_output_0.14.0/alevin/whitelist.txt')

# load raw counts and remove non expressed genes
counts <- load_scseq_counts(test_dir)

get_knee(counts)
counts <- load_scseq(test_dir, project = 'test')

# subset based on alevin whitelist
whitelist <- readLines(whitelist_path)
# scseq <- counts[, colnames(counts) %in% whitelist]
# scseq <- Seurat::CreateSeuratObject(scseq, project = 'test')

# determine empty droplets
# out <- DropletUtils::emptyDrops(counts)
# is.cell <- which(out$FDR < 0.001)
scseq <- qc_scseq_nb(counts)

sce <- srt_to_sce(scseq)
sce <- add_scseq_qc_metrics(sce)
hist(sce$pct_counts_mito)

# standard analysis pipeline for alevin whitelisted
scseq <- preprocess_scseq(scseq)
scseq <- add_scseq_clusters(scseq)
scseq <- run_tsne(scseq)
markers <- get_scseq_markers(scseq)

# explore results
explore_scseq_clusters(scseq, markers)

