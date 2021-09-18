biogps <- readRDS('data-raw/biogps/biogps.rds')
cell_info <- readRDS('data-raw/cell_info/cell_info.rds')
genes <- readRDS('data-raw/genes/genes.rds')
pert_names <- readRDS('data-raw/drug_gene_queries/pert_names.rds')
ensmap <- readRDS('data-raw/ensmap/ensmap.rds')

azimuth_refs <- c('human_pbmc', 'human_lung', 'human_motorcortex', 'mouse_motorcortex')
names(azimuth_refs) <- c(rep('Homo sapiens', 3), 'Mus musculus')

# constants
gray <- '#f5f5f5'
const <- list(
    colors = list(
        n0 = gray,
        ft = c(gray, 'blue'),
        qc = c(gray, 'red')
    ),
    features = list(
        qc = c('ribo_percent', 'mito_percent', 'log10_sum', 'log10_detected', 'doublet_score'),
        reverse = c('ribo_percent', 'log10_sum', 'log10_detected')
    )
)

usethis::use_data(biogps, cell_info, genes, pert_names, azimuth_refs, ensmap, const, internal = TRUE, overwrite = TRUE)
