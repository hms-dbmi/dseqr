biogps <- readRDS('data-raw/biogps/biogps.rds')
cell_info <- readRDS('data-raw/cell_info/cell_info.rds')
genes <- readRDS('data-raw/genes/genes.rds')
pert_names <- readRDS('data-raw/drug_gene_queries/pert_names.rds')
ensmap <- readRDS('data-raw/ensmap/ensmap.rds')


azimuth_refs <- c('human_pbmc', 'human_lung', 'human_bonemarrow', 'human_motorcortex', 'mouse_motorcortex')
azimuth_species <- c(rep('Homo sapiens', 4), 'Mus musculus')
azimuth_labels <- c('PBMC - Human',
                    'Lung - Human',
                    'Bone Marrow - Human',
                    'Motor Cortex - Human',
                    'Motor Cortex - Mouse')

symphony_refs <- c('pbmcs_10x', 'scmuscle')
symphony_species <- c('Homo sapiens', 'Mus musculus')
symphony_labels <- c('10x PBMCs Atlas - Human', 'Cornell scMuscle - Mouse')

refs <- data.frame(
    name = c(azimuth_refs, symphony_refs),
    species = c(azimuth_species, symphony_species),
    label = c(azimuth_labels, symphony_labels),
    type = c(rep('Azimuth', length(azimuth_refs)),
             rep('symphony', length(symphony_refs)))
)

# constants
gray <- '#f5f5f5'
const <- list(
    colors = list(
        n0 = gray,
        ft = c(gray, 'blue'),
        qc = c(gray, 'red')
    ),
    features = list(
        qc = c('ribo_percent', 'mito_percent', 'log10_sum', 'log10_detected', 'doublet_score', 'num_significant'),
        metrics = c('low_lib_size', 'low_n_features', 'high_subsets_mito_percent', 'low_subsets_ribo_percent', 'high_doublet_score'),
        reverse = c('ribo_percent', 'log10_sum', 'log10_detected')
    ),
    max.cells = 80000
)

usethis::use_data(biogps, cell_info, genes, pert_names, refs, ensmap, const, internal = TRUE, overwrite = TRUE)

