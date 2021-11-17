biogps <- readRDS('data-raw/biogps/biogps.rds')
cell_info <- readRDS('data-raw/cell_info/cell_info.rds')
genes <- readRDS('data-raw/genes/genes.rds')
pert_names <- readRDS('data-raw/drug_gene_queries/pert_names.rds')
ensmap <- readRDS('data-raw/ensmap/ensmap.rds')


azimuth_refs <- c('human_pbmc', 'human_lung', 'human_motorcortex', 'mouse_motorcortex')
azimuth_species <- c(rep('Homo sapiens', 3), 'Mus musculus')
azimuth_labels <- c('Human - PBMC', 'Human - Lung', 'Human - Motor Cortex', 'Mouse - Motor Cortex')

symphony_refs <- c('pbmcs_10x')
symphony_species <- rep('Homo sapiens', 1)
symphony_labels <- c('10x PBMCs Atlas')

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
        qc = c('ribo_percent', 'mito_percent', 'log10_sum', 'log10_detected', 'doublet_score'),
        reverse = c('ribo_percent', 'log10_sum', 'log10_detected')
    )
)

usethis::use_data(biogps, cell_info, genes, pert_names, refs, ensmap, const, internal = TRUE, overwrite = TRUE)
