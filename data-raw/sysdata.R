biogps <- readRDS('data-raw/biogps/biogps.rds')
cell_info <- readRDS('data-raw/cell_info/cell_info.rds')
genes <- readRDS('data-raw/genes/genes.rds')
pert_names <- readRDS('data-raw/drug_gene_queries/pert_names.rds')
ensmap <- readRDS('data-raw/ensmap/ensmap.rds')


azimuth_refs <- c('human_pbmc',
                  'human_lung',
                  'human_lung_v2',
                  'human_bonemarrow',
                  'human_differentiated_tcell',
                  'mouse_til_tcells',
                  'mouse_virus_cd8_tcells',
                  'mouse_virus_cd4_tcells',
                  'human_motorcortex',
                  'mouse_motorcortex',
                  'human_stimulated_pbmc',
                  'human_fetus')

azimuth_labels <- c('PBMC - Human',
                    'Lung V1 - Human',
                    'Lung V2 - Human',
                    'Bone Marrow - Human',
                    'Differentiated CD4 T-cells - Human',
                    'Tumor-Infiltrating T-cells - Mouse',
                    'Virus-Specific CD8 T-cells - Mouse',
                    'Virus-Specific CD4 T-cells - Mouse',
                    'Motor Cortex - Human',
                    'Motor Cortex - Mouse',
                    'PBMC Stimulated - Human',
                    'Fetal Development - Human')

azimuth_species <- ifelse(grepl('human_', azimuth_refs), 'Homo sapiens', 'Mus musculus')

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
gray <- '#FFFFFFCC'
const <- list(
  colors = list(
    bool = c(gray, "#0000FF80"),
    qc = c(gray, 'red'),
    ft = viridis::plasma(10, direction = -1)[-1]
  ),
  features = list(
    qc = c('ribo_percent', 'mito_percent', 'log10_sum', 'log10_detected', 'doublet_score', 'mapping.score'),
    metrics = c('low_lib_size', 'low_n_features', 'high_subsets_mito_percent', 'low_subsets_ribo_percent', 'high_doublet_score'),
    reverse = c('ribo_percent', 'log10_sum', 'log10_detected')
  ),
  max.cells = 80000,
  ref = list(
    azimuth_patterns = paste(
      '^celltype',
      '^annotation',
      '^class$',
      '^cluster$',
      '^subclass$',
      '^cross_species_cluster$',
      '^ann_level_[0-9]$',
      '^ann_finest_level$',
      '^condition$',
      sep = '|'
    ),
    initial_resolns = c(
      'human_pbmc' = 'predicted.celltype.l2',
      'human_lung' = 'predicted.annotation.l1',
      'human_lung_v2' = 'predicted.ann_level_4',
      'human_bonemarrow' = 'predicted.celltype.l2',
      'human_differentiated_tcell' = 'predicted.celltype.cytokines',
      'human_motorcortex' = 'predicted.subclass',
      'human_stimulated_pbmc' = 'predicted.condition',
      'mouse_motorcortex' = 'predicted.subclass',
      'human_fetus' = 'annotation.l1',
      'default' = 'predicted.celltype'
    )
  )
)

usethis::use_data(biogps, cell_info, genes, pert_names, refs, ensmap, const, internal = TRUE, overwrite = TRUE)

