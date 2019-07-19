library(drugseqr)

cmap_annot <- get_drugs_table('CMAP02')
l1000_annot <- get_drugs_table('L1000')

cmap_cells <- unique(cmap_annot$cell_line)
l1000_cells <- unique(l1000_annot$cell_line)

# GSE92742_Broad_LINCS_cell_info.txt.gz from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742
cell_info <- fread('data-raw/cell_info/GSE92742_Broad_LINCS_cell_info.txt')

cmap_cells[!cmap_cells %in% cell_info$cell_id]
l1000_cells[!l1000_cells %in% cell_info$cell_id]

# fill in info for missing celllines
ssMCF7 <- cell_info[cell_id == 'MCF7']
ssMCF7$cell_id <- 'ssMCF7'
SKMEL5 <- c('SKMEL5', 'cell line', 'SKMEL5', -666, -666, 'tumor', 'skin', 'melanoma', 'adherent', -666, -666, 24, 'F', 'Caucasian')
SKMEL5 <- sapply(SKMEL5, as.list, USE.NAMES = FALSE)
SNUC4 <- c('SNUC4', 'cell line', 'SNUC4', -666, -666, 'tumor', 'large intestine', 'colorectal adenocarcinoma', -666, -666, -666, 35, 'M', -666)
SNUC4 <- sapply(SNUC4, as.list, USE.NAMES = FALSE)

cell_info <- rbind(cell_info, ssMCF7, SKMEL5, SNUC4)

# get columns for app
cell_info <- cell_info[, .(cell_id, primary_site, sample_type)]

# fill in missing site metadata
cell_info[primary_site == -666, cell_id]
setkey(cell_info, cell_id)
cell_info[c('MNEU.E', 'NEU', 'NEU.KCL'), 'primary_site'] <- 'central nervous system'
cell_info['HUES3', 'primary_site'] <- 'embryonic stem cell'

# fill in missing sample type metadata
cell_info[sample_type == -666, cell_id]
cell_info['H1299', 'sample_type'] <- 'tumor'

cell_info[primary_site == 'haematopoietic and lymphoid tissue', 'primary_site'] <- 'haematopoietic/lymphoid'

# needed for selectizeInput
cell_info$value <- cell_info$label <- cell_info$cell_id


# ordering for selectizeInput dropdown
cell_info$primary_site <- factor(cell_info$primary_site,
                                 levels = c(
                                  'haematopoietic/lymphoid',
                                 'blood',
                                 'bone',
                                 'lung',
                                 'large intestine',
                                 'stomach',
                                 'skin',
                                 'muscle',
                                 'adipose',
                                 'kidney',
                                 'liver',
                                 'pancreas',
                                 'central nervous system',
                                 'autonomic ganglia',
                                 'vascular system',
                                 'breast',
                                 'ovary',
                                 'endometrium',
                                 'prostate',
                                 'embryonic stem cell'))

cell_info <- cell_info[order(primary_site, sample_type), .SD, by = primary_site]

cell_info <- list(
  cmap = cell_info[cell_id %in% cmap_cells],
  l1000 = cell_info[cell_id %in% l1000_cells]
)

saveRDS(cell_info, 'data-raw/cell_info/cell_info.rds')
