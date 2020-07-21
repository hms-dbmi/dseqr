# run pathway analysis for each
library(Biobase)
library(crossmeta)

data_dir <- '/home/alex/Documents/Batcave/zaklab/drugseqr/data-raw/drug_paths'
ebfit_dir <- '/home/alex/Documents/Batcave/zaklab/drugseqr/data-raw/drug_paths/tmp'

# make sure names are identical
cmap_es_names <- colnames(readRDS("~/Documents/Batcave/zaklab/drugseqr.data/inst/extdata/cmap_es_ind.rds"))

# split up fit result for pathway analysis -----
# one time setup
rma_processed <- readRDS("~/Documents/Batcave/GEO/ccdata/data-raw/cmap_es/rma_processed_ind.rds")

ebfit <- rma_processed$ebayes_sv

dir.create(ebfit_dir)
for (i in 1:ncol(ebfit)) {
  cat('Working on', i, 'of', ncol(ebfit), '...\n')
  sig_name <- cmap_es_names[i]
  fname <- fs::path_sanitize(paste0(sig_name, '.rds'))
  ebfit_path <- file.path(ebfit_dir, fname)
  save_path <- file.path(data_dir, 'results', 'cmap', fname)
  if (file.exists(ebfit_path)) next()

  saveRDS(ebfit[,i], ebfit_path)
}

# map PROBE --> ENTREZID in ebfit object
eset <- rma_processed$eset
rm(rma_processed); gc()

ensql <- '/home/alex/Documents/Batcave/GEO/crossmeta/data-raw/entrezdt/ensql.sqlite'
annotation(eset) <- 'GPL96'
fData(eset)$PROBE <- featureNames(eset)
sampleNames(eset) <- paste0('s', 1:ncol(eset))

# map <- fData(symbol_annot(eset, ensql = ensql))
# saveRDS(map, file.path(data_dir, 'map.rds'))

# run pathway analyses ------

map <- readRDS(file.path(data_dir, 'map.rds'))
map <- map[, c('PROBE', 'SYMBOL', 'ENTREZID')]
map <- map[!duplicated(map$PROBE), ]


# things that need for each pathway analysis
species <- 'Hs'
gslist.go <- drugseqr:::get_gslist(species)
gslist.kegg <- drugseqr:::get_gslist(species, type = 'kegg')

gs.names.go <- drugseqr:::get_gs.names(gslist.go, species = species)
gs.names.kegg <- drugseqr:::get_gs.names(gslist.kegg, type = 'kegg', species = species)


# pathway analysis for each coefficient in ebfit
# pathway analysis for each coefficient in ebfit
get_drug_paths <- function(ebfit_path, map, gslist.go, gslist.kegg, gs.names.go, gs.names.kegg, species = 'Hs', min.genes = 4, filter_tt = TRUE) {

  # add  ENTREZID to ebfit
  ebfit <- readRDS(ebfit_path)
  map <- map[row.names(ebfit), ]
  ebfit$genes <- map

  order_or_filt <- function(r) {
    pvals <- pmin(r$P.Up, r$P.Down)
    is.sig <- pvals < 0.1
    order(pvals)[is.sig]
  }

  # goana/kegga needs significant genes
  tt <- limma::topTable(ebfit, coef = 1, n = Inf, sort.by = 'p')
  has.sig <- sum(tt$adj.P.Val < 0.05) != 0


  if (has.sig) {
    goana_res <- limma::goana(ebfit, geneid = 'ENTREZID', species = species)
    goana_res <- goana_res[order_or_filt(goana_res), ]
    goana_res <- drugseqr:::add_path_genes(goana_res, gslist.go, ebfit)

    kegga_res <- limma::kegga(ebfit, geneid = 'ENTREZID', species = species)
    kegga_res <- kegga_res[order_or_filt(kegga_res), ]
    kegga_res <- drugseqr:::add_path_genes(kegga_res, gslist.kegg, ebfit)

  } else {
    goana_res <- NULL
    kegga_res <- NULL
  }

  statistic <- ebfit$t[, 1]
  names(statistic) <- ebfit$genes$ENTREZID

  # get cameraPR GO result
  go <- limma::cameraPR(statistic, index = gslist.go)
  go <- tibble::add_column(go, Term = gs.names.go[row.names(go)], .before = 'NGenes')
  go <- go[go$NGenes >= min.genes, ]
  go$FDR <- p.adjust(go$PValue, 'BH')
  go <- go[go$PValue < 0.1, ]

  # get cameraPR KEGG result
  kg <- limma::cameraPR(statistic, index = gslist.kegg)
  kg <- tibble::add_column(kg, Term = gs.names.kegg[row.names(kg)], .before = 'NGenes')
  kg <- kg[kg$NGenes >= min.genes, ]
  kg$FDR <- p.adjust(kg$PValue, 'BH')
  kg <- kg[kg$PValue < 0.1, ]

  if (filter_tt) tt <- tt[tt$P.Value < 0.5, ]

  return(list(tt = tt, go = go, kg = kg, goana = goana_res, kegga = kegga_res))

}

library(foreach)
library(doParallel)
registerDoParallel(4)

dir.create(file.path(data_dir, 'results', 'cmap'))

res <- foreach(i=1:3742) %dopar% {
  sig_name <- cmap_es_names[i]
  fname <- fs::path_sanitize(paste0(sig_name, '.rds'))
  save_path <- file.path(data_dir, 'results', 'cmap', fname)
  ebfit_path <- file.path(ebfit_dir, fname)

  if (file.exists(save_path)) return(NA)
  path_res <- get_drug_paths(ebfit_path, map, gslist.go, gslist.kegg, gs.names.go, gs.names.kegg)

  saveRDS(path_res, save_path)
}
