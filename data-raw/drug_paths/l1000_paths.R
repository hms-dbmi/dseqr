# run pathway analysis for each
library(Biobase)
library(crossmeta)

data_dir <- '/home/alex/Documents/Batcave/zaklab/dseqr/data-raw/drug_paths'
ebfit_dir <- '/run/media/alex/My Passport/Batcave/GEO/l1000/level1/6-limma/ebayes'

ebfits <- list.files(ebfit_dir)

l1000_genes <- colnames(readRDS("~/Documents/Batcave/zaklab/dseqr.data/inst/extdata/l1000_genes_es.rds"))
l1000_drugs <- colnames(readRDS("~/Documents/Batcave/zaklab/dseqr.data/inst/extdata/l1000_drugs_es.rds"))

# get map from l1000_es colname to ebfit filename
l1000_drug_ebfits <- paste0(l1000_drugs, '.rds')
miss <- !l1000_drug_ebfits %in% ebfits
l1000_drug_ebfits[miss] <- gsub('^TUL-X', 'TUL_X', l1000_drug_ebfits[miss])
miss <- !l1000_drug_ebfits %in% ebfits
l1000_drug_ebfits[miss] <- gsub('(+/-)-7-hydroxy-2-(N,N-di-n-propylamino)', 'X......7.hydroxy.2..N.N.di.n.propylamino.', fixed = TRUE, l1000_drug_ebfits[miss])
miss <- !l1000_drug_ebfits %in% ebfits
table(miss)

l1000_gene_ebfits <- paste0(l1000_genes, '.rds')
l1000_gene_ebfits <- gsub('-oe_', '_oe_', l1000_gene_ebfits)
l1000_gene_ebfits <- gsub('-sh_', '_sh_', l1000_gene_ebfits)
l1000_gene_ebfits <- gsub('-lig_', '_lig_', l1000_gene_ebfits)
table(l1000_gene_ebfits %in% ebfits)


# map PROBE --> ENTREZID same as for cmap_paths.R
map <- readRDS(file.path(data_dir, 'map.rds'))
map <- map[, c('PROBE', 'SYMBOL', 'ENTREZID')]
map <- map[!duplicated(map$PROBE), ]

# things that need for each pathway analysis
species <- 'Hs'
gslist.go <- dseqr:::get_gslist(species)
gslist.kegg <- dseqr:::get_gslist(species, type = 'kegg')

gs.names.go <- dseqr:::get_gs.names(gslist.go, species = species)
gs.names.kegg <- dseqr:::get_gs.names(gslist.kegg, type = 'kegg', species = species)


# pathway analysis for each coefficient in ebfit
get_drug_paths <- function(ebfit_path, map, gslist.go, gslist.kegg, gs.names.go, gs.names.kegg, species = 'Hs', min.genes = 4, filter_tt = FALSE) {

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
    goana_res <- dseqr:::add_path_genes(goana_res, gslist.go, ebfit)

    kegga_res <- limma::kegga(ebfit, geneid = 'ENTREZID', species = species)
    kegga_res <- kegga_res[order_or_filt(kegga_res), ]
    kegga_res <- dseqr:::add_path_genes(kegga_res, gslist.kegg, ebfit)

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
  go$FDR <- stats::p.adjust(go$PValue, 'BH')
  go <- go[go$PValue < 0.1, ]

  # get cameraPR KEGG result
  kg <- limma::cameraPR(statistic, index = gslist.kegg)
  kg <- tibble::add_column(kg, Term = gs.names.kegg[row.names(kg)], .before = 'NGenes')
  kg <- kg[kg$NGenes >= min.genes, ]
  kg$FDR <- stats::p.adjust(kg$PValue, 'BH')
  kg <- kg[kg$PValue < 0.1, ]

  if (filter_tt) tt <- tt[tt$P.Value < 0.5, ]

  return(list(tt = tt, go = go, kg = kg, goana = goana_res, kegga = kegga_res))

}

library(foreach)
library(doParallel)
registerDoParallel(10)

dir.create(file.path(data_dir, 'results', 'l1000_drugs'))

res <- foreach(i=1:length(l1000_drug_ebfits)) %dopar% {
  sig_name <- l1000_drugs[i]
  fname <- fs::path_sanitize(paste0(sig_name, '.rds'))
  save_path <- file.path(data_dir, 'results', 'l1000_drugs', fname)
  ebfit_path <- file.path(ebfit_dir, l1000_drug_ebfits[i])

  if (file.exists(save_path)) return(NA)
  path_res <- get_drug_paths(ebfit_path, map, gslist.go, gslist.kegg, gs.names.go, gs.names.kegg)

  saveRDS(path_res, save_path)
}

dir.create(file.path(data_dir, 'results', 'l1000_genes'))

res <- foreach(i=1:length(l1000_gene_ebfits)) %dopar% {
  sig_name <- l1000_genes[i]
  fname <- fs::path_sanitize(paste0(sig_name, '.rds'))
  save_path <- file.path(data_dir, 'results', 'l1000_genes', fname)
  ebfit_path <- file.path(ebfit_dir, l1000_gene_ebfits[i])

  if (file.exists(save_path)) return(NA)
  path_res <- get_drug_paths(ebfit_path, map, gslist.go, gslist.kegg, gs.names.go, gs.names.kegg)

  saveRDS(path_res, save_path)
}
