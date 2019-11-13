#' Differential expression of KEGG pathways.
#'
#' Performs PADOG pathway analysis using KEGG database.
#' PADOG outperforms other pathway analysis algorithms at prioritizing expected pathways (see references).
#'
#' @param eset Expression set.
#' @param prev_anal Previous result of \code{\link{diff_expr}}.
#'
#' @references Tarca AL, Bhatti G, Romero R. A Comparison of Gene Set Analysis Methods
#'    in Terms of Sensitivity, Prioritization and Specificity. Chen L, ed. PLoS ONE.
#'    2013;8(11):e79217. doi:10.1371/journal.pone.0079217.
#'
#'    Dong X, Hao Y, Wang X, Tian W. LEGO: a novel method for gene set over-representation
#'    analysis by incorporating network-based gene weights. Scientific Reports.
#'    2016;6:18871. doi:10.1038/srep18871.
#'
#' @export
#'
#' @examples
#' eset <- readRDS('data-raw/patient_data/sjia/bulk/mono/eset.rds')
#' prev_anal <- readRDS('data-raw/patient_data/sjia/bulk/mono/diff_expr_symbol_ferritin.rds')
#' data_dir <- 'data-raw/patient_data/sjia/bulk/mono'
#' anal_name <- 'ferritin'
#'
#' path_anal <- diff_path(eset, prev_anal, data_dir, anal_name, NI = 24)
#'
diff_path <- function(eset, prev_anal, data_dir, anal_name, rna_seq = TRUE, browse = FALSE, NI = 1000, type = c('KEGG', 'GO')){

  # subset eset to previously analysed samples
  prev_pdata <- prev_anal$pdata
  eset <- eset[, row.names(prev_pdata)]
  Biobase::pData(eset)$group <- prev_pdata$group

  if (rna_seq) eset <- add_vsd(eset)

  # remove replicates/duplicates and annotate with human ENTREZID
  dups <- iqr_replicates(eset, annot = 'ENTREZID_HS', rm.dup = TRUE)
  eset <- dups$eset

  # contrast levels
  groups <- prev_pdata$group
  ctrl <- 'ctrl'
  test <- 'test'

  # group
  incon <- groups %in% c(ctrl, test)
  group <- ifelse(groups[incon] == ctrl, 'c', 'd')

  # other padog inputs
  esetm  <- Biobase::exprs(eset[, incon])

  if (type[1] == 'KEGG') {
    gs.names <- gs.names.kegg
    gslist <- gslist.kegg
  } else if (type[1] == 'GO') {
    gs.names <- gs.names.go
    gslist <- gslist.go
  }

  # run padog
  padog_table <- PADOG::padog(esetm = esetm, group = group, parallel = TRUE, ncr = 4, gs.names = gs.names, gslist = gslist,
                              verbose = FALSE, rna_seq = rna_seq, pdata = prev_pdata, browse = browse, NI = NI)

  # save results
  fname <- paste0('diff_path_', tolower(type[1]), '_', anal_name, '.rds')
  saveRDS(padog_table, file.path(data_dir, fname))

  return(padog_table)
}
