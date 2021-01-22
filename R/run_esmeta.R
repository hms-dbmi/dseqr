run_esmeta <- function(tts) {
  # only consider genes
  anals <- lapply(tts, function(tt) list(top_tables = list(tt)))
  es <- crossmeta::es_meta(anals)
  es <- es$all$filt
  return(es)
}

extract_enids <- function(tts) {
  # get enids
  syms <- lapply(tts, row.names)
  enids <- lapply(tts, `[[`, 'ENTREZID')
  enids <- unlist(enids)
  names(enids) <- unlist(syms)
  enids[!duplicated(enids)]
}


es_to_tt <- function(es, enids, cols) {
  keep <- c('mu', 'var', 'z', 'pval', 'fdr')
  new <- c('dprime', 'vardprime', 't', 'P.Value', 'adj.P.Val')

  es <- es[, keep]
  colnames(es) <- new
  es <- es[order(es$P.Value), ]

  # make look like top table for app consistency
  es$AveExpr <- es$B <- NA
  es$logFC <- es$dprime
  es$ENTREZID <- enids[row.names(es)]
  es$adj.P.Val.Amb <- es$adj.P.Val
  es$ambient <- FALSE

  return(es[, cols])
}

tt_to_es <- function(tt) {
  keep <- c('ENTREZID', 'dprime', 'vardprime', 't', 'P.Value', 'adj.P.Val')
  new <- c('ENTREZID', 'mu', 'var', 'z', 'pval', 'fdr')
  es <- tt[, keep]
  colnames(es) <- new
  return(es)
}

write.csv.safe <- function(df, fname, tozip) {
  if (!nrow(df)) return(tozip)

  utils::write.csv(df, fname)
  tozip <- c(tozip, fname)
  return(tozip)
}
