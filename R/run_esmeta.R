run_esmeta <- function(tts) {
  # only for clusters with residual dof
  anals <- list()
  for (clust in names(tts)) {
    tt <- tts[[clust]]
    if (!'t' %in% colnames(tt)) next
    anals[[clust]] <- list(top_tables = list('test-ctrl'=tt))
  }

  if (!length(anals)) return(NULL)

  es <- crossmeta::es_meta(anals)
  es <- es$all$filt
  return(es)
}

run_esmeta_logc <- function(tts) {
  # all genes
  genes <- lapply(tts, function(tt) {genes <- row.names(tt); names(genes) <- tt$ENTREZID; genes})
  genes <- unlist(unname(genes))
  dups <- duplicated(genes)
  genes <- genes[!dups]

  logfc <- data.frame(row.names = genes)

  for (clust in names(tts)) {
    tt <- tts[[clust]]
    if (!'logFC' %in% colnames(tt)) next
    idx <- match(row.names(tt), genes)
    logfc[[clust]] <- NA
    logfc[idx, clust] <- tt$logFC
  }

  # frequency of most consistent direction
  freq <- apply(logfc, 1, function(r) {
    pos <- r > 0
    pos <- pos[!is.na(pos)]
    max(sum(pos), sum(!pos))/length(pos)
  })
  ord <- order(freq, decreasing = TRUE)

  logfc$logFC <- apply(logfc, 1, mean, na.rm = TRUE)
  logfc$ENTREZID <- names(genes)
  logfc$dprime <- logfc$logFC
  logfc[ord, c('ENTREZID', 'logFC', 'dprime')]
}

# p value meta-analysis
run_pmeta <- function(tts) {
  genes <- Reduce(union, lapply(tts, row.names))
  pvals <- data.frame(row.names = genes)

  for (clust in names(tts)) {
    tt <- tts[[clust]]
    if (!'P.Value' %in% colnames(tt)) next
    idx <- match(row.names(tt), genes)
    pvals[[clust]] <- NA
    pvals[idx, clust] <- tt$P.Value
  }

  pvals$p.meta <- apply(pvals, 1, sumz)
  pvals <- pvals[!is.na(pvals$p.meta),, drop = FALSE]
  pvals$fdr <- p.adjust(pvals$p.meta, method = 'BH')

  pvals <- pvals[order(pvals$fdr), c('p.meta', 'fdr')]
  return(pvals)
}


# adapted from metap
sumz <- function (p, weights = NULL, data = NULL, subset = NULL, na.action = na.fail) {
  p <- p[!is.na(p)]
  weights <- rep(1, length(p))
  keep <- (p > 0) & (p < 1)
  invalid <- sum(1L * keep) < 2
  if (invalid) {
    p.meta <- NA

  } else {
    if (sum(1L * keep) != length(p)) {
      warning("Some studies omitted")
      omitw <- weights[!keep]
      if ((sum(1L * omitw) > 0) & !noweights)
        warning("Weights omitted too")
    }
    zp <- (qnorm(p[keep], lower.tail = FALSE) %*% weights[keep])/sqrt(sum(weights[keep]^2))
    p.meta <- stats::pnorm(zp, lower.tail = FALSE)
  }
  return(p.meta)
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

  # has replicates
  if (all(keep %in% colnames(tt))) {
    es <- tt[, keep]
    colnames(es) <- new

  } else {
    es <- tt[, c('ENTREZID', 'logFC')]
    colnames(es)[2] <- 'mean.logFC'
  }
  return(es)
}

write.csv.safe <- function(df, fname, tozip) {
  if (!nrow(df)) return(tozip)

  utils::write.csv(df, fname)
  tozip <- c(tozip, fname)
  return(tozip)
}
