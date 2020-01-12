dtangle.genes <- function(top_table, scseq, clusters) {

  scseq <- scseq[, scseq$cluster %in% clusters]
  tests <- pairwise_wilcox(scseq)
  markers <- get_scseq_markers(tests, pval.type = 'all')

  top_table <- top_table[top_table$adj.P.Val < 0.05, ]

  universe <- intersect(row.names(top_table), row.names(markers[[1]]))
  top_table <- top_table[universe, ]

  get_pvals <- function(gene) sapply(markers, function(x) x[gene, 'FDR'])
  pvals <- sapply(row.names(top_table), get_pvals)
  pvals <- t(pvals)

  pheatmap::pheatmap(pvals, show_rownames = FALSE)

}

levels(scseq$cluster) <- annot

clusters <- grep('Mono', levels(scseq$cluster), value = T)

