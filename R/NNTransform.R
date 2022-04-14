# from Azimuth
NNTransform <- function(object, meta.data, neighbor.slot = "query_ref.nn", key = "ori.index") {
  on.exit(expr = gc(verbose = FALSE))
  ind <- Seurat::Indices(object[[neighbor.slot]])
  ori.index <- t(x = sapply(X = seq_len(nrow(ind)), FUN = function(i) {
    return(meta.data[ind[i, ], key])
  }))
  rownames(x = ori.index) <- rownames(x = ind)
  methods::slot(object = object[[neighbor.slot]], name = "nn.idx") <- ori.index
  return(object)
}
