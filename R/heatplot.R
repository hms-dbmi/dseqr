# modified from enrichplot

#' Plot heatmap for genes from goana/kegga analysis
#'
#' @param path_res data.frame, result of \code{limma::goana} or \code{limma::kegga} with
#' significant genes for each pathway added by \code{add_path_genes}.
#' @param foldChange Vector of logFC values with names equal to gene names as in \code{path_res$GeneName}
#'
#' @return \code{ggplot} object. If foldChange provided then tiles colored accordingly.
#' @export
#'
heatplot <- function(path_res, foldChange = NULL) {

  geneSets <- extract_geneSets(path_res)
  d <- list2df(geneSets)

  if (!is.null(foldChange)) {
    d$foldChange <- foldChange[as.character(d[,2])]
    p <- ggplot2::ggplot(d, ggplot2::aes_(~Gene, ~categoryID)) +
      ggplot2::geom_tile(ggplot2::aes_(fill = ~foldChange), color = "white") +
      # ggplot2::scale_fill_continuous(low="blue", high="red", name = "fold change")
      ggplot2::scale_fill_gradient2(
        low = "#00FF00",
        mid = "white",
        high = "#FF0000",
        midpoint = 0)

  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes_(~Gene, ~categoryID)) +
      ggplot2::geom_tile(color = 'white')
  }
  p + ggplot2::xlab(NULL) +
    ggplot2::ylab(NULL) +
    ggplot2::theme_linedraw() +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 60, hjust = 1))
}

extract_geneSets <- function(path_res) {
  geneSets <- path_res$GeneNames
  geneSets <- lapply(geneSets, function(x) strsplit(x, '/')[[1]])
  names(geneSets) <- path_res$Term
  return(geneSets)
}

list2df <- function (inputList) {
  ldf <- lapply(1:length(inputList), function(i) {
    data.frame(categoryID = rep(names(inputList[i]), length(inputList[[i]])),
               Gene = inputList[[i]])
  })
  do.call("rbind", ldf)
}
