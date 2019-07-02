#' Explore Single Cell Clusters
#'
#' @param scseq \code{SingleCellExperiment} or \code{Seurat} object
#' @param markers Named list of \code{data.frame}s where \code{row.names} are the marker genes. One list per cluster in \code{scseq}
#'  with the same name as the cluster.
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' data_dir <- '~/Documents/Batcave/zaklab/drugseqr/data-raw/single-cell/example-anals/sjia'
#'
#' explore_scseq_clusters(data_dir)
#'

explore_scseq_clusters <- function(data_dir, pt.size = 3) {
  # for development so that auto refresh when change file
  options(shiny.autoreload = TRUE)

  app_dir <- '~/Documents/Batcave/zaklab/drugseqr/inst/shiny-apps/scseq'

  # pass arguments to app through options then run
  shiny::shinyOptions(data_dir = data_dir, pt.size = pt.size)
  shiny::runApp(app_dir, launch.browser = TRUE)

}
