% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modules-sc-utils.R
\name{integrate_saved_scseqs}
\alias{integrate_saved_scseqs}
\title{Integrate previously saved SingleCellExperiments}
\usage{
integrate_saved_scseqs(
  sc_dir,
  integration_name,
  dataset_names = NULL,
  scseqs = NULL,
  integration_type = c("harmony", "fastMNN", "Azimuth", "symphony"),
  exclude_clusters = NULL,
  exclude_cells = NULL,
  subset_metrics = NULL,
  is_include = NULL,
  founder = integration_name,
  pairs = NULL,
  hvgs = NULL,
  ref_name = NULL,
  npcs = 30,
  cluster_alg = "leiden",
  resoln = 1,
  progress = NULL,
  value = 0,
  tx2gene_dir = NULL
)
}
\arguments{
\item{sc_dir}{Directory with saved single-cell datasets.}

\item{integration_name}{Name for new integrated analysis.}

\item{integration_type}{Charactor vector of one or more integration types.}

\item{exclude_clusters}{Charactor vector of clusters for excluding cells. Only included to save to args.}

\item{exclude_cells}{Character vector of cell names to exclude.}

\item{subset_metrics}{Metrics to subset based on.}

\item{is_include}{Boolean - are cells that match \code{subset_metrics} included or excluded?}

\item{progress}{optional Shiny \code{Progress} object.}
}
\value{
TRUE is successful, otherwise FALSE
}
\description{
Performs integration and saves as a new analysis.
Used by \code{explore_scseq_clusters} shiny app.
}
\seealso{
\code{\link{run_fastmnn}} \code{\link{run_harmony}}
}
\keyword{internal}
