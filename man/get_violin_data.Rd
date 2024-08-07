% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_scseq.R
\name{get_violin_data}
\alias{get_violin_data}
\title{Get data for single-cell violin plots}
\usage{
get_violin_data(
  feature,
  scseq,
  selected_cluster,
  by.sample = FALSE,
  decreasing = feature \%in\% c("ribo_percent", "log10_sum", "log10_detected"),
  with_all = FALSE,
  h5logs = NULL
)
}
\arguments{
\item{feature}{Feature name to generate violin plot for. Either a row or \code{colData} of \code{scseq}.}

\item{scseq}{\code{SingleCellExperiment}.}

\item{selected_cluster}{Name of the selected cluster.}

\item{by.sample}{if \code{TRUE} plot \code{feature} violin for each \code{scseq$batch}. Default (\code{FALSE})
will plot \code{feature} for each \code{scseq$cluster}.}

\item{decreasing}{if \code{TRUE}, violinlines with smaller mean values of \code{feature} will show up on top.
Used to show features where smaller values indicate potential QC issues.}
}
\value{
list used by \link{VlnPlot}
}
\description{
Get data for single-cell violin plots
}
\keyword{internal}
