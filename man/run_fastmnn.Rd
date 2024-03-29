% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_scseq.R
\name{run_fastmnn}
\alias{run_fastmnn}
\title{Run fastMNN integration}
\usage{
run_fastmnn(logcounts, subset.row, scseqs)
}
\arguments{
\item{logcounts}{dgCMatrix of logcounts.}

\item{subset.row}{A vector specifying which features to use for correction.}

\item{scseqs}{list of \code{SingleCellExperiment} objects to integrate.}
}
\value{
Integrated \code{SingleCellExperiment} with logcounts assay, corrected reducedDim, and batch annotation.
}
\description{
Run fastMNN integration
}
\keyword{internal}
