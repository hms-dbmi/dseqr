% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_scseq.R
\name{load_cellranger_counts}
\alias{load_cellranger_counts}
\title{Load cell ranger counts}
\usage{
load_cellranger_counts(data_dir)
}
\arguments{
\item{data_dir}{Path to folder with cell ranger files.}
}
\value{
dgCMatrix
}
\description{
Mainly to avoid having to download massive datasets that have already been
quantified.
}
\keyword{internal}
