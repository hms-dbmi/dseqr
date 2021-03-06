% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modules-bulk-server.R
\name{bulkPage}
\alias{bulkPage}
\title{Logic for Bulk Data Tab}
\usage{
bulkPage(
  input,
  output,
  session,
  data_dir,
  sc_dir,
  bulk_dir,
  indices_dir,
  gs_dir
)
}
\arguments{
\item{input, output, session}{standard shiny module boilerplate}

\item{data_dir}{path to folder with application name}

\item{sc_dir}{sub folder of \code{data_dir} where single-cell data is stored}

\item{bulk_dir}{sub folder of \code{data_dir} where bulk data is stored}

\item{indices_dir}{Path to directory containing \code{kallisto} indices and whitelists.}
}
\value{
list with reactive \code{new_dataset} that is triggered with a new
  bulk dataset is added.
}
\description{
to be called with \link[shiny]{callModule}
}
