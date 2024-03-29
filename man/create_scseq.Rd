% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_scseq.R
\name{create_scseq}
\alias{create_scseq}
\title{Load kallisto/bustools quantification into a SingleCellExperiment object.}
\usage{
create_scseq(data_dir, tx2gene_dir, project)
}
\arguments{
\item{data_dir}{Directory with raw and kallisto/bustools or CellRanger quantified single-cell RNA-Seq files.}

\item{tx2gene_dir}{Path to directory containing transcript to gene maps
produced by \link[dseqr.data]{load_tx2gene}.}

\item{project}{String identifying sample.}
}
\value{
\code{SingleCellExperiment} object with empty droplets removed and ambience in rowData.
}
\description{
Load kallisto/bustools quantification into a SingleCellExperiment object.
}
\keyword{internal}
