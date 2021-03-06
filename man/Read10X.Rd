% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_scseq.R
\name{Read10X}
\alias{Read10X}
\title{Load in data from 10X}
\usage{
Read10X(data.dir = NULL, gene.column = 2, unique.features = TRUE)
}
\arguments{
\item{data.dir}{Directory containing the matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv
files provided by 10X. A vector or named vector can be given in order to load
several data directories. If a named vector is given, the cell barcode names
will be prefixed with the name.}

\item{gene.column}{Specify which column of genes.tsv or features.tsv to use for gene names; default is 2}

\item{unique.features}{Make feature names unique (default TRUE)}
}
\value{
If features.csv indicates the data has multiple data types, a list
  containing a sparse matrix of the data from each type will be returned.
  Otherwise a sparse matrix containing the expression data will be returned.
}
\description{
Enables easy loading of sparse data matrices provided by 10X genomics.
}
\keyword{internal}
