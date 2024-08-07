% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modules-path-utils.R
\name{get_path_df}
\alias{get_path_df}
\title{Get data.frame for plotting gene expression values of a pathway}
\usage{
get_path_df(top_table, path_id = NULL, pert_signature = NULL, nmax = 200)
}
\arguments{
\item{path_id}{String with KEGG pathway id.}

\item{nmax}{Maximum number of genes to keep from CMAP02/L1000 common and CMAP02 only genes for Drug and genetic query genes. Default is 200
so that all drug and genetic query genes are shown.}
}
\value{
\code{data.frame} with columns: \itemize{
 \item Gene gene names.
 \item Dprime standardized unbiased effect size values.
 \item sd standard deviations of \code{Dprime}.
 \item Link url to GeneCards page for gene.
}
}
\description{
Get data.frame for plotting gene expression values of a pathway
}
\keyword{internal}
