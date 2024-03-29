% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modules-path-utils.R
\name{construct_path_df}
\alias{construct_path_df}
\title{Used by get_path_df to construct the return result}
\usage{
construct_path_df(top_table)
}
\arguments{
\item{top_table}{Filtered result of \code{\link[limma]{toptable}}
to limit number of plotted genes.}
}
\value{
\code{data.frame} used for plotting query genes in Drugs tab
}
\description{
Used by get_path_df to construct the return result
}
\keyword{internal}
