% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff_utils.R
\name{get_mds}
\alias{get_mds}
\title{Get scalings for MDS plots}
\usage{
get_mds(exprs, adj, group)
}
\arguments{
\item{exprs}{\code{matrix} of expression values.}

\item{adj}{\code{matrix} of expression values with surrogate variables/pairs regressed out.}

\item{group}{Character vector with values \code{'control'} and \code{'test'} indicating group membership.}
}
\value{
List of tibbles with MDS scalings with and without SVA
}
\description{
For interactive MDS plot of expression values with and without surrogate variable analysis.
}
\keyword{internal}
