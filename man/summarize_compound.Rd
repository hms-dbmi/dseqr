% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modules-drugs-utils.R
\name{summarize_compound}
\alias{summarize_compound}
\title{Summarize query results and annotations by perturbation}
\usage{
summarize_compound(query_table, is_genetic = FALSE, ntop = 1500)
}
\arguments{
\item{query_table}{\code{data.frame} of perturbation correlations and annotations.}

\item{is_genetic}{is \code{query_table} from L1000 genetic perts?}

\item{ntop}{Integer indicating the number of rows to keep after sorting by \code{arrange_by}.}
}
\value{
\code{data.frame} of perturbation correlations and annotations summarized by perturbation.
}
\description{
Takes a \code{data.frame} with one row per signatures and summarizes to one row per compound.
}
\details{
Variables related to individual signatures (cell line, dose, duration, and sample number) are
pasted together and added as a list to \code{'title'} column. Query correlation values are also added as a list to
the \code{'Correlation'} column.

Clinical status is summarized by keeping the most advanced phase only (e.g. Launched > Phase 3). For all other variables,
all unique entries are paste together seperated by \code{'|'}.
}
\keyword{internal}
