% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modules-sc-server.R
\name{handle_sc_progress}
\alias{handle_sc_progress}
\title{Update Progress from Background Processes}
\usage{
handle_sc_progress(bgs, progs, new_dataset)
}
\arguments{
\item{bgs}{\code{reactivevalues} of \link[callr]{r_bg}}

\item{progs}{\code{reactivevalues} of \link[shiny]{Progress}}

\item{new_dataset}{\code{reactive} that triggers update of available datasets.}
}
\description{
Update Progress from Background Processes
}
\keyword{internal}