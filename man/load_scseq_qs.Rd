% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modules-sc-utils.R
\name{load_scseq_qs}
\alias{load_scseq_qs}
\title{Load SingleCellExperiment from qs file}
\usage{
load_scseq_qs(
  dataset_dir,
  meta = NULL,
  groups = NULL,
  with_logs = FALSE,
  with_counts = FALSE
)
}
\arguments{
\item{dataset_dir}{Path to folder with scseq.qs file}

\item{meta}{data.frame with column \code{group} and \code{row.names} as sample
names corresponding to \code{scseq$batch}. Default (\code{NULL}) loads previous
specification from file.}

\item{groups}{character vector of length two. First value is test group name
and second in control group name.}

\item{with_logs}{should logcounts be loaded? Default is \code{FALSE} to increase
speed and reduce memory usage.}

\item{with_counts}{should counts be loaded? Default is \code{FALSE} to increase
speed and reduce memory usage.}
}
\value{
SingleCellExperiment
}
\description{
Also attaches clusters from last applied leiden resolution and stores resolution.
}
