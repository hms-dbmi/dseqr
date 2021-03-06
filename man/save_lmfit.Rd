% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff_utils.R
\name{save_lmfit}
\alias{save_lmfit}
\title{Save lmfit result to disk}
\usage{
save_lmfit(lm_fit, dataset_dir, numsv = 0, anal_suffix = "")
}
\arguments{
\item{lm_fit}{Result of run_limma or run_limma_scseq}

\item{dataset_dir}{directory to save results in}

\item{numsv}{Number of surrogate variables modeled. Default is 0.}

\item{anal_suffix}{suffix to append to saved name.}
}
\description{
Save lmfit result to disk
}
\keyword{internal}
