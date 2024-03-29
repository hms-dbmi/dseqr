% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_scseq.R
\name{import_scseq}
\alias{import_scseq}
\title{Import raw single cell fastq or cellranger files for app}
\usage{
import_scseq(
  dataset_name,
  uploaded_data_dir,
  sc_dir,
  tx2gene_dir,
  species = NULL,
  progress = NULL,
  value = 0,
  founder = dataset_name,
  npcs = 30,
  cluster_alg = "leiden",
  resoln = 1,
  ref_name = NULL,
  metrics = c("low_lib_size", "low_n_features", "high_subsets_mito_percent",
    "low_subsets_ribo_percent", "high_doublet_score")
)
}
\arguments{
\item{dataset_name}{Name of dataset}

\item{uploaded_data_dir}{Directory with fastq or cellranger files}

\item{sc_dir}{Single cell directory for app. Will store results in \code{dataset_name} subdirectory}

\item{tx2gene_dir}{Path to directory containing transcript to gene maps
produced by \link[dseqr.data]{load_tx2gene}.}

\item{species}{Name of species. Used for R object imports only (detected for cellranger files).}

\item{progress}{Optional shiny \code{Progress} object. Default will print progress.}

\item{value}{Integer indicating step of pipeline.}

\item{founder}{Name of dataset that \code{dataset_name} originates from.}
}
\description{
Import raw single cell fastq or cellranger files for app
}
\keyword{internal}
