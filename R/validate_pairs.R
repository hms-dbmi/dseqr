#' Validate sample pairing for pair-ended RNA seq
#'
#' @param pairs Numeric vector of integers and/or \code{NA}. Positions with the same integer
#'   value indicate samples that are paired in a pair-ended experiment.
#' @param rows Numeric vector of integers indicating selected rows.
#' @param reps Numeric vector of integers and/or \code{NA}. Positions with the same integer
#'   value indicate samples that replicates.
#'
#' @return \code{TRUE} if the pairing is valid, otherwise \code{FALSE}.
#' @export
#' @keywords internal
#'
validate_pairs <- function(pairs, rows, reps) {

  rep_rows <- reps[rows]
  all_rep_rows <- all(rep_rows %in% rows, na.rm = TRUE)

  if (length(rows) < 2) {
    msg <- "Select at least two rows to mark as pairs."

  } else if (!anyNA(rep_rows) && length(unique(rep_rows)) < 2) {
    msg <- "Select at least two non-replicate rows to mark as pairs."

  } else if (!anyNA(rep_rows) && !all_rep_rows) {
    msg <- "All replicates must be included in the same pair."

  } else if (length(unique(rep_rows)) > 2) {
    msg <- "A pair must include exactly two replicate groups or one replicate group and one additional sample."

  } else if (!all(is.na(pairs[rows]))) {
    msg <- "Selected row(s) already belong to a pair. Click 'Reset' if you need to start over."

  } else {
    msg <- NULL
  }

  return(msg)
}

#' Validate sample replicates for RNA seq
#'
#' @inheritParams validate_pairs
#'
#' @return \code{TRUE} if the replicate is valid, otherwise \code{FALSE}.
#' @export
#' @keywords internal
validate_reps <- function(pairs, rows, reps) {

  if (length(rows) < 2) {
    msg <- "Select at least two rows to mark as replicates."

  } else if (!all(is.na(reps[rows]))) {
    msg <- "Selected row(s) already belong to a replicate. Click 'Reset' if you need to start over."

  } else if (!all(is.na(pairs[rows]))) {
    msg <- "Replicates must be specified first for paired samples that include replicates. Click 'Reset' if you need to start over."

  } else {
    msg <- NULL
  }

  return(msg)
}


#' Get first sequence identifiers for fastq.gz files
#'
#' Function is used to determine if an experiment is single or paired.
#'
#' @param fastq_paths Character vector of paths to fastq.gz files
#'
#' @return Named character vector of first sequence id lines (start with @) in \code{fastq_paths}.
#'   Names are the \code{fastq_paths}.
#' @keywords internal
#' @export
#'
#' @examples
#' # required fastq.gz files in data-raw/example-data (e.g. IBD example)
#' fastq_paths <- list.files(file.path('data-raw', 'example-data'), '.fastq.gz$', full.names = TRUE)
#' fastq_id1s <- get_fastq_id1s(fastq_paths)
#'
get_fastq_id1s <- function(fastq_paths) {

  # get first line with @ symbol (sequence identifier)
  fastq_id1s <- sapply(fastq_paths, function(f) {
    incon <- gzfile(f)
    while (TRUE) {
      line = readLines(incon, n = 1)
      if (grepl('^@', line)) {
        break
      }
    }
    close(incon)
    return(line)
  })
  return(fastq_id1s)
}


#' Detect if experiment is pair-ended.
#'
#' @param fastq_id1s Character vector of first sequence identifiers from fastq.gz files. Returned from \code{\link{get_fastq_id1s}}.
#'
#' @return boolean indicating if experiement is pair-ended (\code{TRUE}) or single-ended (\code{FALSE}).
#' @export
#' @keywords internal
detect_paired <- function(fastq_id1s) {

  # older illumina sequence identifiers have 1 part
  # newer illumina sequence identifiers have 2 space-seperated parts
  id_parts <- strsplit(fastq_id1s, ' ')
  older <- all(sapply(id_parts, length) == 1)
  newer <- all(sapply(id_parts, length) == 2)

  if (older) {
    # pair is 1 or 2 at end of sequence id after /
    pairs <- gsub('^.+?/([12])$', '\\1', fastq_id1s)

  } else if (newer) {
    # pair is 1 or 2 followed by : followed by N or Y at beginning of second part
    id_parts2 <- sapply(id_parts, `[`, 2)
    pairs <- gsub('^([12]):[YN]:.+$', '\\1', id_parts2)

  } else {
    stop("fastq.gz files don't appear to be from older/newer Illumina software. Please contact package author.")
  }

  # SRA also accepts /1 and /2 at end of read name
  is_sra <- any(grepl('^@SRR\\d+', fastq_id1s))
  if (is_sra) {
    pairs <- gsub('^.+?/([12]$)', '\\1', pairs)
  }

  # paired experiments will have '1' and '2'
  uniq.pairs <- unique(pairs)
  paired <- setequal(c('1', '2'), uniq.pairs)

  # unpaired will have only '1'
  if (!paired) stopifnot(uniq.pairs == '1')

  return(paired)
}

