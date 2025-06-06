% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{validate_sequence_input_form}
\alias{validate_sequence_input_form}
\title{Validate a sequence vector}
\usage{
validate_sequence_input_form(
  seqs,
  require_alignment = TRUE,
  require_names = FALSE
)
}
\arguments{
\item{seqs}{A (named) character vector of sequences. Each value should be an
amino acid HA sequence, and each name should be the name of the strain or
isolate for that sequence.}

\item{require_alignment}{Logical. If TRUE, checks that all sequences are the
same length.}

\item{require_names}{Logical. If TRUE, checks that all sequences are named.}
}
\value{
`seqs` (invisibly). Errors early if the validation check fails.
}
\description{
The input vector must be a character vector. Additional checks can be imposed
with the require_alignment and require_names arguments.
}
