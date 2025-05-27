#' Validate a sequence vector
#'
#' The input vector must be a character vector. Additional checks can be imposed
#' with the require_alignment and require_names arguments.
#'
#' @param seqs A (named) character vector of sequences. Each value should be an
#'   amino acid HA sequence, and each name should be the name of the strain or
#'   isolate for that sequence.
#' @param require_alignment Logical. If TRUE, checks that all sequences are the
#'   same length.
#' @param require_names Logical. If TRUE, checks that all sequences are named.
#'
#' @returns `seqs` (invisibly). Errors early if the validation check fails.
#' @export
#'
validate_sequence_input_form <- function(seqs,
                                         require_alignment = TRUE,
                                         require_names = FALSE) {
  # Check type
  if (!is.character(seqs)) {
    cli::cli_abort("{.arg seqs} must be a character vector.")
  }

  # Check for names if required
  if (require_names) {
    if (is.null(names(seqs)) || any(names(seqs) == "")) {
      cli::cli_abort(paste0(
        "{.arg seqs} must be a named character vector if",
        "{.arg require_names = TRUE}."
      ))
    }
  }

  # Check that all sequences are non-empty strings
  if (!all(nzchar(seqs))) {
    cli::cli_abort("All sequences in {.arg seqs} must be non-empty strings.")
  }

  # Check alignment if required
  if (require_alignment && length(unique(nchar(seqs))) != 1) {
    cli::cli_abort(paste0(
      "Sequences are not aligned: not all elements have the same",
      "number of characters."
    ))
  }

  invisible(seqs)
}


#' Remove ambiguous residues from aligned sequences
#'
#' This function removes positions containing specified ambiguous residues
#' (e.g., "X") from a character vector of aligned sequences. Removal is
#' performed listwise—i.e., any position containing an ambiguous residue in any
#' sequence will be removed from all sequences. This ensures all sequences
#' remain aligned.
#'
#' @param seqs A character vector of aligned amino acid sequences. Sequences may
#'   optionally be named.
#' @param ambiguous_residues A length-1 character string where each character
#'   represents a residue to treat as ambiguous. Defaults to `"xX?"`.
#'
#' @return A character vector of the same length as `seqs`, with the same names
#'   (if any), where all positions containing ambiguous residues in any sequence
#'   have been removed.
#'
#' @details
#' If `ambiguous_residues` is an empty string `""`, no residues will be removed.
#'
#' A common modification that is required for some distance metrics is
#' adding a gap character to `ambiguous_residies`, i.e., `"xX?-"`.
#'
#' @examples
#' seqs <- c(a = "ACDXFG", b = "AXCXFG", c = "ACDYFG")
#' remove_ambiguous_residues(seqs)
#' # Returns c(a = "ACFG", b = "ACFG", c = "ACFG")
#'
#' @export
remove_ambiguous_residues <- function(seqs, ambiguous_residues = "xX?") {
  # Validate seqs
  validate_sequence_input_form(seqs)

  # Validate ambiguous_residues
  if (is.null(ambiguous_residues) || is.na(ambiguous_residues)) {
    cli::cli_abort("{.arg ambiguous_residues} must not be NULL or NA.")
  }
  if (!is.character(ambiguous_residues) || length(ambiguous_residues) != 1L) {
    cli::cli_abort("{.arg ambiguous_residues} must be a single character string.")
  }
  if (isFALSE(nzchar(ambiguous_residues))) {
    return(seqs)
  }

  # Convert to vector of single characters
  residues <- strsplit(ambiguous_residues, "")[[1]]

  # Split sequences into lists of characters
  seq_split <- strsplit(seqs, "")

  # Find ambiguous residue positions
  ambiguous_idx <- Reduce(
    union,
    lapply(seq_split, function(x) which(x %in% residues)),
    init = integer(0)
  )

  # Remove ambiguous positions
  if (length(ambiguous_idx) > 0) {
    seq_clean <- lapply(seq_split, function(x) x[-ambiguous_idx])
  } else {
    seq_clean <- seq_split
  }

  # Recombine into character vector
  seq_out <- vapply(seq_clean, paste0, collapse = "", FUN.VALUE = character(1))
  names(seq_out) <- names(seqs)

  return(seq_out)
}

