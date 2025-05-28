#' Extract Specific Characters from a String by Position
#'
#' This helper function extracts one or more characters from a string using
#' their byte positions. It converts the input string to raw bytes, selects the
#' specified positions, and converts them back to a character.
#'
#' This is often faster and easier than splitting the string, subsetting,
#' and pasting the string back together.
#'
#' @param str A character string from which characters will be extracted.
#' @param pos An integer vector of positions (1-based) indicating which
#'   characters to extract.
#'
#' @return A character string containing the extracted characters.
#'
#' @examples
#' extract_string_chars("hello", c(1, 2))  # Returns "he"
#'
#' @note This function operates at the raw byte level and may not behave as
#'   expected with multibyte or non-ASCII characters. For the purposes of
#'   this package, all validated sequence strings will only contain
#'   ASCII strings and this will work as expected.
#'
#' @export
extract_string_chars <- function(str, pos) {
  sub <- charToRaw(str)[pos]
  return(rawToChar(sub))
}

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

#' Convert a distance matrix to tidy format
#'
#' Transforms a square numeric distance matrix into a long-format tibble with
#' customizable names for the row and column variables. Optionally, returns only
#' unique row-column combinations (i.e., lower triangle).
#'
#' @param d A square numeric matrix (or coercible to one) with row and column
#'   names.
#' @param rows_to A string (length 1) naming the new column for row labels.
#'   Default is `"Strain1"`.
#' @param cols_to A string (length 1) naming the new column for column labels.
#'   Default is `"Strain2"`.
#' @param unique_pairs Logical (length 1). If `TRUE`, return only the unique
#'   row-column combinations from the lower triangle and diagonal of the matrix.
#'   Defaults to `FALSE`.
#'
#' @return A tibble in long format with columns:
#' \describe{
#'   \item{[rows_to]}{Factor, ordered by the row order of `d`}
#'   \item{[cols_to]}{Factor, ordered and reversed from column order of `d`}
#'   \item{d}{Distance value between row and column entries}
#' }
#'
#' @examples
#' mat <- matrix(1:9, nrow = 3)
#' rownames(mat) <- colnames(mat) <- c("A", "B", "C")
#' tidy_dist_mat(mat)
#' tidy_dist_mat(mat, unique_pairs = TRUE)
#'
#' @export
tidy_dist_mat <- function(d,
                          rows_to = "Strain1",
                          cols_to = "Strain2",
                          unique_pairs = FALSE) {

  # Validate option arguments
  if (!is.character(rows_to) || length(rows_to) != 1) {
    cli::cli_abort("{.arg rows_to} must be a character vector of length 1.")
  }
  if (!is.character(cols_to) || length(cols_to) != 1) {
    cli::cli_abort("{.arg cols_to} must be a character vector of length 1.")
  }
  if (!is.logical(unique_pairs) || length(unique_pairs) != 1) {
    cli::cli_abort((
      "{.arg unique_pairs} must be a logical (TRUE or FALSE) of length 1."
    ))
  }

  # Validate d matrix
  if (is.data.frame(d) || tibble::is_tibble(d)) {
    d <- as.matrix(d)
  }
  if (!is.matrix(d) || !is.numeric(d)) {
    cli::cli_abort(
      "{.arg d} must be a numeric matrix (or a data frame coercible to one)."
    )
  }
  if (nrow(d) != ncol(d)) {
    cli::cli_abort(paste0(
      "{.arg d} must be a square matrix: it has {nrow(d)} rows and {ncol(d)} ",
      "columns."
    ))
  }
  if (is.null(rownames(d)) || is.null(colnames(d))) {
    cli::cli_abort("{.arg d} must have both row names and column names.")
  }

  # Convert to tibble
  df <- tibble::as_tibble(d, rownames = rows_to)

  df_long <- tidyr::pivot_longer(
    df,
    cols = -dplyr::all_of(rows_to),
    names_to = cols_to,
    values_to = "d"
  )

  # Order factor levels
  df_long[[rows_to]] <- forcats::fct_inorder(df_long[[rows_to]])
  df_long[[cols_to]] <- forcats::fct_rev(forcats::fct_inorder(df_long[[cols_to]]))

  # Filter for unique combinations if needed
  if (isTRUE(unique_pairs)) {
    # Get the unique row names and column names
    rn <- levels(df_long[[rows_to]])
    cn <- levels(df_long[[cols_to]])
    # This part finds which rows of the tidy dataframe are in the lower triangle
    # of the distance matrix
    keep <- vapply(
      seq_len(nrow(df_long)),
      function(i) {
        # For each row in the df, get the position of the unique row name
        # and unique column name that are in that row, then compare the order.
        # This makes sure we only include the observations that were in the
        # lower triangle of the original distance matrix.
        r_idx <- match(as.character(df_long[[rows_to]][i]), rn)
        c_idx <- match(as.character(df_long[[cols_to]][i]), cn)
        # check lower triangle (including diagonal)
        # Because we reversed the index of the column names earlier we need to
        # undo the reversal here to make sure the comparison is right
        # But basically this checks if the row name comes before the col name
        return(r_idx >= (length(cn) - c_idx + 1))
      },
      # Template for the returned object
      logical(1)
    )
    # Now filter the dataframe so we only have the lower triangular indices
    # Without dropping unique values
    df_long <- df_long[keep, , drop = FALSE]
  }

  return(df_long)
}


