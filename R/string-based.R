#' Compute Pairwise Normalized String Distances Between Aligned Sequences
#'
#' This function computes a pairwise distance matrix between aligned sequences
#' using a string distance metric from the `stringdist` package. Ambiguous
#' residues (e.g., "X", "x", "?") are removed before distance computation.
#' By default, the hamming distance is compouted, but the
#' \code{stringdist::stringdist()} documentation lists other possible methods
#' and arguments.
#'
#' @param seqs A character vector or list of aligned sequences.
#' @param ambiguous_residues A character vector of ambiguous residues to remove
#'   from each sequence before comparison. Defaults to \code{c("x", "X", "?")}.
#' @param ... Additional arguments passed to \code{stringdist::stringdist()}
#'   (e.g., \code{method = "hamming"}).
#'
#' @return A symmetric numeric matrix of pairwise normalized string distances.
#'   Each distance is normalized by the length of the cleaned first sequence in
#'   the pair.
#'
#' @details Only the lower triangle is explicitly computed, and the upper
#'   triangle is filled in by symmetry. This function assumes that all sequences
#'   are of equal length and aligned.
#'
#'   All distances are normalized by dividing by the aligned sequence length.
#'
#' @examples
#' seqs <- c(
#'   "A/H1N1/South Carolina/1/1918" = "mktiialsyifclvlgqdfpgndnstat",
#'   "A/H3N2/Darwin/9/2021" = "mktiialsnilclvfaqkipgndnstat",
#'   "B/Sichuan/379/1999" = "drictgitssnsphvvktatqgevnvtg"
#' )
#' dist_string(seqs, method = "hamming")
#'
#' @export
dist_string <- function(seqs, ambiguous_residues = "xX?", ...) {
  # Validate sequence input
  validate_sequence_input_form(
    seqs,
    require_alignment = TRUE,
    require_names = FALSE
  )

  n <- length(seqs)
  res <- matrix(
    0,
    nrow = n, ncol = n,
    dimnames = list(names(seqs), names(seqs))
  )

  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      # Remove ambiguous residues
      cleaned_seqs <- remove_ambiguous_residues(
        c(seqs[[i]], seqs[[j]]),
        ambiguous_residues
      )

      # Compute normalized distance
      dist <- stringdist::stringdist(
        a = cleaned_seqs[[1]],
        b = cleaned_seqs[[2]],
        ...
      ) / nchar(cleaned_seqs[[1]])

      res[i, j] <- dist
      res[j, i] <- dist  # symmetric assignment
    }
  }

  return(res)
}

