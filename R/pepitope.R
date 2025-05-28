#' Retrieve Epitope Residue Sites for Influenza HA Subtypes
#'
#' This function returns epitope residue positions (sites A–E and optionally all
#' sites combined) for a given influenza HA subtype. It supports H1N1, H3N2,
#' and influenza B lineages (Victoria, Yamagata, or pre-split) with options for
#' harmonizing B lineages and adjusting for signal peptide inclusion.
#'
#' @param subtype Character. HA subtype or lineage, e.g., `"H1N1"`, `"A(H3N2)"`,
#'   `"B/Victoria"`, etc.
#' @param sites Character vector. Which epitope sites to return. Must be a subset of
#'   `c("a", "b", "c", "d", "e", "all_epitopes")`.
#'   Defaults to `c("a", "b", "c", "d", "e")`.
#' @param harmonize_b_lineages Logical. If `TRUE`, merges all B lineages into a
#'   unified "B" category. Defaults to `TRUE`.
#' @param seq_includes_signal_peptide Logical. If `TRUE`, epitope sites will be
#'   offset by the length of the HA signal peptide. Defaults to `TRUE`.
#'
#' @return A named list of integer vectors, with names corresponding to the
#'   requested `sites`.
#'
#' @examples
#' get_pepitope_sites("A(H1N1)", sites = c("a", "b"))
#' get_pepitope_sites("B/Yamagata", harmonize_b_lineages = FALSE)
#'
#' @export
get_pepitope_sites <- function(
    subtype,
    sites = c("a", "b", "c", "d", "e"),
    harmonize_b_lineages = TRUE,
    seq_includes_signal_peptide = TRUE
) {
  # Now make sure the subtype is one of the allowed ones.
  allowed_subtypes <- c(
    "A(H1N1)", "A(H3N2)", "B/Yamagata", "B/Victoria",
    "B/Presplit", "B"
  )
  if (!(subtype %in% allowed_subtypes)) {
    cli::cli_abort(
      c(
        "x" = "Invalid {.arg subtype} value: {.val {subtype}}.",
        "!" = "Must be one of: {.val {allowed_subtypes}}."
      )
    )
  }

  # Validate the harmonize_b_lineages argument
  if (!is.logical(harmonize_b_lineages) || length(harmonize_b_lineages) != 1) {
    cli::cli_abort(
      "{.arg harmonize_b_lineages} must be a single TRUE or FALSE value."
    )
  }

  # Validate the seq_includes_signal_peptide argument
  if (!is.logical(seq_includes_signal_peptide) ||
      length(seq_includes_signal_peptide) != 1) {
    cli::cli_abort(
      "{.arg seq_includes_signal_peptide} must be a single TRUE or FALSE value."
    )
  }

  # Validate the allowed_sites argument
  allowed_sites <- c("a", "b", "c", "d", "e", "all_epitopes")
  if (!all(sites %in% allowed_sites)) {
    cli::cli_abort(
      c(
        "x" = "Invalid value{?s} in {.arg sites}.",
        "!" = "Must be one or more of: {.val {allowed_sites}}"
      )
    )
  }

  # If harmonize_b_lineages is TRUE and the subtype is a b subtype, change
  # the subtype to b_unified
  b_subtypes <- allowed_subtypes[startsWith(allowed_subtypes, "B")]
  if (isTRUE(harmonize_b_lineages) && (subtype %in% b_subtypes)) {
    subtype <- "B/Unified"
  }

  # Now based on the subtype, output the correct list of pepitope sites
  epitope_sites <- switch(
    subtype,
    "A(H1N1)" = list(
      a = c(118, 120:122, 126:129, 132:135, 137, 139:143, 146, 147, 149, 165,
            252, 253),
      b = c(124, 125, 152:157, 160, 162, 183:187, 189:191, 193:196),
      c = c(34:38, 40, 41, 43:45, 269:274, 276:278, 283, 288, 292, 295, 297,
            298, 302, 303, 305:310),
      d = c(89, 94:96, 113, 117, 163, 164, 166:174, 176, 179, 198, 200, 202,
            204:216, 222:227, 235, 237, 241, 243:245),
      e = c(47, 48, 50, 51, 53, 54, 56:58, 66, 68:75, 78:80, 82:86, 102,
            257:261, 263, 267)
    ),
    "A(H3N2)" = list(
      a = c(122, 124, 126, 130:133, 135, 137:138, 140, 142:146, 150, 152, 168),
      b = c(128, 129, 155:160, 163, 165, 186:190, 192:194, 196:198),
      c = c(44:48, 50, 51, 53, 54, 273, 275:276, 278:280, 294, 297, 299, 300,
            304:305, 307:312),
      d = c(96, 102, 103, 117, 121, 167, 170:177, 179, 182, 201, 203, 207:209,
            212:219, 226:230, 238, 240, 242, 244, 246:248),
      e = c(57, 59, 62, 63, 67, 75, 78, 80:83, 86:88, 91, 92, 94, 109, 260:262,
            265)
    ),
    "B/Unified" = list(
      a = c(121, 122, 123, 125, 126, 134, 135, 136, 137, 139, 141, 142, 144,
            146:151, 155, 157, 176, 177),
      b = c(127, 129, 133, 160:168, 171:174, 195:209),
      c = c(34:40, 288:294, 308, 309, 314:327),
      d = c(93, 101, 102, 116, 120, 175:190, 211:214, 217:220, 222:230, 232,
            233, 241:246, 253:258),
      e = c(42, 44, 48, 56, 58, 59, 63, 71, 73, 75, 77:80, 83:85, 88, 89, 91,
            108, 272, 273, 275:277, 279, 280)
    ),
    "B/Victoria" = list(
      a =  c(121, 122, 123, 125, 126, 134, 135, 136, 137, 139, 141, 142, 144,
             146, 147, 148, 149, 150, 151, 155, 157, 177),
      b = c(127, 129, 133, 160, 161, 162, 163, 164, 165, 166, 168, 172, 174,
            196, 197, 198, 199, 200, 202, 203, 204, 206, 207, 208, 209),
      c = c(34, 35, 36, 37, 38, 39, 40, 289, 291, 292, 293, 294, 309, 315,
            317, 318, 320, 321, 323, 324, 325, 326, 327),
      d = c(93, 101, 102, 116, 120, 176, 179, 180, 182, 183, 184, 185, 186,
            187, 188, 190, 212, 214, 218, 219, 220, 223, 224, 225, 226,
            227, 228, 229, 230, 233, 242, 243, 244, 245, 246, 254, 255,
            256, 257, 258),
      e = c(42, 44, 48, 56, 58, 59, 63, 71, 73, 75, 77, 78, 79, 80, 83, 84, 85,
            88, 89, 91, 108, 273, 276, 277, 280)
    ),
    "B/Yamagata" = ,
    "B/Presplit" = list(
      a = c(121, 122, 123, 125, 126, 134, 135, 136, 137, 139, 141, 142, 144,
            146, 147, 148, 149, 150, 151, 155, 157, 176),
      b = c(127, 129, 133, 160, 161, 162, 163, 164, 165, 167, 171, 173, 195,
            196, 197, 198, 199, 201, 202, 203, 205, 206, 207, 208),
      c = c(34, 35, 36, 37, 38, 39, 40, 288, 290, 291, 292, 293, 308, 314,
            316, 317, 319, 320, 322, 323, 324, 325, 326),
      d = c(93, 101, 102, 116, 120, 175, 178, 179, 181, 182, 183, 184, 185,
            186, 187, 189, 211, 213, 217, 218, 219, 222, 223, 224,
            225, 226, 227, 228, 229, 232, 241, 242, 243, 244, 245, 253, 254,
            255, 256, 257),
      e = c(42, 44, 48, 56, 58, 59, 63, 71, 73, 75, 77, 78, 79, 80, 83, 84, 85,
            88, 89, 91, 108, 272, 275, 276, 279)
    ),
    {cli::cli_abort(c(
      "x" = "Subtype {.val {subtype}} is not supported.",
      "!" = paste0(
        "Use one of: {.val A(H1N1)}, {.val A(H3N2)}, {.val B/Victoria},",
        "{.val B/Yamagata}, or {.val B/Presplit}."
      )
    ))}
  )

  # Next if the signal peptide is part of the sequences, we need to adjust for
  # that, since it looks like the original numbers in the Gupta and Pan
  # papers leave out the signal peptides.
  if (isTRUE(seq_includes_signal_peptide)) {
    signal_peptide_length <- switch(
      subtype,
      "A(H1N1)" = 17L,
      "A(H3N2)" = 16L,
      "B/Victoria" = ,
      "B/Yamagata" = ,
      "B/Presplit" = ,
      "B/Unified" = ,
      "B" = 14L,
      TRUE ~ NA_integer_
    )

    epitope_sites <- lapply(
      epitope_sites,
      function(site) site + signal_peptide_length
    )
  }

  # Now construct the overall epitope residue list and we're done
  epitope_sites$all_epitopes <- sort(unique(unlist(epitope_sites)))

  return(epitope_sites[sites])
}

#' Compute p-Epitope Distance Between Two Sequences
#'
#' Calculates the p-epitope antigenic distance between two aligned amino acid
#' sequences of influenza HA using epitope residue sites defined for the
#' specified subtype. Several summary modes are supported, including the
#' dominant (maximum), mean, median, and Anderson's average definition.
#'
#' @param seq_1,seq_2 Character strings. Aligned amino acid sequences.
#' @param subtype Character. HA subtype or lineage (e.g., `"H1N1"`, `"H3N2"`,
#'   `"B/Victoria"`).
#' @param mode Character. Summary mode to apply across epitope sites. Options:
#'   `"dominant"` (default), `"mean"`, `"median"`, `"anderson"`, or `epitope` to
#'   return per-epitope distances.
#' @param harmonize_b_lineages Logical. Whether to treat B lineages as unified.
#'   Defaults to `TRUE`.
#' @param ambiguous_residues Optional argument (currently unused).
#'
#' @return A single numeric value (summary distance) or a numeric vector of
#'   normalized distances per epitope site if `mode` is `NULL`.
#'
#' @references
#' - Gupta et al. (2006), PMID: 16460844
#' - Pan et al. (2010), PMID: 21123189
#' - Anderson et al. (2018), PMID: 29433425
#'
pepitope <- function(
    seq_1,
    seq_2,
    subtype,
    mode = "dominant",
    harmonize_b_lineages = TRUE,
    ambiguous_residues = "xX?"
  ) {
  # We don't need to validate the sequences cause they get validated
  # as an entire vector

  # Get the numbers for the residues in each epitope
  p_epi_sites <- get_pepitope_sites(
    subtype,
    sites = c('a', 'b', 'c', 'd', 'e'),
    harmonize_b_lineages = TRUE
  )

  # This part calculates the p_epitope distance between the two seqs
  epi_dists_raw <- purrr::map_dbl(
    p_epi_sites,
    \(current_sites) {
      aa_1 <- extract_string_chars(seq_1, current_sites)
      aa_2 <- extract_string_chars(seq_2, current_sites)
      stringdist::stringdist(aa_1, aa_2, method = "hamming")
    }
  )

  # Normalize the distances by the length of each site
  epi_dists <- epi_dists_raw / lengths(p_epi_sites)

  # Now decide which summary measure to return
  if (mode %in% c("dominant", "max", "maximum")) {
    # Dominant p-epitope = maximum p-epitope
    # See Gupta 2006 (PMID 16460844).
    p_epi <- max(epi_dists)
  } else if (mode == "anderson") {
    # In Anderson 2018 (PMID 29433425) they take the average of all of the
    # epitope differences. They call this "p-all-epitope" and it is equivalent
    # in practice but provided here just in case. We leave out the scaling by
    # 20 but otherwise the formula is reproduced exactly from the paper.
    p_epi <- mean(epi_dists / lengths(p_epi_sites))
  } else if (mode  %in% c("all", "average", "mean")) {
    # THis is the original definition of p-all-epitope from Pan et al (2010),
    # PMID: 21123189. They basically look at the Hamming distance but only
    # considering the sites that are in the epitopes.
    # This is the same as the mean of the epitope distances, but we use
    # the exact formula of the paper for posterity.
    p_epi <- sum(epi_dists) / sum(lengths(p_epi_sites))
  } else if (mode == "median") {
    # Median instead of mean because Andreas likes the median.
    p_epi <- stats::median(epi_dists)
  } else if (is.null(mode) || is.na(mode) || mode == "" || mode == "epitope") {
    # If there's no summary mode specified, return the entire vector of
    # differences.
    p_epi <- epi_dists
  } else {
    cli::cli_warn(c(
      "!" = "Unknown {.arg mode} supplied to {.fun pepitope}.",
      "i" = paste0(
        "Available values are {.val dominant}, {.val anderson}, ",
        "{.val all}, {.val median}, or {.val epitope}."
      ),
      "i" = "Returning all normalized epitope distances."
    ))
    p_epi <- epi_dists
  }

  return(p_epi)
}

#' Compute Pairwise p-Epitope Distance Matrix
#'
#' Calculates the pairwise p-epitope distances between all sequences in a
#' character vector, using defined epitope sites for a given influenza subtype.
#'
#' @param seqs Named character vector of amino acid sequences. Each string must
#'   be of equal length and represent a full-length HA protein sequence.
#' @param subtype Character. Influenza subtype or lineage. Currently allowed
#'   values are `"H1N1"`, `"A(H3N2)"`, `"B/Yamagata"`, `"B/Victoria`, and
#'   `"B/Presplit`. If `harmonize_b_lineages` is `TRUE` you can also specify
#'   just `B`.
#' @param mode Character. How to summarize the epitope-wise distances. Options:
#'   - `"dominant"`/`"max"`: return the maximum epitope distance (Gupta et al. 2006).
#'   - `"anderson"`: average of normalized distances (Anderson et al. 2018).
#'   - `"all"`/`"average"`/`"mean"`: mean over all epitope residues
#'   (Pan et al. 2010).
#'   - `"median"`: median of per-epitope distances.
#'   - `NULL` or empty string: returns the full vector of distances.
#' @param harmonize_b_lineages Logical. If `TRUE`, harmonizes B lineages using a
#'   unified epitope definition. Defaults to `TRUE`.
#' @param ambiguous_residues Character vector. Residue symbols to exclude from
#'   comparison (e.g. `"xX?"`, the default).
#'
#' @return A symmetric numeric matrix of p-epitope distances with sequence names
#'   as row and column names.
#'
#' @details The function computes the lower triangle of a distance matrix using
#'   the `pepitope()` function, then mirrors it to fill the full matrix.
#'   Ambiguous residues are removed from each pairwise comparison before
#'   computing distance.
#'
#' @references
#' - Gupta et al. (2006), PMID: 16460844
#' - Pan et al. (2010), PMID: 21123189
#' - Anderson et al. (2018), PMID: 29433425
#'
#' @seealso [`pepitope()`]
#'
#' @export
dist_pepi <- function(
    seqs,
    subtype,
    mode = "dominant",
    harmonize_b_lineages = TRUE,
    ambiguous_residues = "xX?"
) {
  # Set up an empty matrix to hold results
  n <- length(seqs)
  res <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(names(seqs), names(seqs))
  )

  # Calculate the lower triangle of the distance matrix for all unique
  # combinations of two strains
  for (i in 2:nrow(res)) {
    for (j in 1:(i-1)) {
      cleaned_sequences <- remove_ambiguous_residues(
        c(seqs[[i]], seqs[[j]]),
        ambiguous_residues
      )
      res[i, j] <- pepitope(
        cleaned_sequences[[1]],
        cleaned_sequences[[2]],
        subtype = subtype,
        mode = mode,
        harmonize_b_lineages = harmonize_b_lineages,
        ambiguous_residues = ambiguous_residues
      )
      res[j, i] <- res[i, j]
    }
  }

  return(res)
}

