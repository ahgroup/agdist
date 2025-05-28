#' Generate Grantham Distance Matrix
#'
#' Computes a 20x20 symmetric matrix of pairwise Grantham distances between
#' amino acids.
#'
#' The matrix uses the amino acid order from the Grantham distance table on
#' Wikipedia, based on Grantham (1984) [DOI:10.1126/science.185.4154.862]. These
#' distances reflect differences in composition, polarity, and molecular volume
#' between amino acid side chains.
#'
#' @details
#' This function is optimized for speed and runs quickly enough that
#' regenerating the matrix on each call is faster than saving/loading from disk.
#'
#' @references
#' - Grantham (1984) [DOI:10.1126/science.185.4154.862]
#'
#' @return A 20x20 numeric matrix with row and column names corresponding to the
#'   one-letter codes for amino acids (lowercase), where each cell contains the
#'   Grantham distance between the corresponding pair of amino acids.
#'
#' @examples
#' grantham_matrix <- generate_grantham_matrix()
#' grantham_matrix["a", "v"]
#'
#' @export
generate_grantham_matrix <- function() {
  # Column order of amino acids from wikipedia table
  aa_order <- c(
    "s", "r", "l", "p", "t", "a", "v", "g", "i", "f", "y",
    "c", "h", "q", "n", "k", "d", "e", "m", "w"
  )

  # Raw Grantham matrix values from wikipedia
  grantham_values <- c(
    110, 145, 74, 58, 99, 124, 56, 142, 155, 144, 112, 89, 68, 46, 121, 65,
    80, 135, 177, 102, 103, 71, 112, 96, 125, 97, 97, 77, 180, 29, 43,
    86, 26, 96, 54, 91, 101, 98, 92, 96, 32, 138, 5, 22, 36, 198, 99, 113,
    153, 107, 172, 138, 15, 61, 38, 27, 68, 42, 95, 114, 110, 169, 77, 76,
    91, 103, 108, 93, 87, 147, 58, 69, 59, 89, 103, 92, 149, 47, 42, 65, 78,
    85, 65, 81, 128, 64, 60, 94, 113, 112, 195, 86, 91, 111, 106, 126, 107,
    84, 148, 109, 29, 50, 55, 192, 84, 96, 133, 97, 152, 121, 21, 88, 135,
    153, 147, 159, 98, 87, 80, 127, 94, 98, 127, 184, 21, 33, 198, 94, 109,
    149, 102, 168, 134, 10, 61, 22, 205, 100, 116, 158, 102, 177, 140, 28,
    40, 194, 83, 99, 143, 85, 160, 122, 36, 37, 174, 154, 139, 202, 154, 170,
    196, 215, 24, 68, 32, 81, 40, 87, 115, 46, 53, 61, 29, 101, 130, 94,
    23, 42, 142, 174, 101, 56, 95, 110, 45, 160, 181, 126, 152, 67
  )

  # Convert to an R matrix and fill in the upper traingle and diagonal
  mat <- matrix(0, nrow = 20, ncol = 20, dimnames = list(aa_order, aa_order))
  mat[lower.tri(mat)] <- grantham_values
  mat <- mat + t(mat)
  return(mat)
}

#' Generate FLU Amino Acid Substitution Matrix
#'
#' Constructs the FLU substitution matrix for amino acid changes, based on
#' data from the FLU model developed by Dang et al. (2010). The FLU model
#' is specifically tuned for influenza protein evolution and is formatted in
#' PAML/PhyML-style: the lower triangle of the rate matrix followed by
#' equilibrium frequencies.
#'
#' @return A 20x20 numeric matrix representing the FLU amino acid substitution rates,
#'         with amino acids ordered alphabetically by their one-letter codes.
#'
#' @details
#' This function hardcodes the substitution matrix values from the original
#' dataset provided by Dang et al. (2010) to avoid repeated file I/O during usage.
#' The final matrix is symmetric and suitable for phylogenetic and protein
#' evolution analyses using the FLU model.
#'
#' @references
#' - Dang, C.C., Le, Q.S., Gascuel, O., & Le, V.S. (2010). FLU, an amino acid substitution
#' model for influenza proteins. *BMC Evolutionary Biology*, 10, 99.
#' \doi{10.1186/1471-2148-10-99}
#' - Source file available from: \url{ftp://ftp.sanger.ac.uk/pub/1000genomes/lsq/FLU}
#'
#' @export
generate_FLU_matrix <- function() {
  input_data <- c(
    # Substitution rates (lower triangle of 20x20 matrix, excluding diagonal)
    0.138658764751059, 0.0533665787145181, 0.161000889039552, 0.584852305649886,
    0.00677184253227681, 7.73739287051356, 0.0264470951166826, 0.16720700818221,
    1.30249856764315e-005, 0.014132062548787, 0.353753981649393, 3.29271694159791,
    0.530642655337477, 0.145469388422239, 0.00254733397966779, 1.4842345032161,
    0.124897616909194, 0.0616521921873234, 5.37051127867923, 3.91106992668137e-011,
    1.19562912226203, 1.13231312248046, 1.19062446519178, 0.322524647863997,
    1.93483278448943, 0.116941459124876, 0.108051341246072, 1.59309882471598,
    0.214757862168721, 1.87956993845887, 1.38709603234116, 0.887570549414031,
    0.0218446166959521, 5.33031341222104, 0.256491863423002, 0.0587745274250666,
    0.149926734229061, 0.246117171830255, 0.21857197541607, 0.0140859174993809,
    0.00111215807314139, 0.0288399502994541, 0.0142107118685268, 1.62662283098296e-005,
    0.243190142026506, 0.0231169515264061, 0.296045557460629, 0.000835873174542931,
    0.00573068208525287, 0.00561362724916376, 1.02036695531654, 0.016499535540562,
    0.00651622937676521, 0.321611693603646, 3.51207228207807, 0.474333610192982,
    15.3000966197798, 2.6468479652886, 0.290042980143818, 3.83228119049152e-006,
    2.559587177122, 3.88148880863814, 0.264148929349066, 0.347302791211758,
    0.227707997165566, 0.129223639195248, 0.0587454231508643, 0.890162345593224,
    0.00525168778853117, 0.0417629637305017, 0.111457310321926, 0.190259181297527,
    0.313974351356074, 0.00150046692269255, 0.00127350890508147, 9.01795420287895,
    6.74693648486614, 1.33129161941264, 0.0804909094320368, 0.0160550314767596,
    0.000836445615590923, 1.0600102849456e-006, 0.10405366623526, 0.0326806570137471,
    0.00100350082518749, 0.00123664495412902, 0.119028506158521, 1.46335727834648,
    2.98680003596399, 0.319895904499071, 0.279910508981581, 0.659311477863896,
    0.154027179890711, 0.0364417719063219, 0.188539456415654, 1.59312060172652e-013,
    0.712769599068934, 0.319558828428154, 0.0386317614553493, 0.924466914225534,
    0.0805433268150369, 0.634308520867322, 0.195750631825315, 0.0568693216513547,
    0.0071324304661639, 3.01134451903854, 0.950138410087378, 3.88131053061457,
    0.338372183381345, 0.336263344504404, 0.487822498528951, 0.307140298031341,
    1.58564657669139, 0.580704249811294, 0.290381075260226, 0.570766693213698,
    0.283807671568883, 0.00702658828739369, 0.996685669575839, 2.08738534433198,
    5.4182981753166, 0.183076905018197, 2.14033231636063, 0.135481232622983,
    0.011975265782196, 0.60234096342392, 0.2801248951174, 0.0188080299490973,
    0.368713573381758, 2.90405228596936, 0.0449263566753846, 1.52696419998775,
    2.03151132062208, 0.000134906239484254, 0.54225109402693, 2.2068599339404,
    0.195966354027106, 1.36942940801512, 0.000536284040016542, 1.4893873721753e-005,
    0.0941066800969967, 0.0440205200833047, 0.155245492137294, 0.196486447133033,
    0.0223729191088972, 0.0321321499585514, 0.431277662888057, 4.97641445484395e-005,
    0.0704600385245663, 0.814753093809928, 0.000431020702277328, 0.0998357527014247,
    0.207066205546908, 0.0182892882245349, 0.0998554972524385, 0.373101926513925,
    0.525398542949365, 0.601692431136271, 0.0722059354079545, 0.104092870343653,
    0.0748149970972622, 6.44895444648517, 0.273934263183281, 0.340058468374384,
    0.0124162215506117, 0.874272174533394, 5.39392424532822, 0.000182294881489116,
    0.392552239890831, 0.124898020409882, 0.42775543040588, 3.53200526987468,
    0.103964386383736, 0.0102575172450253, 0.297123975243582, 0.0549045639492389,
    0.406697814049488, 0.285047948309311, 0.337229618868315, 0.0986313546653266,
    14.3940521944257, 0.890598579382591, 0.0731279296372675, 4.90484223478739,
    0.592587985458668, 0.0589719751511691, 0.0882564232979724, 0.654109108255219,
    0.256900461407996, 0.167581646770807,

    # Equilibrium frequencies (last 20 values)
    0.0470718, 0.0509102, 0.0742143, 0.0478596, 0.0250216,
    0.0333036, 0.0545874, 0.0763734, 0.0199642, 0.0671336,
    0.0714981, 0.0567845, 0.0181507, 0.0304961, 0.0506561,
    0.0884091, 0.0743386, 0.0185237, 0.0314741, 0.0632292
  )

  # Separate rates and equilibrium frequencies
  equilibrium_frequencies <- tail(input_data, 20)
  matrix_data <- head(input_data, length(input_data) - 20)

  aa_order <- c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i",
                "l", "k", "m", "f", "p", "s", "t", "w", "y", "v")

  # Fill upper triangle and transpose
  # We fill the upper triangle instead of the lower
  mat <- matrix(0, nrow = 20, ncol = 20, dimnames = list(aa_order, aa_order))
  mat[upper.tri(mat)] <- matrix_data
  final_matrix <- mat + t(mat)

  return(final_matrix)
}

#' Calculate Substitution Distance Between Amino Acid Sequences
#'
#' Computes the average substitution distance between two aligned amino acid
#' sequences based on a specified substitution matrix. Currently supports the
#' Grantham and FLU matrices.
#'
#' @param seq1 A character string representing the first aligned amino acid
#'   sequence.
#' @param seq2 A character string representing the second aligned amino acid
#'   sequence.
#' @param method Character string specifying the substitution matrix to use.
#'   Supported values are \code{"grantham"} and \code{"flu"} (case-insensitive).
#' @param ambiguous_residues A character string of ambiguous residues to remove
#'   before computing distance.
#'
#' @return A numeric scalar representing the mean substitution distance between
#'   \code{seq1} and \code{seq2}.
#'
#' @details This function first removes ambiguous residues from both sequences
#'   using \code{\link{remove_ambiguous_residues}}, validates the remaining
#'   residues, and computes pairwise distances using a substitution matrix. The
#'   result is normalized by sequence length.
#'
#'   Eventually we plan to support more matrices like BLOSUM and Sneath's index.
#'   If you want to use a specific substitution matrix please let us know.
#'
#' @seealso \code{\link{generate_grantham_matrix}},
#'   \code{\link{generate_FLU_matrix}}
#'
#' @references
#' - Grantham, R. (1974). Amino acid difference formula to help
#' explain protein evolution. \emph{Science}, 185(4154), 862-864.
#' \doi{10.1126/science.185.4154.862}
#'
#' - Dang, C.C., Le, Q.S., Gascuel, O., & Lartillot, N. (2010). FLU, an amino acid
#' substitution model for influenza proteins. \emph{BMC Evolutionary Biology},
#' 10, 99. \doi{10.1186/1471-2148-10-99}
#'
#' @export
substitution <- function(seq1,
                         seq2,
                         method = "grantham",
                         ambiguous_residues = "xX?-") {
  # Assume the sequences get validated before coming into this function
  # Convert to lowercase and remove ambiguous residues
  seqs <- c(seq1, seq2) |>
    tolower() |>
    remove_ambiguous_residues(ambiguous_residues)

  seq1_clean <- strsplit(seqs[[1]], "")[[1]]
  seq2_clean <- strsplit(seqs[[2]], "")[[1]]

  # Select substitution matrix
  method_sanitized <- tolower(method)
  substitution_matrix <- switch(
    method_sanitized,
    "grantham" = generate_grantham_matrix(),
    "flu" = generate_FLU_matrix(),
    cli::cli_abort(c(
      "Invalid {.arg method} specified: {.val {method}}.",
      "i" = "Currently supported options are {.val Grantham} and {.val FLU}.",
      "i" = "(Options are case-insensitive.)"
    ))
  )

  # Validate amino acid characters
  valid_aa <- rownames(substitution_matrix)
  used_aa <- union(seq1_clean, seq2_clean)
  invalid_aa <- setdiff(used_aa, valid_aa)

  if (length(invalid_aa) > 0) {
    cli::cli_art(paste0(
      "Sequences contain invalid amino acid characters: ",
      paste(invalid_aa, collapse = ", ")
    ))
  }

  # Compute pairwise substitution distances
  distances <- substitution_matrix[cbind(seq1_clean, seq2_clean)]

  # Return mean distance
  return(mean(distances))
}

#' Pairwise Substitution Distances Between Amino Acid Sequences
#'
#' Computes pairwise substitution distances between aligned amino acid sequences
#' in a vector, using a specified substitution matrix. Only the lower triangle
#' of the distance matrix is calculated to reduce redundant computations.
#'
#' @param seqs A named character vector of aligned amino acid sequences. All
#'   sequences must be of equal length and named.
#' @param method Character string specifying the substitution matrix to use.
#'   Supported values are \code{"grantham"} and \code{"flu"} (case-insensitive).
#' @param ambiguous_residues A character string of ambiguous residues to remove
#'   before computing distances.
#'
#' @return A symmetric numeric matrix of pairwise mean substitution distances.
#'
#' @details Only the lower triangle of the matrix is computed to avoid redundant
#' calculations. The diagonal is set to zero. The matrix is then symmetrized
#' before being returned.
#'
#' @seealso \code{\link{generate_grantham_matrix}},
#'   \code{\link{generate_FLU_matrix}}
#'
#' @examples
#' seqs <- c(
#'   "A/H1N1/South Carolina/1/1918" = "mktiialsyifclvlgqdfpgndnstat",
#'   "A/H3N2/Darwin/9/2021" = "mktiialsnilclvfaqkipgndnstat",
#'   "B/Sichuan/379/1999" = "drictgitssnsphvvktatqgevnvtg"
#' )
#' dist_substitution(seqs, method = "grantham")
#'
#' @export
dist_substitution <- function(seqs,
                              method = "grantham",
                              ambiguous_residues = "xX?") {
  # Validate input
  validate_sequence_input_form(
    seqs,
    require_alignment = TRUE,
    require_names = FALSE
  )

  n <- length(seqs)
  res <- matrix(0, nrow = n, ncol = n, dimnames = list(names(seqs), names(seqs)))

  # Compute pairwise distances only for lower triangle and then copy to the
  # upper triangle. Diagonals are pre-filled with zero.
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      res[i, j] <- substitution(
        seqs[[i]], seqs[[j]],
        method = method,
        ambiguous_residues = ambiguous_residues
      )
      res[j, i] <- res[i, j]
    }
  }

  return(res)
}

