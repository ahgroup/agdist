#' Calculate pairwise distances between sequences using various methods
#'
#' This function provides a unified interface for computing pairwise distances
#' between sequences using a variety of distance metrics. It dispatches to
#' appropriate internal functions based on the specified method.
#'
#' @param x A character vector of sequences (typically amino acid or nucleotide
#'   sequences). Currently, you can pass NULL if you want to use the
#'   `"cartographic"` or `"cophenetic"` methods without specifying a sequence
#'   vector.
#' @param method A character string specifying the distance method to use.
#'   Available options include:
#'   \itemize{
#'     \item `"hamming"`: Computes Hamming distance.
#'     \item `"dl"`: Computes Damerau-Levenshtein distance.
#'     \item `"p-epitope"`: Computes p-epitope distances for a specified subtype, by default the dominant p-peitope distance is returned but you can access all arguments of \code{\link{dist_pepi}}.
#'     \item `"p-all-epitope"`: Computes p-all-epitope distances for a subtype.
#'     \item `"cophenetic"` or `"tree"`: Computes distances based on a phylogenetic tree.
#'     \item `"cartographic"` or `"cartography"`: Computes cartographic distances from a Racmacs map.
#'     \item `"temporal"`: Computes temporal distances based on sequence names.
#'   }
#' @param subtype (Optional) A character string indicating the subtype (e.g.,
#'   `"H1"`, `"B"`) to use when computing p-epitope distances. Required for
#'   methods `"p-epitope"` and `"p-all-epitope"`.
#' @param tree (Optional) A `phylo` object used when `method = "cophenetic"` or
#'   `"tree"`.
#' @param map (Optional) A Racmacs map object used when `method =
#'   "cartographic"` or `"cartography"`.
#' @param ... Additional arguments passed to the corresponding distance
#'   functions.
#'
#' @return A numeric matrix of pairwise distances between the input sequences.
#'
#' @details This function acts as a dispatcher, validating the `method` input
#'   and routing the data to the correct underlying distance function. Required
#'   arguments for specific methods (e.g., `subtype`, `tree`, `map`) must be
#'   supplied when applicable or an error will be thrown.
#'
#' @seealso
#'   \code{\link{dist_string}} for `"hamming"` and `"dl"` distances,
#'   \code{\link{dist_pepi}} for `"p-epitope"` and `"p-all-epitope"` distances,
#'   \code{\link{tree_to_distances}} for `"cophenetic"` or `"tree"` distances,
#'   \code{\link{racmaps_map_to_distances}} for `"cartographic"` or `"cartography"` distances,
#'   \code{\link{dist_temporal}} for `"temporal"` distances.
#'
#' @examples
#' \dontrun{calculate_distance(my_seq, "hamming")}
#' \dontrun{calculate_distance(my_seq, "p-epitope", subtype = "A(H1N1)")}
#' \dontrun{calculate_distance(my_seq, "cartographic", map = my_ac_map)}
#'
#' @export
calculate_distance <- function(x,
                               method,
                               subtype = NULL,
                               tree = NULL,
                               map = NULL,
                               ...) {
  dotslist <- list(...)
  args <- names(dotslist)

  # Assume each method validates the sequences with the alignment and names
  # parameters that it needs, so we only need to check the type of method
  if ((!is.character(method)) && (length(method) == 1)) {
    cli::cli_abort("{.arg method} should be a length-one character vector.")
  }

  # Dispatch the correct method, error if it's not a match
  method <- tolower(method)
  if (method == "hamming") {
    out <- dist_string(x, method = "hamming")
  } else if (method == "dl") {
    out <- dist_string(x, method = "dl")
  } else if (method == "p-epitope") {
    # Validate needed arguments without defaults
    if (!missing(subtype)) {
      out <- dist_pepi(x, subtype = subtype, ...)
    } else {
      cli::cli_abort(paste0(
        "If {.arg method} is {.val {method}}, you must also specify the ",
        "{.arg subtype} argument."
      ))
    }
  } else if (method == "p-all-epitope") {
    # Validate needed arguments without defaults
    if (!missing(subtype)) {
      out <- dist_pepi(x, subtype = subtype, ..., mode = "all")
    } else {
      cli::cli_abort(paste0(
        "If {.arg method} is {.val {method}}, you must also specify the ",
        "{.arg subtype} argument."
      ))
    }
  } else if (method == "cophenetic" | method == "tree") {
    # TODO add a step that validates that the tree matches the sequence
    # vector that we passed
    if (!missing(tree)) {
      out <- tree_to_distances(tree)
    } else {
      cli::cli_abort(paste0(
        "If {.arg method} is {.val {method}}, you must also specify the ",
        "{.arg tree} argument."
      ))
    }
  } else if (method == "cartographic" | method == "cartography") {
    # TODO add a step that validates that the map matches the sequence
    # vector that we passed
    if (!missing(map)) {
      out <- racmaps_map_to_distances(map)
    } else {
      cli::cli_abort(paste0(
        "If {.arg method} is {.val {method}}, you must also specify the ",
        "{.arg map} argument."
      ))
    }
  } else if (method == "temporal") {
    out <- dist_temporal(x, ...)
  } else {
    cli::cli_abort(c(
      "The {.var method} {.val {method}} is not supported.",
      "i" = "See the documentation for a list of allowed {.var method} values."
    ))
  }

  return(out)
}
