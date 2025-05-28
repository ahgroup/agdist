#' Convert Racmacs Map to Pairwise Distance Matrix
#'
#' This utility function converts the output of a fitted Racmacs map into a
#' pairwise Euclidean distance matrix between all antigenic strains.
#'
#' @param map A fitted Racmacs map object, typically the result of a Racmacs
#'   optimization run. The map must contain 2D coordinates for each antigen
#'   (strain).
#'
#' @return A symmetric matrix of pairwise Euclidean distances between strains,
#'   with row and column names corresponding to strain names.
#'
#' @details This function extracts the 2D coordinates of antigens from the
#'   Racmacs map and computes Euclidean distances between all pairs. Only the
#'   lower triangle is explicitly computed, and the upper triangle is filled by
#'   symmetry.
#'
#' @examples
#' \dontrun{
#'   map <- Racmacs::read.acmap("my_map_file.ace")
#'   distance_matrix <- racmaps_map_to_distances(map)
#'   print(distance_matrix)
#' }
#' @export
racmaps_map_to_distances <- function(map) {
  requireNamespace("Racmacs", quietly = TRUE)

  # Get the coordinates from the map
  coords <- Racmacs::agCoords(map)
  strains <- rownames(coords)

  # Holder for results, sets diagonals to zero
  res <- stats::dist(coords)

  # fill in the upper triangle
  out <- as.matrix(res)

  return(out)
}

#' Extract distances from a phylogenetic tree
#'
#' Computes a distance matrix from a tree object.
#'
#' @param x An object representing a phylogenetic tree or similar structure.
#' @param ... Additional arguments passed to methods.
#'
#' @return A numeric matrix with pairwise distances between tips.
#' @export
tree_to_distances <- function(x, ...) {
  UseMethod("tree_to_distances")
}

#' Extract distances from a phylo object
#'
#' Computes the pairwise cophenetic distance matrix from a phylogenetic tree of
#' class `phylo`. The output is sorted by tip labels.
#'
#' @param x A phylogenetic tree of class `phylo` (from the `ape` package).
#' @param ... Not used.
#'
#' @return An object of class `dist`, with the distances sorted by tip labels.
#' @export
#' @method tree_to_distances phylo
tree_to_distances.phylo <- function(x, ...) {
  requireNamespace("ape", quietly = TRUE)

  tree_dist_matrix <- ape::cophenetic.phylo(x)

  # Sort the distance matrix by row and column names
  sorted_labels <- sort(rownames(tree_dist_matrix))
  tdm_sorted <- tree_dist_matrix[sorted_labels, sorted_labels]

  return(tdm_sorted)
}

#' Extract distances from a pml object
#'
#' Computes the pairwise cophenetic distance matrix from a phylogenetic tree of
#' class `pml` produced by the `phangorn` package.
#' The output is sorted by tip labels.
#'
#' @param x A phylogenetic tree of class `pml` (from the `phangorn` package).
#' @param ... Not used.
#'
#' @return An object of class `dist`, with the distances sorted by tip labels.
#'
#' @export
#'
#' @method tree_to_distances pml
tree_to_distances.pml <- function(x, ...) {
  requireNamespace("ape", quietly = TRUE)
  requireNamespace("phangorn", quietly = TRUE)

  tree_dist_matrix <- ape::cophenetic.phylo(x$tree)

  # Sort the distance matrix by row and column names
  sorted_labels <- sort(rownames(tree_dist_matrix))
  tdm_sorted <- tree_dist_matrix[sorted_labels, sorted_labels]

  return(tdm_sorted)
}

#' Default method for `tree_to_distances`
#'
#' This method is called when `tree_to_distances()` is used on an unsupported
#' object class. It prints an informative message listing available methods.
#'
#' @param x An object of unsupported class.
#' @param ... Not used.
#'
#' @return No return value; prints a message to the console.
#' @export
#' @method tree_to_distances default
tree_to_distances.default <- function(x, ...) {
  available_methods <- utils::methods("tree_to_distances")
  method_classes <- available_methods |>
    as.character() |>
    sub(pattern = "tree_to_distances\\.", replacement = "")

  method_classes <- method_classes[!method_classes == "default"]

  cli::cli_abort(c(
    "!" = "No method available for objects of class {.cls {class(x)}}.",
    "i" = "Methods are available for these classes: {.cls {method_classes}}."
  ))

  msg <- c(
    "!" = "No method available for objects of class {.cls {class(x)}}.",
    "i" = "Methods are available for these classes:",
    stats::setNames(
      sprintf("- {.cls %s}", method_classes),
      rep(" ", length(method_classes))
    )
  )

  cli::cli_abort(msg)

  invisible(x)
}

