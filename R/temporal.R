#' Extract 4-digit Years from Strings
#'
#' This function extracts 4-digit years (in the range 1900–2099) from a vector
#' of strings.
#'
#' @param full_name A character vector containing strings with embedded 4-digit
#'   years.
#'
#' @return An integer named vector of years, with names corresponding to the
#'   input strings. If no year is found in a string, the corresponding value
#'   will be `NA`.
#'
#' @examples
#' extract_year(c(
#'   "A/H1N1/South Carolina/1/1918",
#'   "A/H3N2/Darwin/9/2021",
#'   "B/Sichuan/379/1999"
#' ))
#' seqs <- c(
#'   "A/H1N1/South Carolina/1/1918" = "mktiialsyifclvlgqdfpgndnstatlclgh",
#'   "A/H3N2/Darwin/9/2021" = "mktiialsnilclvfaqkipgndnstat",
#'   "B/Sichuan/379/1999" = "drictgitssnsphvvktatqgevnvtgai"
#' )
#' extract_year(names(seqs))
#'
#' @export
extract_year <- function(full_name) {
  # Regular expression for a 4-digit year starting with 19 or 20
  year_pattern <- "(19|20)\\d{2}"

  # Extract matching substrings and convert to integer
  yrs_match <- regexpr(year_pattern, full_name)

  # If not every sequence has a match, we can stop now.
  if (any(yrs_match == -1)) {
    nmx <- which(yrs_match == -1)
    nonmatches <- full_name[nmx[1:min(length(nmx), 5)]]
    cli::cli_abort(paste0(
      "Some sequence names do not contain a valid 4-digit year: ",
      paste0("{.val ", nonmatches, "} (index ", nmx, ")", collapse = "; "),
      ifelse(length(nmx) > 5, "; and more.", ".")
    ))
  }

  # If every sequence has a match, get the matches and convert to numbers
  years <- regmatches(full_name, yrs_match)
  years <- as.integer(ifelse(years == "", NA, years))
  names(years) <- full_name

  return(years)
}

#' Compute Temporal Distance Matrix Between Sequences
#'
#' Computes a pairwise distance matrix based on the year embedded in sequence
#' names. Supported methods are `"absolute"` (symmetric), `"forward"` (x - y),
#' and `"backward"` (y - x).
#'
#' @param seqs A named character vector. The names must contain 4-digit years.
#' @param temp_dir The method used for computing temporal distance. One of
#'   `"absolute"` (default), `"forward"`, or `"backward"`.
#'
#' @return A square numeric matrix of pairwise temporal distances.
#'
#' @details This function extracts years from the names of `seqs` using a
#'   regular expression, then computes pairwise temporal distances using the
#'   specified method.
#'
#'   Which distance is "backwards" or "forwards" is semantic and depends on the
#'   order of your strain names. They are both provided for convenience.
#'
#' @examples
#' seqs <- c(
#'   "A/H1N1/South Carolina/1/1918" = "mktiialsyifclvlgqdfpgndnstatlclgh",
#'   "A/H3N2/Darwin/9/2021" = "mktiialsnilclvfaqkipgndnstat",
#'   "B/Sichuan/379/1999" = "drictgitssnsphvvktatqgevnvtgai"
#' )
# dist_temporal(seqs, "absolute")
#'
#' @export
dist_temporal <- function(seqs, temp_dir = "absolute") {
  # Validate input format
  validate_sequence_input_form(
    seqs,
    require_alignment = FALSE,
    require_names = TRUE
  )

  # Extract years from names
  years <- extract_year(names(seqs))

  # Normalize method input
  method_san <- tolower(as.character(temp_dir))

  # Dispatch function
  method_fun <- switch(
    method_san,
    "absolute" = function(x, y) abs(x - y),
    "forward"  = function(x, y) x - y,
    "backward" = function(x, y) y - x,
    cli::cli_abort(c(
      "Invalid {.arg temp_dir} value: {.val {temp_dir}}.",
      "x" = "Allowed values are: {.val absolute}, {.val forward}, {.val backward}."
    ))
  )

  # Compute distance matrix
  # Use the outer product function to do the operation pairwise, if n is the
  # length of years then the outer product returns an n x n square matrix where
  # the ij entry is f(years[i], years[j]) where f() is the custom function
  # passed below.
  dist <- outer(years, years, method_fun)

  return(dist)
}



