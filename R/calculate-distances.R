calculate_distance <- function(x, method, ...) {
  # Validate the sequences input

  # Dispatch the correct method, error if it's not a match
  method <- tolower(method)
  if (method == "hamming") {
    cli::cli_abort("Sorry, the {.val {method}} method is still WIP!")
  } else if (method == "stringdist") {
    cli::cli_abort("Sorry, the {.val {method}} method is still WIP!")
  } else if (method == "p-epitope") {
    cli::cli_abort("Sorry, the {.val {method}} method is still WIP!")
  } else if (method == "p-all-epitope") {
    cli::cli_abort("Sorry, the {.val {method}} method is still WIP!")
  } else if (method == "cophenetic") {
    cli::cli_abort("Sorry, the {.val {method}} method is still WIP!")
  } else if (method == "cartographic") {
    cli::cli_abort("Sorry, the {.val {method}} method is still WIP!")
  } else if (method == "temporal") {
    cli::cli_abort("Sorry, the {.val {method}} method is still WIP!")
  } else {
    cli::cli_abort(c(
      "The {.var method} {.val {method}} is not supported.",
      "i" = "See the documentation for a list of allowed {.var method} values."
    ))
  }
}
