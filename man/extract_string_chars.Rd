% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{extract_string_chars}
\alias{extract_string_chars}
\title{Extract Specific Characters from a String by Position}
\usage{
extract_string_chars(str, pos)
}
\arguments{
\item{str}{A character string from which characters will be extracted.}

\item{pos}{An integer vector of positions (1-based) indicating which
characters to extract.}
}
\value{
A character string containing the extracted characters.
}
\description{
This helper function extracts one or more characters from a string using
their byte positions. It converts the input string to raw bytes, selects the
specified positions, and converts them back to a character.
}
\details{
This is often faster and easier than splitting the string, subsetting,
and pasting the string back together.
}
\note{
This function operates at the raw byte level and may not behave as
  expected with multibyte or non-ASCII characters. For the purposes of
  this package, all validated sequence strings will only contain
  ASCII strings and this will work as expected.
}
\examples{
extract_string_chars("hello", c(1, 2))  # Returns "he"

}
