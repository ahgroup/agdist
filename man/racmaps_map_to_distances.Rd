% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-based.R
\name{racmaps_map_to_distances}
\alias{racmaps_map_to_distances}
\title{Convert Racmacs Map to Pairwise Distance Matrix}
\usage{
racmaps_map_to_distances(map)
}
\arguments{
\item{map}{A fitted Racmacs map object, typically the result of a Racmacs
optimization run. The map must contain 2D coordinates for each antigen
(strain).}
}
\value{
A symmetric matrix of pairwise Euclidean distances between strains,
  with row and column names corresponding to strain names.
}
\description{
This utility function converts the output of a fitted Racmacs map into a
pairwise Euclidean distance matrix between all antigenic strains.
}
\details{
This function extracts the 2D coordinates of antigens from the
  Racmacs map and computes Euclidean distances between all pairs. Only the
  lower triangle is explicitly computed, and the upper triangle is filled by
  symmetry.
}
\examples{
\dontrun{
  map <- Racmacs::read.acmap("my_map_file.ace")
  distance_matrix <- racmaps_map_to_distances(map)
  print(distance_matrix)
}
}
