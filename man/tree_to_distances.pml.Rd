% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-based.R
\name{tree_to_distances.pml}
\alias{tree_to_distances.pml}
\title{Extract distances from a pml object}
\usage{
\method{tree_to_distances}{pml}(x, ...)
}
\arguments{
\item{x}{A phylogenetic tree of class `pml` (from the `phangorn` package).}

\item{...}{Not used.}
}
\value{
An object of class `dist`, with the distances sorted by tip labels.
}
\description{
Computes the pairwise cophenetic distance matrix from a phylogenetic tree of
class `pml` produced by the `phangorn` package.
The output is sorted by tip labels.
}
