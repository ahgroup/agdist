% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/substitution.R
\name{generate_FLU_matrix}
\alias{generate_FLU_matrix}
\title{Generate FLU Amino Acid Substitution Matrix}
\usage{
generate_FLU_matrix()
}
\value{
A 20x20 numeric matrix representing the FLU amino acid substitution rates,
        with amino acids ordered alphabetically by their one-letter codes.
}
\description{
Constructs the FLU substitution matrix for amino acid changes, based on
data from the FLU model developed by Dang et al. (2010). The FLU model
is specifically tuned for influenza protein evolution and is formatted in
PAML/PhyML-style: the lower triangle of the rate matrix followed by
equilibrium frequencies.
}
\details{
This function hardcodes the substitution matrix values from the original
dataset provided by Dang et al. (2010) to avoid repeated file I/O during usage.
The final matrix is symmetric and suitable for phylogenetic and protein
evolution analyses using the FLU model.
}
\references{
- Dang, C.C., Le, Q.S., Gascuel, O., & Le, V.S. (2010). FLU, an amino acid substitution
model for influenza proteins. *BMC Evolutionary Biology*, 10, 99.
\doi{10.1186/1471-2148-10-99}
- Source file available from: \url{ftp://ftp.sanger.ac.uk/pub/1000genomes/lsq/FLU}
}
