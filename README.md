
<!-- README.md is generated from README.Rmd. Please edit that file -->

# agdist <img src="man/figures/logo.png" align="right" alt="" width="120" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/ahagroup/agdist/workflows/R-CMD-check/badge.svg)](https://github.com/ahgroup/agdist/actions)
[![test-coverage](https://github.com/ahgroup/agdist/workflows/test-coverage/badge.svg)](https://github.com/ahagroup/agdist/actions)
<!-- badges: end -->

**This package is still under active development. If you see this
message, assume that several parts of the package don’t quite work
yet.**

## Description

**agdist** allows easy computation of antigenic distances.

## Quick example

These few lines of code produce a list of antigenic distances.

``` r
# install if needed - package is not yet on CRAN
# remotes::install_github('ahgroup/agdist')
library(agdist)
# load (and check) data
sequencedata <- read.csv("sequence-file.csv")
# compute distances
agdistances <- compute_distances(
  sequencedata$aligned_sequences,
  method = "p-epitope"
)
# look at results
agdistances
```

## Getting Started

If you think something like this is of interest to you, hop on over to
the [*Quick Start* tutorial
(vignette)](https://ahgroup.github.io/agdist/articles/quickstart.html).

## Further information

To dig deeper, check out all available vignettes/tutorials for
**agdist**.

The package website (if you are not already on it) is
<https://ahgroup.github.io/agdist/>. You can access all the tutorials
from there, as well as some additional information.

## Citation and Contributors

This R package was developed by [Zane Billings](https://wzbillings.com/)
and [Andreas Handel](https://www.andreashandel.com/). A full list of
contributors and a Bibtex entry for the citation [can be found
here](https://ahgroup.github.io/agdist/authors.html).

This project was/is partially supported by NIH grants XXX.
