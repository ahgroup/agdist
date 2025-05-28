# Zane's notes

* Data importing
  - right now I think this should be the user's job. All of our functions should
  accept a character vector of aligned sequences as the argument.
  - Importing and aligning the sequences shouldn't be a job for this package,
  there are already multiple specialized packages for that.
* Simple distances
  - right now all of the sequences accept a character vector and they compute
  distances pairwise between each entry. I think ideally each function should
  have a `cor()` style `x` and `y` argument where the default for `y` is `NULL`,
  and if you specify both `x` and `y` the function computes the distances with
  elementwise comparisons (i.e. `dist(x[1], y[1]), dist(x[2], y[2]), ...`)
  instead of only having the pairwise option.
  - string distances: basically just an interface to `stringdist` that does
  some preprocessing, it's written so that anything in the `...` gets passed
  to the stringdist function. Right now they get automatically length-normalized
  but it should probably be an argument for whether they get normalized.
  - specialized sequence distances
     * I didn't write the p-epitopes the crazy way Andreas wanted to, so if he
     wants users to be able to supply their own sites he needs to add
     that functionality. It doesn't even really make sense to do
     p-epitope on B strains either.
     * dominant pepitope
     * p-all-epitope
     * grantham
     * FLU sub
  - temporal distance: will extract the years from a named sequence vector.
    * also probably needs to have a method for strain names, and a method
  for a vector of years.
    * It should probably try to autodetect, but also
  allow users to manually specify which. Or we can expose the extract_years()
  function to users and only allow actual numbers into the temporal distance
  calculation function, and require users to do it as a processing step.
    * but the main method for `compute_distances()` will look at the names.
* Model-based distances
  - these functions need
  - trees: this package shouldn't build trees for users, there's already a lot
  of packages to do that. unfortunately there's not as much standardization
  for trees, so we probably just have to write multiple methods to support
  multiple different tree formats. we can start by supporting the R package
  phangorn and provide instructions for a decent default tree, but building the
  tree isn't really our job and we aren't experts.
  - cartography: our package shouldn't build cartographies for users. We should
  tell them to use Racmacs. Instead, our distance function will accept a
  Racmacs output file (typically .acmap) or racmacs object as the input.
* Wrapper: `compute_distances()` should take a named sequence vector `x` as the
first input. The names should correspond to strain names, and the values should
be the ALIGNED protein sequences. This will also have to have a `...` argument
to pass specific other details to the other distance functions. It will only
be able to compute one distance at a time unless we can come up with a different
API that works for all distances without overloading the `...`.
