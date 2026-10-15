# Distribution functions for Multinomial Processing Tree (MPT) models

Density and random generation functions for multinomial processing tree
models specified with
[`mpt()`](https://venpopov.com/bmm/reference/mpt.md). Please note that
these functions are currently not vectorized.

## Usage

``` r
dmpt(x, pars, mpt_model, tree = NULL, log = TRUE, covariates = NULL, ...)

rmpt(
  n,
  size,
  pars,
  mpt_model,
  tree = NULL,
  unpack = FALSE,
  covariates = NULL,
  ...
)
```

## Arguments

- x:

  Integer vector of length `K`, where `K` is the number of response
  categories, giving the number of observed responses per category. If
  `x` is a named vector, its names are matched to the response
  categories and the vector is reordered accordingly; an unnamed vector
  is read in the order of the model's response categories,
  `mpt_model$resp_vars$resp_cats`, which is the order of the first tree
  given to [`mpt()`](https://venpopov.com/bmm/reference/mpt.md)
  ([`mpt()`](https://venpopov.com/bmm/reference/mpt.md) stores the
  branches of every tree in that order, whatever order they were written
  in). `rmpt()` returns its columns in the same order.

- pars:

  A named vector or list with the values of at least the latent
  parameters appearing in the branch expressions of the selected tree,
  on the probability scale and strictly between 0 and 1, as for fixed
  parameters in [`bmm()`](https://venpopov.com/bmm/reference/bmm.md).
  Values for parameters the tree does not use are ignored, so one vector
  can serve every tree of the model. Parameters of a simplex group must
  sum to 1.

- mpt_model:

  A `bmmodel` object created with
  [`mpt()`](https://venpopov.com/bmm/reference/mpt.md) specifying the
  model that densities or random samples should be generated for.

- tree:

  Character. For models with multiple trees, the name of the tree to
  compute probabilities for. Can be omitted for single-tree models.

- log:

  Logical; if `TRUE` (default), densities are returned on the log scale.

- covariates:

  A named vector or list with the values of all covariates declared in
  the model that appear in the branch expressions of the selected tree.

- ...:

  can be used to pass additional values used in the branch expressions,
  as an alternative to the `covariates` argument.

- n:

  Integer. Number of observations to generate data for.

- size:

  The total number of observations across all response categories.

- unpack:

  Logical; if `TRUE` and `n = 1`, returns a named vector instead of a
  matrix. This allows automatic unpacking of response categories into
  separate columns when used with
  [`dplyr::reframe()`](https://dplyr.tidyverse.org/reference/reframe.html).
  Default is `FALSE`.

## Value

`dmpt` gives the multinomial density of a vector of response counts, and
`rmpt` generates random response counts for the response categories.

## Note

Unlike the densities of the circular models in this package
([`dsdm()`](https://venpopov.com/bmm/reference/SDMdist.md),
[`dmixture2p()`](https://venpopov.com/bmm/reference/mixture2p_dist.md),
[`dmixture3p()`](https://venpopov.com/bmm/reference/mixture3p_dist.md),
[`dimm()`](https://venpopov.com/bmm/reference/IMMdist.md)) and unlike
[`stats::dmultinom()`](https://rdrr.io/r/stats/Multinom.html), `dmpt()`
returns the log density by default, matching
[`dm3()`](https://venpopov.com/bmm/reference/m3dist.md).

## References

Batchelder, W. H., & Riefer, D. M. (1999). Theoretical and empirical
review of multinomial process tree modeling. Psychonomic Bulletin &
Review, 6(1), 57-86. https://doi.org/10.3758/BF03210812

## Examples

``` r
tree_target <- mpt_tree("target", list(
  old = "D + (1 - D) * g",
  new = "(1 - D) * (1 - g)"
))
tree_lure <- mpt_tree("lure", list(
  old = "(1 - D) * g",
  new = "D + (1 - D) * (1 - g)"
))
model <- mpt(list(tree_target, tree_lure), tree_id = "item_type")

dmpt(
  x = c(35, 15), pars = c(D = 0.7, g = 0.5),
  mpt_model = model, tree = "target"
)
#> [1] -5.702642
rmpt(
  n = 10, size = 50, pars = c(D = 0.7, g = 0.5),
  mpt_model = model, tree = "lure"
)
#>       old new
#>  [1,]   6  44
#>  [2,]  13  37
#>  [3,]   3  47
#>  [4,]  10  40
#>  [5,]   6  44
#>  [6,]   9  41
#>  [7,]   6  44
#>  [8,]   4  46
#>  [9,]   9  41
#> [10,]   4  46
```
