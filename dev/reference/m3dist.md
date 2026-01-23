# Distribution functions for the Memory Measurement Model (M3)

Density and random generation functions for the memory measurement
model. Please note that these functions are currently not vectorized.

## Usage

``` r
dm3(
  x,
  pars,
  m3_model,
  act_funs = construct_m3_act_funs(m3_model, warnings = FALSE),
  log = TRUE,
  ...
)

rm3(
  n,
  size,
  pars,
  m3_model,
  act_funs = construct_m3_act_funs(m3_model, warnings = FALSE),
  ...
)
```

## Arguments

- x:

  Integer vector of length `K` where K is the number of response
  categories and each value is the number of observed responses per
  category

- pars:

  A named vector of parameters of the memory measurement model

- m3_model:

  A `bmmodel` object specifying the m3 model that densities or random
  samples should be generated for

- act_funs:

  A `bmmformula` object specifying the activation functions for the
  different response categories for the "custom" version of the M3. The
  default will attempt to construct the standard activation functions
  for the "ss" and "cs" model version. For a custom m3 model you need to
  specify the act_funs argument manually

- log:

  Logical; if `TRUE` (default), values are returned on the log scale.

- ...:

  can be used to pass additional variables that are used in the
  activation functions, but not parameters of the model

- n:

  Integer. Number of observations to generate data for

- size:

  The total number of observations in all categories

## Value

`dm3` gives the density of the memory measurement model, and `rm3` gives
the random generation function for the memory measurement model.

## References

Oberauer, K., & Lewandowsky, S. (2019). Simple measurement models for
complex working-memory tasks. Psychological Review, 126(6), 880–932.
https://doi.org/10.1037/rev0000159

## Examples

``` r
model <- m3(
  resp_cats = c("corr", "other", "npl"),
  num_options = c(1, 4, 5),
  choice_rule = "simple",
  version = "ss"
)
dm3(x = c(20, 10, 10), pars = c(a = 1, b = 1, c = 2), m3_model = model)
#> [1] -10.93765
```
