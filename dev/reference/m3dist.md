# Distribution functions for the Memory Measurement Model (M3)

Density and random generation functions for the memory measurement
model. Please note that these functions are currently not vectorized.

## Usage

``` r
dm3(x, pars, m3_model, act_funs = NULL, log = TRUE, ...)

rm3(n, size, pars, m3_model, act_funs = NULL, unpack = FALSE, ...)
```

## Arguments

- x:

  Integer vector of length `K` where K is the number of response
  categories and each value is the number of observed responses per
  category

- pars:

  A named vector of parameters of the memory measurement model. Note:
  The fixed parameter `b` does not need to be provided - it will be
  automatically added from the model specification if missing.

- m3_model:

  A `bmmodel` object specifying the m3 model that densities or random
  samples should be generated for

- act_funs:

  A `bmmformula` object specifying the activation functions for the
  different response categories. This can be either: (1) Just the
  activation formulas (one for each response category), or (2) A full
  bmmformula including both activation formulas and other parameters. If
  a full formula is provided, only the formulas matching response
  categories will be extracted. The default will attempt to construct
  the standard activation functions for the "ss" and "cs" model version.
  For a custom m3 model you need to specify the act_funs argument
  manually.

- log:

  Logical; if `TRUE` (default), values are returned on the log scale.

- ...:

  can be used to pass additional variables that are used in the
  activation functions, but not parameters of the model

- n:

  Integer. Number of observations to generate data for

- size:

  The total number of observations in all categories

- unpack:

  Logical; if `TRUE` and `n = 1`, returns a named vector instead of a
  matrix. This allows automatic unpacking of response categories into
  separate columns when used with
  [`dplyr::reframe()`](https://dplyr.tidyverse.org/reference/reframe.html).
  Default is `FALSE` for backward compatibility.

## Value

`dm3` gives the density of the memory measurement model, and `rm3` gives
the random generation function for the memory measurement model.

## References

Oberauer, K., & Lewandowsky, S. (2019). Simple measurement models for
complex working-memory tasks. Psychological Review, 126(6), 880–932.
https://doi.org/10.1037/rev0000159

## Examples

``` r
# Basic usage - b parameter is added automatically
model <- m3(
  resp_cats = c("corr", "other", "npl"),
  num_options = c(1, 4, 5),
  choice_rule = "simple",
  version = "ss"
)

# No need to provide b parameter
dm3(x = c(20, 10, 10), pars = c(a = 1, c = 2), m3_model = model)
#> [1] -14.88885
rm3(n = 10, size = 100, pars = c(a = 1, c = 2), m3_model = model)
#>       corr other npl
#>  [1,]   36    58   6
#>  [2,]   41    55   4
#>  [3,]   38    58   4
#>  [4,]   44    51   5
#>  [5,]   28    60  12
#>  [6,]   35    59   6
#>  [7,]   39    58   3
#>  [8,]   36    61   3
#>  [9,]   36    56   8
#> [10,]   33    61   6

# Can also use full formula (activation formulas are extracted automatically)
full_formula <- bmf(
  corr ~ b + a + c,
  other ~ b + a, 
  npl ~ b,
  a ~ 1,
  c ~ 1
)
rm3(n = 10, size = 100, pars = c(a = 1, c = 2), 
    m3_model = model, act_funs = full_formula)
#>       corr other npl
#>  [1,]   34    61   5
#>  [2,]   42    55   3
#>  [3,]   35    57   8
#>  [4,]   33    59   8
#>  [5,]   48    49   3
#>  [6,]   37    56   7
#>  [7,]   39    56   5
#>  [8,]   30    64   6
#>  [9,]   49    48   3
#> [10,]   49    43   8
    
if (FALSE) { # \dontrun{
# Use with dplyr::reframe() for automatic unpacking into columns
library(dplyr)
library(tibble)
param_grid <- expand.grid(a = c(0.5, 1, 1.5), c = c(1, 2, 3))

simulated_data <- param_grid %>%
  rowwise() %>%
  reframe(
    a = a,
    c = c,
    # unpack=TRUE returns named vector; wrap in as_tibble_row for auto-unpacking
    as_tibble_row(rm3(n = 1, size = 100, pars = c(a = a, c = c), 
                      m3_model = model, unpack = TRUE))
  )
# Result has columns: a, c, corr, other, npl
} # }
```
