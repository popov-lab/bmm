# Category log-probability under the Gaussian choice rule of `m3()`

R companion to the Stan function `m3_gauss_logp` that
`m3(choice_rule = "gaussian")` places in each category's activation
formula. `brms` evaluates the non-linear formula in R for `log_lik()`,
`posterior_predict()` and `posterior_epred()`, looking the function up
on the search path; it is exported for that reason and is not meant to
be called directly.

## Usage

``` r
m3_gauss_logp(k, ...)
```

## Arguments

- k:

  Integer index of the response category.

- ...:

  The K category activations followed by the K option counts, as numbers
  or draws-by-observation matrices (as supplied by brms).

## Value

The log probability of category `k`, with the shape of the first
activation.
