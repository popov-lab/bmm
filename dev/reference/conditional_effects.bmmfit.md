# Conditional Effects for BMM Models

Compute conditional effects for parameters of a bmmfit object. This
method provides a more intuitive interface than directly calling
[`brms::conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)
on bmmfit objects, by:

- Accepting model parameter names directly (e.g., `"kappa"`, `"thetat"`)

- Automatically determining whether parameters are distributional or
  non-linear

- Optionally applying inverse link transformations to show parameters on
  their natural scale

## Usage

``` r
# S3 method for class 'bmmfit'
conditional_effects(x, par = NULL, scale = c("native", "sampling"), ...)
```

## Arguments

- x:

  A bmmfit object (created by
  [`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md))

- par:

  Character string. Name of the model parameter to compute effects for.
  This should be one of the parameter names from the original model
  specification (see `names(x$bmm$model$parameters)`). If `NULL` (the
  default), conditional effects are computed for all estimated
  (non-fixed) parameters.

- scale:

  Character. Scale on which to show the parameter:

  `"native"` (default)

  :   Show on natural scale using inverse link transformation. For
      example, `kappa` with log link shown on exp scale, `thetat` with
      logit (mixture2p) or softmax (mixture3p) link shown on the
      probability scale.

  `"sampling"`

  :   Show on the sampling scale (as used during MCMC). For example,
      `kappa` with log link shown on log scale.

- ...:

  Additional arguments passed to
  [`brms::conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html).
  Common arguments include:

  - `effects`: Character vector specifying which predictor effects to
    plot

  - `conditions`: Named list for setting values of covariates

  - `int_conditions`: Conditions for interactions

  - `prob`: Probability mass to include in credible intervals (default
    0.95)

  - `spaghetti`: Logical, whether to add spaghetti lines

  - `method`: Method for computing effects ("posterior_predict" or
    "posterior_epred"; "posterior_predict" is not available for models
    with a multinomial family, see Details)

## Value

A `brms_conditional_effects` object (from brms), which can be:

- Plotted directly using
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

- Converted to a data frame for custom plotting

- Combined with other conditional effects plots

## Details

### Parameter Types

bmm models use two types of parameters internally:

- **Non-linear parameters (`nlpar`)**: Model parameters that enter the
  likelihood through a non-linear formula, like `kappa` and `thetat` of
  [`mixture3p()`](https://popov-lab.github.io/bmm/dev/reference/mixture3p.md)
  or `c` and `a` of
  [`m3()`](https://popov-lab.github.io/bmm/dev/reference/m3.md)

- **Distributional parameters (`dpar`)**: Parameters of the response
  distribution itself, like `kappa` and `c` of
  [`sdm()`](https://popov-lab.github.io/bmm/dev/reference/sdm.md)

Users should not need to know this distinction -
`conditional_effects.bmmfit()` automatically routes to the correct
parameter type.

### Scale Transformations

By default (`scale = "native"`), parameters are shown on their natural
scale by applying inverse link transformations:

- `log` link → exp transformation

- `logit` link → inverse logit (probability scale)

- `tan_half` link → 2\*atan transformation (radians)

- `identity` link → no transformation

Use `scale = "sampling"` to see parameters on the scale used during MCMC
sampling.

`estimate__`, `lower__`, `upper__` and the spaghetti lines follow
`scale`. `se__` does not: it is the spread of the draws on the scale
`brms` or bmm summarises them on, so for a link other than identity it
can sit on a different scale than `estimate__`. For example, `se__` of
`kappa` is on the sampling scale in
[`mixture3p()`](https://popov-lab.github.io/bmm/dev/reference/mixture3p.md)
and on the native scale in
[`sdm()`](https://popov-lab.github.io/bmm/dev/reference/sdm.md), at
either `scale`. Use `lower__` and `upper__` to describe uncertainty on
the requested scale. With `robust = FALSE`, `estimate__` is the mean and
`se__` the standard deviation of the draws on that scale, and
`estimate__` is then mapped to the requested scale: at
`scale = "native"` the estimate of a non-linear parameter with a log
link is the exponential of the mean log value, not the mean of the
exponentiated draws.

The default of `robust` is `TRUE` (median and MAD), as in
[`brms::conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html),
with one exception: `thetat` and `thetant` of
[`mixture3p()`](https://popov-lab.github.io/bmm/dev/reference/mixture3p.md)
at `scale = "native"` default to `robust = FALSE` (mean and standard
deviation).

### Models with a multinomial family

In models fitted with a multinomial response
([`m3()`](https://popov-lab.github.io/bmm/dev/reference/m3.md),
[`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md),
[`sdt_cdp()`](https://popov-lab.github.io/bmm/dev/reference/sdt_cdp.md),
[`sdt_ranking()`](https://popov-lab.github.io/bmm/dev/reference/sdt_ranking.md)),
the parameters are latent quantities: they have no predictive
distribution, no response categories and no response points. For these
models:

- The effect is the posterior of the parameter itself, summarised by the
  median and MAD of the draws; use `robust = FALSE` for their mean and
  standard deviation.

- The grid of the plotted predictor and the values of all other
  predictors are set by the rules of
  [`brms::conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html).

- Without `effects`, only the effects `brms` lists for the parameter's
  own formula are returned. A parameter without predictors has no
  effects.

- `categorical`, `ordinal`, `select_points` and `transform` are errors,
  and `method` must be `"posterior_epred"` or `"posterior_linpred"` (or
  an alias
  [`brms::conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)
  reads as one of them, such as `"fitted"`), which agree for a
  parameter.

## See also

[`brms::conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)
for the underlying brms function

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit a mixture model with set size effect on kappa
fit <- bmm(
  formula = bmf(kappa ~ 0 + setsize, thetat ~ 1),
  data = zhang_luck_2008,
  model = mixture3p(
    resp_error = "response_error",
    nt_features = paste0("col_lure", 1:5),
    set_size = "setsize"
  )
)

# Get conditional effects for kappa on natural scale (exp of log)
ce_kappa <- conditional_effects(fit, par = "kappa")
plot(ce_kappa)

# Get conditional effects for kappa on log scale (sampling scale)
ce_kappa_log <- conditional_effects(fit, par = "kappa", scale = "sampling")
plot(ce_kappa_log)

# Get effects for thetat (memory probability)
ce_thetat <- conditional_effects(fit, par = "thetat")
plot(ce_thetat)

# Specify which effects to plot
ce_specific <- conditional_effects(fit, par = "kappa", effects = "setsize")

# Combine with other brms options
ce_detailed <- conditional_effects(
  fit,
  par = "kappa",
  effects = "setsize",
  spaghetti = TRUE,
  ndraws = 100
)
} # }
```
