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
  [`bmm()`](https://venpopov.com/bmm/reference/bmm.md))

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
    "posterior_epred")

## Value

A `brms_conditional_effects` object (from brms), which can be:

- Plotted directly using
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

- Converted to a data frame for custom plotting

- Combined with other conditional effects plots

## Details

### Parameter Types

bmm models use two types of parameters internally:

- **Non-linear parameters (`nlpar`)**: Core model parameters like
  `kappa`, `c`, `a`, `thetat`

- **Distributional parameters (`dpar`)**: Derived parameters used in
  brms mixture distributions

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
