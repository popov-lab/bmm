# Posterior draws of model parameters on the native scale

`bmm` samples all parameters on their link scale (`log` for precision
and boundary parameters, `logit` or `softmax` for mixture weights,
`tan_half` for circular locations). `native_parameters()` returns
posterior draws of those parameters on the scale they are interpreted
and reported on, evaluated over a grid of predictor values.

Because draws are returned rather than summaries, any contrast is
ordinary arithmetic on the draws, with the baseline made explicit by the
grid.

## Usage

``` r
native_parameters(
  x,
  newdata = NULL,
  pars = NULL,
  re_formula = NULL,
  scale = c("native", "sampling"),
  ndraws = NULL,
  draw_ids = NULL,
  summary = FALSE,
  prob = 0.95,
  robust = FALSE,
  ...
)
```

## Arguments

- x:

  A `bmmfit` object returned by
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md).

- newdata:

  Data to evaluate the parameters on. If `NULL` (the default), the
  unique observed combinations of the parameters' predictors are used,
  taken as rows of the model's own data. Columns you do not supply are
  filled from the first row of the model data.

- pars:

  Character vector of parameters to return. If `NULL` (the default), all
  parameters of the model are returned, including those fixed to a
  constant. This filters the output only; parameters transformed jointly
  are always computed together.

- re_formula:

  Which group-level effects to include, as in
  [`brms::posterior_linpred()`](https://mc-stan.org/rstantools/reference/posterior_linpred.html).
  See Details.

- scale:

  Either `"native"` (the default) to apply the model's inverse link
  functions, or `"sampling"` to return the untransformed linear
  predictor.

- ndraws:

  Number of posterior draws to use. If `NULL` (the default), all draws
  are used. The same draws are used for every parameter.

- draw_ids:

  Indices of the draws to use. Overrides `ndraws` and makes the result
  reproducible.

- summary:

  Logical. If `FALSE` (the default), return draws. If `TRUE`, return
  posterior summaries computed *after* the transformation.

- prob:

  Probability mass of the credible interval when `summary = TRUE`.

- robust:

  Logical. If `TRUE`, `summary = TRUE` reports the median and median
  absolute deviation instead of the mean and standard deviation.

- ...:

  Further arguments passed to
  [`brms::posterior_linpred()`](https://mc-stan.org/rstantools/reference/posterior_linpred.html),
  such as `allow_new_levels` and `sample_new_levels`.

## Value

If `summary = FALSE`, a `data.frame` with one row per draw, grid cell
and parameter, with columns `.chain`, `.iteration`, `.draw`, the grid
variables, `parameter` and `value`. The draw indices refer to the draws
of `x`, so they can be joined with the output of
`posterior::as_draws_df(x)`. If `summary = TRUE`, a `data.frame` with
one row per grid cell and parameter, with columns `Estimate`,
`Est.Error` and the interval bounds.

## Transform first, then summarise

Quantiles are preserved by monotone inverse links but means are not, so
the median and the credible interval of a native-scale parameter are
exact while its mean is not the inverse link of the linear predictor's
mean. This function transforms the draws and only then summarises, which
is why `summary = TRUE` is not the same as transforming the output of
[`summary.bmmfit()`](https://venpopov.com/bmm/dev/reference/summary.bmmfit.md).

The same order is what makes `summary = TRUE` correct for the mixture
weights, whose softmax is not an elementwise map at all: summarising the
transformed draws reports quantiles of the weight's own marginal
posterior, whereas transforming the summaries would report the softmax
of three separate quantiles, which is not a quantile of anything.

## Transforming coefficients is not the same thing

The inverse link applies to the **linear predictor** of a grid cell, not
to an individual regression coefficient. For a model with
`kappa ~ condition`, `exp(b_kappa_conditionB)` is a multiplicative
factor, not `kappa` in condition B; `kappa` in condition B is
`exp(b_kappa_Intercept + b_kappa_conditionB)`. Building the grid is
exactly what removes this step, which is why contrasts are taken between
rows of the output rather than read off the coefficients.

## Group-level effects

`re_formula = NULL` (the default) returns subject-specific parameters
and adds the grouping variables to the grid. `re_formula = NA` sets all
group-level effects to zero, which under a non-identity link gives the
*median* subject, not the population mean. The population mean requires
marginalising over the distribution of group-level effects, which is
done by predicting for a new level:

    native_parameters(fit, newdata = transform(nd, id = "new"),
                      allow_new_levels = TRUE, sample_new_levels = "gaussian")

These three quantities differ substantially at between-subject standard
deviations typical for working memory data.

## What is returned

Parameters fixed to a constant are returned at that constant,
transformed to the native scale. Because `bmm` fixes parameters on the
*link* scale, the native value can differ from the value shown by
[`parameters()`](https://venpopov.com/bmm/dev/reference/parameters.md):
the `ddm` relative starting point `zr` is fixed at `0` under a `logit`
link and is therefore reported as `0.5`, and the `ezdm` and `cswald`
diffusion constant `s` is fixed at `0` under a `log` link and is
reported as `1`.

Trial-level derived quantities, such as the response probabilities of
the mixture models or the category probabilities of `m3`, are not
returned; use
[`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
for those.

The default grid contains the observed combinations of the predictors,
not their full crossing, so cells that were never presented do not
appear. A continuous predictor therefore produces one grid cell per
observed value; supply `newdata` for such models.

## Circular location parameters

The circular models sample their location parameter (`mu` for `sdm`,
`mu1` for `mixture2p`, `mixture3p` and `imm`) through a `tan_half` link,
and it is returned in radians in `(-pi, pi)` — a response bias relative
to the target, since the response variable is the angular error. It is
fixed to `0` unless the `bmmformula` predicts it explicitly, so an
all-zero `mu1` means the model never estimated one.

The inverse link is `2 * atan()`, and the caveat above applies to it
with particular force: `2 * atan(b_mu1_conditionB)` is neither the bias
in condition B nor the difference between conditions. The bias in
condition B is `2 * atan(b_mu1_Intercept + b_mu1_conditionB)`, and the
difference between the conditions is that value minus
`2 * atan(b_mu1_Intercept)`, taken draw by draw.

## Mixture weights

For `mixture3p` the returned `thetat` and `thetant` are probabilities
but they **do not sum to 1**: `brms` holds the linear predictor of one
mixture component at zero, and that component is the guessing
distribution, so the remaining mass `1 - thetat - thetant` is the
probability of a guess. `thetant` is the *total* probability of a
non-target response, summed over the lures, not the probability per
lure.

At set sizes where no non-target was presented the model switches its
non-target components off, so `thetant` is reported as exactly `0` and
`thetat` as [`plogis()`](https://rdrr.io/r/stats/Logistic.html) of its
linear predictor, rather than as the value the set-size regression
extrapolates to. This correction reads the lure indicators out of the
prediction grid; if you supply `newdata`, columns you do not supply are
filled from the first row of the model data rather than recomputed, and
`native_parameters()` warns when the two disagree.

## See also

[`parameters()`](https://venpopov.com/bmm/dev/reference/parameters.md),
[`native_transform()`](https://venpopov.com/bmm/dev/reference/native_transform.md),
[`conditional_effects.bmmfit()`](https://venpopov.com/bmm/dev/reference/conditional_effects.bmmfit.md)

## Examples

``` r
if (FALSE) { # isTRUE(Sys.getenv("BMM_EXAMPLES"))
fit <- bmm(
  bmf(c ~ 0 + set_size, kappa ~ 1),
  data = oberauer_lin_2017,
  model = sdm(resp_error = "dev_rad")
)

# draws of c and kappa for every set size, on the native scale
np <- native_parameters(fit, re_formula = NA)
head(np)

# a contrast is plain arithmetic on the draws
c_draws <- subset(np, parameter == "c")
quantile(
  c_draws$value[c_draws$set_size == 1] - c_draws$value[c_draws$set_size == 4],
  probs = c(0.025, 0.5, 0.975)
)

# posterior summaries of the transformed draws
native_parameters(fit, re_formula = NA, summary = TRUE)
}
```
