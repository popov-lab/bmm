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
  population_summary = c("none", "mean", "median"),
  scale = c("native", "sampling"),
  ndraws = NULL,
  draw_ids = NULL,
  ndraws_population = 1000,
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
  With `population_summary`, the group-level effects to integrate over.
  See Details.

- population_summary:

  `"none"` (the default) returns the parameters at the group levels in
  the grid. `"mean"` returns, for every draw and grid cell, the mean and
  SD of the parameter across the levels of the grouping factors, e.g.
  across participants; `"median"` returns their median and quartiles.
  See *Population summaries*.

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

- ndraws_population:

  Number of group levels simulated per posterior draw for the population
  summaries that are computed by Monte Carlo. See *Population
  summaries*.

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
  such as `allow_new_levels` and `sample_new_levels`. These two cannot
  be combined with `population_summary`.

## Value

If `summary = FALSE`, a `data.frame` with one row per draw, grid cell
and parameter, with columns `.chain`, `.iteration`, `.draw`, the grid
variables, `parameter` and `value`. The draw indices refer to the draws
of `x`, so they can be joined with the output of
`posterior::as_draws_df(x)`. If `summary = TRUE`, a `data.frame` with
one row per grid cell and parameter, with columns `Estimate`,
`Est.Error` and the interval bounds.

With `population_summary`, both formats have a `statistic` column after
`parameter`, and one row per statistic. The attribute
`"population_method"` is a `data.frame` with one row per parameter and
statistic: the `link`, the `method` (`"closed form"`, `"Gauss-Hermite"`,
`"root search"`, `"Monte Carlo"` or `"no group-level effects"`), `n`
(the quadrature nodes or the simulated levels per draw) and the `groups`
integrated over.

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
and adds the grouping variables to the grid. For a group-level summary,
three quantities are easily confused under a non-identity link:

1.  **The median subject.** `re_formula = NA` sets all group-level
    effects to zero. Because the inverse link is monotone, this is the
    parameter of the *median* subject, not the mean across subjects.

2.  **The population mean.** `population_summary = "mean"` averages the
    parameter over the distribution of the group-level effects,
    separately for every posterior draw. Its credible interval describes
    the uncertainty about the population mean, and the `sd` statistic is
    the spread across subjects on the native scale.

3.  **A new subject.** Predicting for a new level,

        native_parameters(fit, newdata = transform(nd, id = "new"),
                          allow_new_levels = TRUE, sample_new_levels = "gaussian")

    draws one new subject per posterior draw. Its interval mixes the
    posterior uncertainty with the variability between subjects: it
    predicts where one more subject would fall. It is *not* a credible
    interval for the population mean.

The three differ substantially at between-subject standard deviations
typical for working memory data. `population_summary = "median"` returns
(1) together with the quartiles across subjects.

## Population summaries

With `population_summary = "mean"` or `"median"`, each parameter is
integrated, for every posterior draw and grid cell, over the
distribution of the group-level effects that `re_formula` selects: all
of them for `NULL`, the terms it lists for a formula (with the same
meaning as in
[`brms::posterior_linpred()`](https://mc-stan.org/rstantools/reference/posterior_linpred.html)),
and none for `NA`, which returns the median subject with zero spread.
Grouping variables never span the grid, and values given for them in
`newdata` are ignored.

`"mean"` returns the statistics `mean` and `sd`, `"median"` the
statistics `median`, `q25` and `q75`. They are computed per posterior
draw and summarised only afterwards, so `summary = TRUE` reports, for
example, the posterior mean and credible interval of the population
mean. Two summaries are therefore involved: `population_summary`
summarises *across subjects* within a draw, `robust` *across draws*. A
difference between two conditions is the difference of their population
means, taken draw by draw.

**Computation.** A parameter is computed exactly or by quadrature when
its `brms` formula is linear, its group-level terms are plain terms such
as `(1 + x | id)`, `(1 + x || id)`, `(1 | id:condition)` or terms of
several grouping factors, and its native value is the inverse link of
its own linear predictor. Its linear predictor is then normal across
subjects, with mean `eta` (the median subject) and standard deviation
`s = sqrt(z' Sigma z)`, where `z` is the group-level design of the grid
cell and `Sigma` the covariance matrix of the group-level effects in
that draw, summed over independent terms. Then:

- identity link: `eta` and `s`;

- log link: `exp(eta + s^2 / 2)` and that value times
  `sqrt(exp(s^2) - 1)`;

- probit link: the mean is `pnorm(eta / sqrt(1 + s^2))`;

- all other links, and the SD under the probit link: Gauss-Hermite
  quadrature with 100 nodes, which for the logit and probit links is
  accurate to better than `1e-6` up to `s = 3`;

- the median and quartiles are exact for every link: the inverse link of
  `eta` and of `eta -/+ 0.674 * s`.

All other parameters are computed by Monte Carlo. These include mixture
weights transformed through a softmax, parameters defined by non-linear
formulas, parameters whose transformation depends on the design, and
terms such as `mm()`, `cs()` or `gr(by = )`. For every draw, `brms`
simulates `ndraws_population` new levels of each grouping factor from
that draw's covariance matrix. The same simulated levels enter every
parameter and every grid cell, and the parameters are transformed
jointly. The Monte Carlo error of a population mean is about
`sd / sqrt(ndraws_population)` within a draw, and it widens the
posterior SD of the summary. `native_parameters()` warns when it widens
by more than 2% and names the `ndraws_population` needed. This path uses
the random number generator, so set a seed before the call for
reproducible output.

**Circular location parameters** (`tan_half` link): `"mean"` returns the
circular mean direction and the circular SD `sqrt(-2 * log(rho))`, where
`rho` is the mean resultant length, both by Gauss-Hermite quadrature
with 200 nodes, whose error is at most `2e-7` up to `s = 1.5` and grows
to `1.3e-5` at `s = 2`, `3.8e-4` at `s = 3` and `6.1e-3` at `s = 5`. A
between-participant SD above about `1.5` on the `tan(theta / 2)` scale
is rare. `"median"` returns the circular median, the direction with the
smallest expected arc distance to the parameter, without quartiles. Once
enough mass wraps around `pi` the circular median is no longer
`2 * atan(eta)`.

## What is returned

Parameters fixed to a constant are returned at that constant,
transformed to the native scale. Because `bmm` fixes parameters on the
*link* scale, the native value can differ from the value shown by
[`parameter_info()`](https://venpopov.com/bmm/dev/reference/parameter_info.md):
the `ddm` relative starting point `zr` is fixed at `0` under a `logit`
link and is therefore reported as `0.5`, and the `ezdm` and `cswald`
diffusion constant `s` is fixed at `0` under a `log` link and is
reported as `1`.

Trial-level derived quantities, such as the category probabilities of
`m3`, are not returned; use
[`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
for those. The circular mixture models have no expected response, so
their `posterior_epred()` stops; the probabilities of their response
types follow from the returned parameters, as the section on mixture
weights describes.

The default grid contains the observed combinations of the predictors,
not their full crossing, so cells that were never presented do not
appear. A continuous predictor therefore produces one grid cell per
observed value; supply `newdata` for such models.

## Circular location parameters

The circular models sample their location parameter `mu` through a
`tan_half` link, and it is returned in radians in `(-pi, pi)` — a
response bias relative to the target, since the response variable is the
angular error. It is fixed to `0` unless the `bmmformula` predicts it
explicitly, so an all-zero `mu` means the model never estimated one.

The inverse link is `2 * atan()`, and the caveat above applies to it
with particular force: `2 * atan(b_conditionB)` is neither the bias in
condition B nor the difference between conditions. The bias in condition
B is `2 * atan(b_Intercept + b_conditionB)`, and the difference between
the conditions is that value minus `2 * atan(b_Intercept)`, taken draw
by draw.

## Mixture weights

For `mixture2p`, `thetat` is the probability that a response comes from
memory rather than from guessing.

For `mixture3p` the returned `thetat` and `thetant` are probabilities
but they **do not sum to 1**: the likelihood weighs both against the
guessing distribution, whose weight is the fixed reference, so the
remaining mass `1 - thetat - thetant` is the probability of a guess.
`thetant` is the *total* probability of a non-target response, summed
over the lures, not the probability per lure.

At set sizes where no non-target was presented the model switches its
non-target components off, so `thetant` is reported as exactly `0` and
`thetat` as [`plogis()`](https://rdrr.io/r/stats/Logistic.html) of its
linear predictor, rather than as the value the set-size regression
extrapolates to. This correction reads the set size from the prediction
grid: if you supply `newdata` without the set-size variable, every row
takes the set size of the first row of the model data.

In the `"slot"` and `"slot_averaging"` versions of `mixture3p`,
`ptarget` is a probability *given* that a response comes from memory.
The probability of a memory response, `min(1, K / set_size)`, depends on
the set size and is not returned. If the formula gives set size 1 a
coefficient of its own, `ptarget` at set size 1 reflects only its prior,
because a single item cannot be swapped.

With `variable_precision = TRUE`, `tau` is returned on its native scale,
and `kappa` is the concentration of the average trial.

## See also

[`parameter_info()`](https://venpopov.com/bmm/dev/reference/parameter_info.md),
[`native_transform()`](https://venpopov.com/bmm/dev/reference/native_transform.md),
[`conditional_effects.bmmfit()`](https://venpopov.com/bmm/dev/reference/conditional_effects.bmmfit.md)

## Examples

``` r
if (FALSE) { # identical(Sys.getenv("BMM_EXAMPLES"), "true")
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

# with group-level effects: the mean and SD across subjects of kappa
fit_re <- bmm(
  bmf(c ~ 0 + set_size, kappa ~ 1 + (1 | ID)),
  data = oberauer_lin_2017,
  model = sdm(resp_error = "dev_rad")
)
native_parameters(fit_re, pars = "kappa", population_summary = "mean",
                  summary = TRUE)
}
```
