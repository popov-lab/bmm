# Latent decision thresholds from a fitted rating SDT model

Extracts the \\K-1\\ latent decision thresholds of a rating signal
detection model fit with
[`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md).
Thresholds are not posterior parameters: they are reconstructed per draw
from the `criterion` and the parameterization parameters (`spacing`
and/or `delta*`), whose mapping to ordered thresholds depends on the
model's `threshold_type`. This function returns those draws on the
latent decision-variable scale together with a posterior summary, so
threshold estimates are accessible without knowing the parameterization.
Which threshold separates "noise" from "signal" responses depends on
whether the number of categories is even or odd; see the section "Where
`criterion` sits" in
[`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md).

## Usage

``` r
sdt_thresholds(fit, conditions = NULL, probs = c(0.025, 0.975), ...)
```

## Arguments

- fit:

  A `bmmfit` object returned by
  [`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md) from
  an SDT model.

- conditions:

  Optional data frame of predictor values at which to evaluate the
  model. Column names must match predictor variables used in the
  formula. If `NULL` (default), unique predictor combinations are
  derived from the data.

- probs:

  Numeric vector of length 2. Lower and upper quantiles for the credible
  interval (default `c(0.025, 0.975)`).

- ...:

  Additional arguments passed to
  [`brms::posterior_linpred()`](https://mc-stan.org/rstantools/reference/posterior_linpred.html),
  such as `draw_ids` to use a subset of the posterior draws. `ndraws` is
  refused: each parameter is drawn by its own call, so random subsets
  would not match across parameters.

## Value

A data frame of class `"bmm_sdt_thresholds"` with columns `marker`
(threshold label `t1`, `t2`, ...), `position` (latent location),
`.draw`, and any condition columns. The object carries a `summary`
attribute (`marker`, `position` posterior mean, `lower`, `upper`, plus
condition columns). The `position`/`lower`/`upper` naming matches the
`lines` attribute of
[`latent_sdt()`](https://popov-lab.github.io/bmm/dev/reference/latent_sdt.md),
which visualises the same quantities. The object also carries `probs`,
`model_class`, `dist` and `conditions`.

## See also

[`latent_sdt()`](https://popov-lab.github.io/bmm/dev/reference/latent_sdt.md),
[`roc_sdt()`](https://popov-lab.github.io/bmm/dev/reference/roc_sdt.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- expand.grid(id = 1:20, stimulus = c(0L, 1L))
dat <- cbind(dat, rsdt_rating(nrow(dat), 200, dat$stimulus, d = 1.5,
                              thresholds = c(-0.5, 0, 0.5), sdratio = 1.3))

fit <- bmm(
  formula = bmf(d ~ 1, criterion ~ 1, spacing ~ 1, sdratio ~ 1),
  data = dat,
  model = sdt_rating(response = c("r1", "r2", "r3", "r4"),
                     stimulus = "stimulus"),
  cores = 4,
  backend = "cmdstanr"
)

sdt_thresholds(fit)
} # }
```
