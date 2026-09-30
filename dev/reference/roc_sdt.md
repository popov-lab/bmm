# Model-implied ROC curve from a fitted SDT model

Computes the receiver operating characteristic (ROC) curve implied by
the posterior distribution of a signal detection model fit with
[`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md). ROC
curves require a response criterion, so they are defined for
[`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)
and
[`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md)
only (not the criterion-free
[`sdt_mafc()`](https://popov-lab.github.io/bmm/dev/reference/sdt_mafc.md)
and
[`sdt_ranking()`](https://popov-lab.github.io/bmm/dev/reference/sdt_ranking.md)).

## Usage

``` r
roc_sdt(
  fit,
  conditions = NULL,
  n_points = 100,
  probs = c(0.025, 0.975),
  criterion_points = NULL,
  ...
)
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

- n_points:

  Integer. Number of equally-spaced points on the smooth model-implied
  ROC curve (default 100); used for both binary and rating models. The
  rating data frame itself still holds K+1 points per draw (K-1
  threshold points plus the (0,0) and (1,1) endpoints).

- probs:

  Numeric vector of length 2. Lower and upper quantiles for the credible
  band (default `c(0.025, 0.975)`).

- criterion_points:

  Optional. Control of the binary multi-criteria behaviour. `NULL`
  (default) auto-detects predictors that vary the criterion only (not
  `d`/`sdratio`) and treats their levels as operating points on one
  curve. Pass a character vector of column names to force that
  classification, or `FALSE` to disable it (one separate curve per
  predictor combination). Ignored for rating models.

- ...:

  Additional arguments passed to
  [`brms::posterior_linpred()`](https://mc-stan.org/rstantools/reference/posterior_linpred.html),
  such as `draw_ids` to use a subset of the posterior draws. `ndraws` is
  refused: each parameter is drawn by its own call, so random subsets
  would not match across parameters.

## Value

A data frame of class `"bmm_sdt_roc"` with columns `FA`, `Hit`, `.draw`,
and any condition columns. The object carries a `summary` attribute
(`FA`, `Hit_mean`, `Hit_lower`, `Hit_upper`) with the smooth
model-implied curve, and a `points` attribute with the model-implied
operating points: one per criterion level for binary multi-criteria
fits, or the K-1 confidence thresholds (labelled `t1`..`t(K-1)`) for
rating fits. It also carries the attributes `probs`, `model_class`,
`dist`, `is_rating` and `conditions`, which the
[`print()`](https://rdrr.io/r/base/print.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods read.
The `Hit_lower` and `Hit_upper` columns in the `summary` attribute are
the pointwise posterior quantiles of the hit rate at each criterion
value, plotted at that criterion's posterior-mean false-alarm rate; this
band does not include uncertainty in the false-alarm rate and is
narrower than a credible band at a fixed false-alarm rate.

## Details

For **binary** models the ROC is traced analytically from the posterior
of `d` (and `sdratio` for unequal-variance SDT) over a grid of criterion
values. When the criterion varies across conditions but `d`/`sdratio` do
not (e.g. a base-rate manipulation, as in
[broeder_schuetz_2009_e3](https://popov-lab.github.io/bmm/dev/reference/broeder_schuetz_2009_e3.md)),
the several criteria are operating points on a single curve: `roc_sdt()`
returns one smooth curve and attaches the model-implied points (one per
criterion level) as a `points` attribute.

For **rating** models the K-1 confidence thresholds define K-1 empirical
ROC points per posterior draw (returned as the data frame). The smooth
model-implied curve is traced over a virtual cut from the posterior of
`d` (and `sdratio`) and attached as the `summary` attribute, with the
K-1 thresholds attached as the `points` attribute (labelled
`t1`..`t(K-1)`) so they fall on the curve.

## See also

[`auc_sdt()`](https://popov-lab.github.io/bmm/dev/reference/auc_sdt.md),
[`roc_observed()`](https://popov-lab.github.io/bmm/dev/reference/roc_observed.md),
[`plot.bmm_sdt_roc()`](https://popov-lab.github.io/bmm/dev/reference/plot.bmm_sdt_roc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Three base-rate conditions shift the criterion, which identifies sdratio
dat <- expand.grid(id = 1:20, stimulus = c(0L, 1L),
                   condition = c("liberal", "neutral", "strict"))
dat$n_trials <- 100L
criteria <- c(liberal = -0.5, neutral = 0, strict = 0.5)
dat$n_old <- rsdt_yn(nrow(dat), dat$n_trials, dat$stimulus, d = 1.5,
                     criterion = criteria[as.character(dat$condition)],
                     sdratio = 1.3)

fit <- bmm(
  formula = bmf(d ~ 1, criterion ~ 0 + condition, sdratio ~ 1),
  data = dat,
  model = sdt_yn(response = "n_old", stimulus = "stimulus",
                 n_trials = "n_trials"),
  cores = 4,
  backend = "cmdstanr"
)

roc <- roc_sdt(fit)
roc
plot(roc, observed = roc_observed(fit))
plot(roc, observed = roc_observed(fit), scale = "quantile")
} # }
```
